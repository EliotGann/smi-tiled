"""Benchmark: Compare CPU (_SplitBinPlan) vs GPU (TorchBinPlan) histogram.

Run on a machine with PyTorch + GPU (CUDA or MPS):
    pip install torch
    python scripts/benchmark_gpu_histogram.py

This generates synthetic data mimicking a real SAXS detector frame and
compares correctness and throughput of the two backends.
"""

import time
import numpy as np


def make_synthetic_geometry(ny=1679, nx=1475):
    """Create synthetic q/chi maps mimicking pil2M SAXS geometry."""
    # Beam center near middle
    bc_row, bc_col = ny * 0.55, nx * 0.5
    rows = np.arange(ny, dtype=float)[:, None] - bc_row
    cols = np.arange(nx, dtype=float)[None, :] - bc_col
    # q proportional to distance from beam (simplified)
    pixel_size_m = 0.172e-3
    dist_m = 5.0
    wavelength_nm = 0.077  # 16 keV
    r_m = np.sqrt((rows * pixel_size_m) ** 2 + (cols * pixel_size_m) ** 2)
    two_theta = np.arctan2(r_m, dist_m)
    q2d = (4.0 * np.pi / wavelength_nm) * np.sin(two_theta / 2.0)
    chi2d = np.degrees(np.arctan2(rows * pixel_size_m, cols * pixel_size_m))
    return q2d, chi2d


def make_synthetic_frame(shape, rng):
    """Generate a synthetic detector frame with realistic features."""
    ny, nx = shape
    img = rng.exponential(scale=10.0, size=shape).astype(np.float64)
    # Add a few hot pixels
    hot_idx = rng.integers(0, ny * nx, size=20)
    img.ravel()[hot_idx] = 1e5
    # Mask: ~5% invalid
    valid = rng.random(shape) > 0.05
    return img, valid


def main():
    print("=" * 70)
    print("GPU Histogram Benchmark")
    print("=" * 70)

    # Import backends
    from smi_tiled.integrator import _SplitBinPlan, _bin_indices
    from smi_tiled.gpu_histogram import TorchBinPlan, get_torch_device, is_available

    device = get_torch_device()
    print(f"\n  PyTorch device: {device}")
    print(f"  GPU available:  {is_available()}")

    if device == "cpu":
        print("\n  WARNING: No GPU detected. Running on CPU (for correctness check).")
        print("  For real speedups, run on a machine with CUDA or Apple Metal.\n")

    # Setup
    ny, nx = 1679, 1475
    n_q, n_chi = 500, 180
    print(f"  Detector: {ny}×{nx} = {ny*nx:,} pixels")
    print(f"  Output grid: {n_q}×{n_chi} = {n_q*n_chi:,} bins")

    rng = np.random.default_rng(42)
    q2d, chi2d = make_synthetic_geometry(ny, nx)

    # Bin edges
    q_valid = q2d[np.isfinite(q2d)]
    chi_valid = chi2d[np.isfinite(chi2d)]
    q_edges = np.linspace(q_valid.min(), q_valid.max(), n_q + 1)
    chi_edges = np.linspace(chi_valid.min(), chi_valid.max(), n_chi + 1)

    # Build plans
    print("\n  Building CPU plan...")
    t0 = time.perf_counter()
    cpu_plan = _SplitBinPlan(q2d, chi2d, q_edges, chi_edges, pixel_splitting=1)
    t_cpu_build = time.perf_counter() - t0
    print(f"    CPU plan build: {t_cpu_build:.3f}s")

    print("  Building GPU plan...")
    t0 = time.perf_counter()
    gpu_plan = TorchBinPlan(q2d, chi2d, q_edges, chi_edges, pixel_splitting=1, device=device)
    t_gpu_build = time.perf_counter() - t0
    print(f"    GPU plan build: {t_gpu_build:.3f}s")

    # Correctness check
    print("\n  Correctness check (single frame)...")
    img, valid = make_synthetic_frame((ny, nx), rng)
    I_cpu, N_cpu = cpu_plan.integrate_frame(img, valid)
    I_gpu, N_gpu = gpu_plan.integrate_frame(img, valid)

    I_err = np.max(np.abs(I_cpu - I_gpu))
    N_err = np.max(np.abs(N_cpu - N_gpu))
    print(f"    Max |I_cpu - I_gpu|: {I_err:.2e}")
    print(f"    Max |N_cpu - N_gpu|: {N_err:.2e}")
    # float32 on MPS gives ~1e-4 relative error on large sums; float64 gives ~1e-10
    rtol = 1e-3 if device == "mps" else 1e-10
    I_max = np.max(np.abs(I_cpu)) or 1.0
    N_max = np.max(np.abs(N_cpu)) or 1.0
    assert I_err / I_max < rtol, f"Intensity mismatch: rel={I_err/I_max:.2e}"
    assert N_err / N_max < rtol, f"Counts mismatch: rel={N_err/N_max:.2e}"
    print("    \u2713 Results match (within precision)")

    # Throughput: single-frame
    n_warmup = 3
    n_bench = 20
    print(f"\n  Throughput benchmark ({n_bench} frames, single-frame mode)...")

    frames = [(make_synthetic_frame((ny, nx), rng)) for _ in range(n_warmup + n_bench)]

    # Warmup GPU
    for img, valid in frames[:n_warmup]:
        gpu_plan.integrate_frame(img, valid)

    # Benchmark CPU
    t0 = time.perf_counter()
    for img, valid in frames[n_warmup:]:
        cpu_plan.integrate_frame(img, valid)
    t_cpu = time.perf_counter() - t0

    # Benchmark GPU
    t0 = time.perf_counter()
    for img, valid in frames[n_warmup:]:
        gpu_plan.integrate_frame(img, valid)
    t_gpu = time.perf_counter() - t0

    fps_cpu = n_bench / t_cpu
    fps_gpu = n_bench / t_gpu
    print(f"    CPU: {t_cpu:.3f}s ({fps_cpu:.1f} frames/s, {t_cpu/n_bench*1000:.1f} ms/frame)")
    print(f"    GPU: {t_gpu:.3f}s ({fps_gpu:.1f} frames/s, {t_gpu/n_bench*1000:.1f} ms/frame)")
    print(f"    Speedup: {fps_gpu/fps_cpu:.1f}×")

    # Throughput: batch mode (if GPU)
    if device != "cpu":
        print(f"\n  Batch benchmark ({n_bench} frames, one kernel launch)...")
        batch_imgs = np.stack([f[0] for f in frames[n_warmup:]])
        batch_vals = np.stack([f[1] for f in frames[n_warmup:]])

        # Warmup
        gpu_plan.integrate_batch(batch_imgs[:3], batch_vals[:3])

        t0 = time.perf_counter()
        I_batch, N_batch = gpu_plan.integrate_batch(batch_imgs, batch_vals)
        t_batch = time.perf_counter() - t0

        fps_batch = n_bench / t_batch
        print(f"    Batch GPU: {t_batch:.3f}s ({fps_batch:.1f} frames/s, "
              f"{t_batch/n_bench*1000:.1f} ms/frame)")
        print(f"    Batch speedup vs CPU: {fps_batch/fps_cpu:.1f}×")

        # Verify batch matches single-frame
        I_single, N_single = gpu_plan.integrate_frame(frames[n_warmup][0], frames[n_warmup][1])
        batch_err = np.max(np.abs(I_batch[0] - I_single))
        assert batch_err < 1e-6, f"Batch vs single mismatch: {batch_err}"
        print(f"    ✓ Batch results match single-frame")

    # Projection for full scan
    print(f"\n{'=' * 70}")
    print("  PROJECTION FOR 3721-FRAME SCAN")
    print(f"{'=' * 70}")
    n_scan = 3721
    est_cpu = n_scan / fps_cpu
    est_gpu = n_scan / fps_gpu
    print(f"    CPU estimate: {est_cpu:.1f}s")
    print(f"    GPU estimate: {est_gpu:.1f}s")
    print(f"    Projected savings: {est_cpu - est_gpu:.1f}s")
    if device != "cpu":
        est_batch = n_scan / fps_batch
        print(f"    GPU batch estimate: {est_batch:.1f}s")

    print("\n=== DONE ===")


if __name__ == "__main__":
    main()
