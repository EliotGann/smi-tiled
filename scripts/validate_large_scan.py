"""One-time validation: run reduce_smi_combined on a large scan and report
progress granularity, timing, output shape, and peak memory consumption."""

import sys
import time
import tracemalloc
from smi_tiled.integrator import reduce_smi_combined

UID = "ce94c000-369d-444a-8078-a9ed3c36b872"
TILED_URI = "https://tiled.nsls2.bnl.gov"
CATALOG = "smi/migration"

# Track progress calls
_progress_log: list[tuple[float, str, int, int]] = []
_t0 = time.perf_counter()
_bar_width = 40


def _progress_cb(stage: str, current: int, total: int) -> None:
    elapsed = time.perf_counter() - _t0
    _progress_log.append((elapsed, stage, current, total))
    # tqdm-style progress bar on every callback
    frac = current / total if total > 0 else 0
    filled = int(_bar_width * frac)
    bar = "█" * filled + "░" * (_bar_width - filled)
    pct = frac * 100
    # Estimate remaining time
    eta = ""
    if frac > 0.01:
        remaining = elapsed / frac * (1 - frac)
        eta = f" ETA {remaining:.0f}s"
    sys.stderr.write(
        f"\r  {stage:20s} |{bar}| {current:>5d}/{total} "
        f"[{elapsed:.0f}s{eta}]  "
    )
    sys.stderr.flush()
    # Print a newline when a stage completes (first call or stage change)
    if len(_progress_log) >= 2 and _progress_log[-2][1] != stage:
        sys.stderr.write("\n")
        sys.stderr.flush()


def main():
    global _t0
    tracemalloc.start()
    _t0 = time.perf_counter()

    print(f"=== Validating reduce_smi_combined on {UID} ===\n")

    result = reduce_smi_combined(
        uid=UID,
        tiled_uri=TILED_URI,
        catalog=CATALOG,
        image_cache_path="auto",
        dezinger_threshold=None,  # disable dezinger for this run
        progress=_progress_cb,
    )

    elapsed_total = time.perf_counter() - _t0
    peak_mem_mb = tracemalloc.get_traced_memory()[1] / 1024**2
    tracemalloc.stop()

    sys.stderr.write("\n")  # newline after progress bar

    # --- Report ---
    print("\n" + "=" * 60)
    print("PROGRESS REPORT")
    print("=" * 60)
    n_calls = len(_progress_log)
    print(f"  Total progress callbacks: {n_calls}")
    if n_calls > 0:
        final = _progress_log[-1]
        print(f"  Final state: {final[1]} {final[2]}/{final[3]}")
        # Compute intervals between calls
        intervals = [_progress_log[i][0] - _progress_log[i-1][0]
                     for i in range(1, n_calls)]
        if intervals:
            print(f"  Callback interval — min: {min(intervals):.3f}s, "
                  f"max: {max(intervals):.3f}s, "
                  f"mean: {sum(intervals)/len(intervals):.3f}s")
        # Stage breakdown
        stages: dict[str, int] = {}
        for _, stage, _, _ in _progress_log:
            stages[stage] = stages.get(stage, 0) + 1
        print(f"  Per-stage call counts: {stages}")

    print(f"\n{'=' * 60}")
    print("TIMING")
    print("=" * 60)
    timing = getattr(result, "timing", None) or {}
    if hasattr(result, "__dict__"):
        timing = getattr(result, "timing", None) or {}
    if isinstance(timing, dict):
        for k, v in timing.items():
            print(f"  {k:25s}: {v:.2f}s")
    print(f"  {'TOTAL (wall)':25s}: {elapsed_total:.2f}s")

    print(f"\n{'=' * 60}")
    print("MEMORY")
    print("=" * 60)
    print(f"  Peak traced memory: {peak_mem_mb:.1f} MB")

    print(f"\n{'=' * 60}")
    print("OUTPUT SHAPES")
    print("=" * 60)
    # Access result attributes (may be dataclass/namedtuple or dict)
    _get = getattr(result, "get", None) or (lambda k: getattr(result, k, None))
    for key in ["q_chi", "iq", "q_chi_frames", "iq_frames"]:
        obj = _get(key)
        if obj is None:
            print(f"  {key:20s}: None")
        elif hasattr(obj, "shape"):
            print(f"  {key:20s}: {obj.shape}")
        elif hasattr(obj, "dims"):
            print(f"  {key:20s}: dims={dict(obj.sizes)}")
        else:
            print(f"  {key:20s}: type={type(obj).__name__}")

    print("\n=== DONE ===")


if __name__ == "__main__":
    main()
