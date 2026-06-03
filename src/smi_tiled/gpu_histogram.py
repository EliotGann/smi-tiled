"""Optional PyTorch GPU-accelerated histogram binning backend.

This module provides a drop-in replacement for the sparse-mat-vec histogram
used in ``_SplitBinPlan.integrate_frame()``.  It accelerates the per-frame
scatter-add operation on CUDA (NVIDIA), MPS (Apple Metal), or CPU via PyTorch.

Usage
-----
The GPU plan is used transparently when ``torch`` is available and a device
is detected.  Import and instantiate:

    from smi_tiled.gpu_histogram import TorchBinPlan, get_torch_device

    device = get_torch_device()  # "cuda", "mps", or "cpu"
    plan = TorchBinPlan(q2d, chi2d, q_edges, chi_edges,
                        pixel_splitting=1, device=device)
    I_hist, N_hist = plan.integrate_frame(img, valid)

The output is identical to ``_SplitBinPlan`` (numpy arrays on CPU).

Requirements
------------
- ``torch >= 2.0`` (for MPS support)
- Install via: ``pip install torch`` or ``conda install pytorch -c pytorch``
"""

from __future__ import annotations

from typing import Any

import numpy as np


def is_available() -> bool:
    """Check if PyTorch is importable and a useful device exists."""
    try:
        import torch
        return (
            torch.cuda.is_available()
            or torch.backends.mps.is_available()
        )
    except ImportError:
        return False


def get_torch_device() -> str:
    """Return the best available PyTorch device string.

    Returns "cuda" (NVIDIA), "mps" (Apple Metal), or "cpu" (fallback).
    """
    try:
        import torch
    except ImportError:
        return "cpu"
    if torch.cuda.is_available():
        return "cuda"
    if torch.backends.mps.is_available():
        return "mps"
    return "cpu"


def _bin_indices(values: np.ndarray, edges: np.ndarray):
    """Assign each value to a bin index. Returns (indices, in_range_mask)."""
    idx = np.searchsorted(edges, values, side="right") - 1
    ok = (idx >= 0) & (idx < len(edges) - 1)
    return idx, ok


class TorchBinPlan:
    """GPU-accelerated precomputed pixel→bin mapping for histogram integration.

    This is a drop-in replacement for ``_SplitBinPlan`` that keeps the
    bin-index arrays on GPU and uses ``torch.scatter_add_`` for the
    per-frame integration.

    The construction is done on CPU (same logic as ``_SplitBinPlan``),
    then the index tensors are moved to the target device once.  Each
    ``integrate_frame()`` call uploads the image, does the scatter on
    device, and downloads the result.

    For best throughput with many frames, use ``integrate_batch()`` to
    process an entire block of frames in one kernel launch.
    """

    def __init__(
        self,
        q2d: np.ndarray,
        chi2d: np.ndarray,
        q_edges: np.ndarray,
        chi_edges: np.ndarray,
        pixel_splitting: int = 1,
        device: str | None = None,
    ) -> None:
        import torch

        self.n_q = len(q_edges) - 1
        self.n_chi = len(chi_edges) - 1
        self.n_bins = self.n_q * self.n_chi
        self.n_pixels = q2d.size
        self.shape = q2d.shape

        if device is None:
            device = get_torch_device()
        self.device = torch.device(device)

        n = max(int(pixel_splitting), 1)
        self.weight = 1.0 / (n * n)

        # Compute sub-pixel positions (same as _SplitBinPlan)
        if n == 1:
            sub_positions = [(q2d, chi2d)]
        else:
            dq_dr = np.gradient(q2d, axis=0)
            dq_dc = np.gradient(q2d, axis=1)
            dchi_dr = np.gradient(chi2d, axis=0)
            dchi_dc = np.gradient(chi2d, axis=1)
            offsets = np.linspace(-0.5 + 0.5 / n, 0.5 - 0.5 / n, n)
            sub_positions = []
            for dr in offsets:
                for dc in offsets:
                    q_sub = q2d + dr * dq_dr + dc * dq_dc
                    chi_sub = chi2d + dr * dchi_dr + dc * dchi_dc
                    sub_positions.append((q_sub, chi_sub))

        # For each sub-pixel offset, store:
        #   - bin_idx: int64 tensor of output bin indices (length = n_valid_pixels)
        #   - pix_idx: int64 tensor of source pixel indices (length = n_valid_pixels)
        # These are the "ok" pixels that fall within the grid.
        self._mappings: list[tuple[Any, Any]] = []
        for q_sub, chi_sub in sub_positions:
            qb, q_ok = _bin_indices(q_sub.ravel(), q_edges)
            cb, c_ok = _bin_indices(chi_sub.ravel(), chi_edges)
            ok = (
                q_ok & c_ok
                & np.isfinite(q_sub.ravel())
                & np.isfinite(chi_sub.ravel())
            )
            bin_idx = (qb[ok] * self.n_chi + cb[ok]).astype(np.int64)
            pix_idx = np.where(ok)[0].astype(np.int64)
            self._mappings.append((
                torch.from_numpy(bin_idx).to(self.device),
                torch.from_numpy(pix_idx).to(self.device),
            ))

        self._torch = torch

    def integrate_frame(
        self, img: np.ndarray, valid: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Integrate one frame. Returns (I_hist, N_hist) as numpy arrays.

        Parameters
        ----------
        img : 2-D array
            Intensity image (may contain NaN).
        valid : 2-D bool array
            Per-frame validity mask (True = include).
        """
        torch = self._torch

        # Prepare flat vectors on device
        img_flat = img.ravel().astype(np.float64)
        valid_flat = valid.ravel()
        # Zero invalid pixels so NaN doesn't propagate
        wI = np.where(valid_flat, img_flat, 0.0)
        wV = valid_flat.astype(np.float64)

        wI_t = torch.from_numpy(wI).to(self.device)
        wV_t = torch.from_numpy(wV).to(self.device)

        I_hist = torch.zeros(self.n_bins, dtype=torch.float64, device=self.device)
        N_hist = torch.zeros(self.n_bins, dtype=torch.float64, device=self.device)

        for bin_idx, pix_idx in self._mappings:
            # Gather pixel values for valid-in-grid pixels
            src_I = wI_t[pix_idx]
            src_V = wV_t[pix_idx]
            I_hist.scatter_add_(0, bin_idx, src_I)
            N_hist.scatter_add_(0, bin_idx, src_V)

        I_out = I_hist.cpu().numpy() * self.weight
        N_out = N_hist.cpu().numpy() * self.weight
        return I_out.reshape(self.n_q, self.n_chi), N_out.reshape(self.n_q, self.n_chi)

    def integrate_batch(
        self, images: np.ndarray, valids: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Integrate a batch of frames in one shot.

        Parameters
        ----------
        images : 3-D array, shape (n_frames, rows, cols)
            Stack of intensity images.
        valids : 3-D bool array, shape (n_frames, rows, cols)
            Per-frame validity masks.

        Returns
        -------
        I_hist : (n_frames, n_q, n_chi) array
        N_hist : (n_frames, n_q, n_chi) array
        """
        torch = self._torch
        n_frames = images.shape[0]

        # Flatten spatial dims: (n_frames, n_pixels)
        imgs_flat = images.reshape(n_frames, -1).astype(np.float64)
        vals_flat = valids.reshape(n_frames, -1).astype(np.float64)
        # Zero invalid pixels
        imgs_flat = np.where(vals_flat > 0, imgs_flat, 0.0)

        imgs_t = torch.from_numpy(imgs_flat).to(self.device)
        vals_t = torch.from_numpy(vals_flat).to(self.device)

        I_all = torch.zeros(n_frames, self.n_bins, dtype=torch.float64, device=self.device)
        N_all = torch.zeros(n_frames, self.n_bins, dtype=torch.float64, device=self.device)

        for bin_idx, pix_idx in self._mappings:
            # Gather columns for all frames: (n_frames, n_valid_pix)
            src_I = imgs_t[:, pix_idx]
            src_V = vals_t[:, pix_idx]
            # Expand bin_idx to match batch: (1, n_valid_pix) -> broadcast
            bin_exp = bin_idx.unsqueeze(0).expand(n_frames, -1)
            I_all.scatter_add_(1, bin_exp, src_I)
            N_all.scatter_add_(1, bin_exp, src_V)

        I_out = I_all.cpu().numpy() * self.weight
        N_out = N_all.cpu().numpy() * self.weight
        return (
            I_out.reshape(n_frames, self.n_q, self.n_chi),
            N_out.reshape(n_frames, self.n_q, self.n_chi),
        )
