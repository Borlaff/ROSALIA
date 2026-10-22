"""
PyTorch acceleration module for the Roman SUPERBACK PSF pipeline.

This module provides high-performance GPU (CUDA, Apple MPS) and multi-threaded CPU
implementations of the SUPERBACK analytical PSF model for the Nancy Grace Roman Space Telescope.
It includes:
- Equal-stamp-size batched multiprocessing with size-ceiling quantizations (always >= required).
- Dynamic batch-sizing to optimize GPU core occupancy and available memory.
- Spatial chunking for gigantic stamps (e.g. m <= 2 bright stars) to bound memory and eliminate OOM errors.
- Inference-mode execution with zero gradient tracking.
- Compatibility across CUDA, MPS, and CPU.
"""

import os
import sys
import gc
import numpy as np
import torch
import astropy.units as u
from astropy.io import fits
# Remove tqdm to prevent thread conflict / SIGSEGV on macOS with PyTorch
# from tqdm import tqdm

# Automatically ensure MPS fallback is enabled on macOS if unassigned
if "PYTORCH_ENABLE_MPS_FALLBACK" not in os.environ:
    os.environ["PYTORCH_ENABLE_MPS_FALLBACK"] = "1"

# Ensure SUPERBACK_PATH points to ROSALIACACHE/CORE/superback if available
if "SUPERBACK_PATH" not in os.environ and "ROSALIACACHE" in os.environ:
    cache_sb = os.path.join(os.environ["ROSALIACACHE"], "CORE", "superback")
    if os.path.isdir(cache_sb):
        os.environ["SUPERBACK_PATH"] = cache_sb

import rosalia as rs
import rosalia.superback as superback

# Optical parameters for the Roman Space Telescope WFI
DEFAULT_FLENGTH = 18.813
DEFAULT_AINNER = 0.783 / 2.0
DEFAULT_AOUTER = 2.553 / 2.0
DEFAULT_PIXSIZE = 10e-6
ARCSEC_TO_RAD = np.pi / (180.0 * 3600.0)


def get_default_device(device=None):
    """
    Resolve the computation device across CUDA, MPS, and CPU.

    Parameters
    ----------
    device : str or torch.device, optional
        Target device ('auto', 'cuda', 'mps', 'cpu', or torch.device).
        If None or 'auto', automatically prioritizes CUDA -> MPS -> CPU.

    Returns
    -------
    torch.device
        Selected PyTorch device.
    """
    if device is None or device == "auto":
        if torch.cuda.is_available():
            return torch.device("cuda")
        elif torch.backends.mps.is_available():
            return torch.device("mps")
        return torch.device("cpu")
    if isinstance(device, str):
        return torch.device(device)
    return device


def resolve_device_and_dtype(device=None, dtype=None):
    """
    Resolve device and appropriate precision.
    Apple MPS does not support float64, so it defaults/falls back to float32.
    CUDA and CPU default to float64 for maximum astronomical precision.
    """
    device = get_default_device(device)
    if dtype is None:
        dtype = torch.float32 if device.type == "mps" else torch.float64
    elif device.type == "mps" and dtype == torch.float64:
        dtype = torch.float32
    return device, dtype


def get_device_memory_info(device=None, memory_fraction=0.70):
    """
    Query available free memory for the target hardware and calculate adaptive budgets.

    Parameters
    ----------
    device : str or torch.device, optional
        Target hardware device ('cuda', 'mps', 'cpu', or None).
    memory_fraction : float, default 0.70
        Fraction of available memory to allocate for the operational budget (typically 0.60 to 0.80).

    Returns
    -------
    dict
        Adaptive memory parameters:
        - 'free_bytes': int, raw available free memory in bytes
        - 'total_bytes': int, total device/system memory in bytes
        - 'budget_bytes': int, operational memory budget in bytes
        - 'memory_fraction': float, clamped fraction used (0.10 to 0.90)
        - 'chunk_size': int, adaptive spatial tile size for large stamps (e.g. 1024, 2048, 4096)
        - 'chunk_threshold': int, stamp size boundary above which stars are chunked
        - 'budget_scale': float, relative scaling factor for star batch sizes
    """
    device, _ = resolve_device_and_dtype(device, None)
    fraction = max(0.10, min(0.90, float(memory_fraction)))

    free_bytes = None
    total_bytes = None

    # 1. Query CUDA VRAM
    if device.type == "cuda" and torch.cuda.is_available():
        try:
            dev_idx = device.index if device.index is not None else torch.cuda.current_device()
            free_bytes, total_bytes = torch.cuda.mem_get_info(dev_idx)
        except Exception:
            pass

    # 2. Query Apple Silicon unified memory / CPU system RAM
    if free_bytes is None:
        try:
            import psutil
            vm = psutil.virtual_memory()
            free_bytes = vm.available
            total_bytes = vm.total
            if device.type == "mps" and hasattr(torch, "mps") and hasattr(torch.mps, "recommended_max_memory"):
                rec_max = torch.mps.recommended_max_memory()
                free_bytes = min(free_bytes, rec_max)
        except Exception:
            free_bytes = 4 * 1024**3
            total_bytes = 8 * 1024**3

    budget_bytes = int(free_bytes * fraction)

    # Determine adaptive chunk size and chunk threshold based on budget
    # 1 tile of size C requires ~ 250 * C^2 bytes peak transient memory
    if budget_bytes >= 8 * 1024**3:
        chunk_size = 4096
        chunk_threshold = 2048
    elif budget_bytes >= 2 * 1024**3:
        chunk_size = 2048
        chunk_threshold = 2048
    elif budget_bytes >= 600 * 1024**2:
        chunk_size = 1024
        chunk_threshold = 1024
    else:
        chunk_size = 512
        chunk_threshold = 512

    # Relative scaling for batch sizing compared to 2 GB baseline
    budget_scale = max(0.5, min(4.0, budget_bytes / (2.0 * 1024**3)))

    return {
        "free_bytes": free_bytes,
        "total_bytes": total_bytes,
        "budget_bytes": budget_bytes,
        "memory_fraction": fraction,
        "chunk_size": chunk_size,
        "chunk_threshold": chunk_threshold,
        "budget_scale": budget_scale,
    }


def create_psf_grid(size_pixels, scale_arcsec=0.108, oversample_factor=1,
                    device=None, dtype=torch.float64):
    """
    Create spatial coordinate grids directly as PyTorch tensors.

    Parameters
    ----------
    size_pixels : int
        Grid dimension (size_pixels x size_pixels).
    scale_arcsec : float, default 0.108
        Pixel scale in arcsec/pixel.
    oversample_factor : int or float, default 1
        Sub-pixel oversampling factor.
    device : torch.device, optional
        Target PyTorch device.
    dtype : torch.dtype, default torch.float64
        Tensor floating-point precision.

    Returns
    -------
    x, y : torch.Tensor
        2D coordinate grids of shape (size_pixels, size_pixels).
    """
    fine_size = int(size_pixels * oversample_factor)
    fine_scale = scale_arcsec / oversample_factor
    half_field = (size_pixels * scale_arcsec) / 2.0
    coords = torch.linspace(
        -half_field + fine_scale / 2.0,
        half_field - fine_scale / 2.0,
        fine_size,
        device=device,
        dtype=dtype,
    )
    x, y = torch.meshgrid(coords, coords, indexing="xy")
    return x, y


def bessel_j1_safe(z, device=None):
    """
    Evaluate J1(z) / z safely across CUDA, MPS, and CPU.

    At z = 0, lim J1(z)/z = 0.5.
    On Apple MPS, Bessel J1 is evaluated via CPU transfer to avoid missing operator warnings.
    """
    dev = z.device if device is None else (torch.device(device) if isinstance(device, str) else device)
    half_const = torch.tensor(0.5, dtype=z.dtype, device=dev)
    if dev.type == "mps":
        z_cpu = z.cpu()
        j1_cpu = torch.special.bessel_j1(z_cpu)
        j1 = j1_cpu.to(dev)
    else:
        j1 = torch.special.bessel_j1(z)

    return torch.where(z == 0, half_const, j1 / z)


def barpsf_pytorch(x_m, y_m, k, w, l, x0, y0, rotation=0.0):
    """
    Fourier transform of a rectangular mirror strut in PyTorch.

    Parameters
    ----------
    x_m, y_m : torch.Tensor
        Spatial coordinates in meters.
    k : torch.Tensor or float
        Wavenumber 2 * pi / lambda / flength.
    w, l : float
        Width and length of the strut in meters.
    x0, y0 : float
        Centroid offset of the strut in meters.
    rotation : float
        Strut orientation angle in radians.

    Returns
    -------
    torch.Tensor
        Complex Fourier transform pattern.
    """
    cos_r = torch.as_tensor(np.cos(rotation), dtype=x_m.dtype, device=x_m.device)
    sin_r = torch.as_tensor(np.sin(rotation), dtype=x_m.dtype, device=x_m.device)

    xrot = x_m * cos_r - y_m * sin_r
    yrot = x_m * sin_r + y_m * cos_r

    s1 = torch.sinc(k * (w / 2.0 / np.pi) * xrot)
    s2 = torch.sinc(k * (l / 2.0 / np.pi) * yrot)

    phase = -k * (x0 * x_m + y0 * y_m)
    complex_phase = torch.complex(torch.cos(phase), torch.sin(phase))

    return s1 * s2 * (l * w) * complex_phase


def trappsf_pytorch(x_m, y_m, k, a, b, c, x0, y0, rotation=0.0):
    """
    Fourier transform of a right trapezoidal strut in PyTorch.
    """
    cos_r = torch.as_tensor(np.cos(rotation), dtype=x_m.dtype, device=x_m.device)
    sin_r = torch.as_tensor(np.sin(rotation), dtype=x_m.dtype, device=x_m.device)

    xrot = x_m * cos_r - y_m * sin_r
    yrot = x_m * sin_r + y_m * cos_r

    s1 = torch.sinc(k * (c / 2.0 / np.pi) * xrot) * c
    h = a * (c - xrot) + b * xrot
    s2 = torch.sinc(k * h * yrot / (2.0 * c * np.pi)) * (h / c)

    c_phase = (c**2 * xrot + a * c * yrot + (b - a) * xrot * yrot / (2.0 * c)) + x0 * xrot + yrot * y0
    phase = -k * c_phase
    complex_phase = torch.complex(torch.cos(phase), torch.sin(phase))

    return s1 * s2 * complex_phase


def normalize_barparam(barparam):
    """
    Normalize strut dictionary or list into a standard sequence of bar dicts.
    """
    if barparam is None:
        barparam = superback.defaultbarparam

    if isinstance(barparam, dict):
        bars = []
        for k in ["1", "2", "3", "4", "5", "6"]:
            if k in barparam and isinstance(barparam[k], dict):
                b = dict(barparam[k])
                b.setdefault("type", "bar")
                bars.append(b)
        if not bars:
            bars = [v for k, v in barparam.items() if isinstance(v, dict) and "w" in v]
        return bars
    elif isinstance(barparam, (list, tuple)):
        return list(barparam)
    return []


def compute_anchor_weights(sample_wave_nm, flux_weights, wave_anchors_nm,
                           device=None, dtype=torch.float64):
    """
    Compute 1D projection weights beta for anchor wavelengths.

    Because linear interpolation is a linear operator:
        sum_j w_j * interp(lambda_j, anchors) = sum_k beta_k * anchor[k]
    This eliminates repetitive 2D array interpolations and accelerates chromatic integration.
    """
    s_wave = np.asarray(sample_wave_nm, dtype=float)
    f_weights = np.asarray(flux_weights, dtype=float)
    anchors = np.asarray(wave_anchors_nm, dtype=float)

    n_anchors = len(anchors)
    idx = np.searchsorted(anchors, s_wave) - 1
    idx = np.clip(idx, 0, n_anchors - 2)

    denom = anchors[idx + 1] - anchors[idx]
    t = (s_wave - anchors[idx]) / denom

    beta = np.zeros(n_anchors, dtype=np.float64)
    np.add.at(beta, idx, (1.0 - t) * f_weights)
    np.add.at(beta, idx + 1, t * f_weights)

    return torch.as_tensor(beta, dtype=dtype, device=device)


@torch.inference_mode()
def romanpsf_pytorch(xoffset, yoffset, wave, barparam=None, x0=0.0, y0=0.0,
                     flength=DEFAULT_FLENGTH, ainner=DEFAULT_AINNER,
                     aouter=DEFAULT_AOUTER, pixsize=DEFAULT_PIXSIZE,
                     device=None, dtype=torch.float64):
    """
    Calculate the Roman WFI complex PSF pattern using PyTorch.

    Supports batched wavelength evaluation across 1D tensor of wavelengths.
    """
    bars = normalize_barparam(barparam)
    device, dtype = resolve_device_and_dtype(device, dtype)

    if not isinstance(xoffset, torch.Tensor):
        xoffset_t = torch.as_tensor(xoffset, dtype=dtype, device=device)
    else:
        xoffset_t = xoffset.to(device=device, dtype=dtype)

    if not isinstance(yoffset, torch.Tensor):
        yoffset_t = torch.as_tensor(yoffset, dtype=dtype, device=device)
    else:
        yoffset_t = yoffset.to(device=device, dtype=dtype)

    is_scalar = not hasattr(wave, "__len__") or isinstance(wave, (str, bytes))
    if is_scalar:
        wave_t = torch.tensor([float(wave)], dtype=dtype, device=device).view(1, 1, 1)
    else:
        if not isinstance(wave, torch.Tensor):
            wave_t = torch.as_tensor(wave, dtype=dtype, device=device).view(-1, 1, 1)
        else:
            wave_t = wave.to(device=device, dtype=dtype).view(-1, 1, 1)

    # Physical coordinates in meters
    x_meters = (xoffset_t * ARCSEC_TO_RAD * flength).unsqueeze(0)
    y_meters = (yoffset_t * ARCSEC_TO_RAD * flength).unsqueeze(0)

    theta = torch.sqrt(x_meters**2 + y_meters**2)
    k = (2.0 * np.pi) / wave_t / flength

    z_out = k * aouter * theta
    z_in = k * ainner * theta

    j1_out = bessel_j1_safe(z_out, device)
    j1_in = bessel_j1_safe(z_in, device)

    mirror_out = j1_out * (2.0 * np.pi * aouter**2)
    mirror_in = j1_in * (2.0 * np.pi * ainner**2)

    complex_dtype = torch.complex128 if dtype == torch.float64 else torch.complex64
    final = (mirror_out - mirror_in).to(complex_dtype)

    bararea = []
    for bar in bars:
        amp = bar.get("amp", 1.0)
        bx0 = bar.get("x0", 0.0)
        by0 = bar.get("y0", 0.0)
        rot = bar.get("theta", 0.0)

        if bar.get("type", "bar") == "bar":
            w = bar["w"]
            l = bar["l"]
            b_val = barpsf_pytorch(x_meters, y_meters, k, w, l, bx0, by0, rotation=rot)
            final -= amp * b_val
            bararea.append(w * l)
        else:
            a = bar["a"]
            b_h = bar["b"]
            c = bar["base"]
            t_val = trappsf_pytorch(x_meters, y_meters, k, a, b_h, c, bx0, by0, rotation=rot)
            final -= amp * t_val
            bararea.append(c * (a + b_h) / 2.0)

    if x0 != 0.0 or y0 != 0.0:
        phase0 = -k * (x0 * x_meters + y0 * y_meters)
        final *= torch.complex(torch.cos(phase0), torch.sin(phase0))

    open_pupil_area = np.pi * (aouter**2 - ainner**2) - np.sum(bararea)
    pixscale_arcsec = pixsize / flength / ARCSEC_TO_RAD
    norm_factor = (k / 2.0 / np.pi / np.sqrt(open_pupil_area)) * (pixsize / pixscale_arcsec)
    final *= norm_factor

    if is_scalar:
        return final.squeeze(0)
    return final


@torch.inference_mode()
def superback_roman_psf_pytorch(position, SCA, bandpass, shape=None, SED=None,
                                magnitude_weights=None, wave_weights=None,
                                oversample_factor=1, flux_weights=None, zp=23.9,
                                pixel_scale=0.11, device=None, dtype=torch.float64,
                                num_anchors=10, return_tensor=False,
                                memory_fraction=0.70, chunk_size=None, **kwargs):
    """
    Render a flux-calibrated, chromatic Roman/WFI PSF with the SUPERBACK model in PyTorch.

    Returns an astropy PrimaryHDU matching ROSALIA conventions.
    """
    device, dtype = resolve_device_and_dtype(device, dtype)

    if magnitude_weights is not None and flux_weights is not None:
        raise ValueError("Provide either magnitude_weights or flux_weights, not both.")
    if (magnitude_weights is not None or flux_weights is not None) and wave_weights is None:
        raise ValueError("wave_weights must be provided with magnitude_weights or flux_weights.")
    if SED is not None and (magnitude_weights is not None or flux_weights is not None):
        raise ValueError("SED cannot be combined with explicit magnitude/flux weights.")

    # 1. Bandpass filter response
    bp = rs.telescopes.Roman.get_filter(instrument="WFI", filter_name=bandpass)
    filter_wave_nm = np.asarray(bp["wavelength_bins"].to(u.nm).value, dtype=float)
    filter_throughput = np.asarray(bp["transmission_bins"], dtype=float)
    order = np.argsort(filter_wave_nm)
    filter_wave_nm = filter_wave_nm[order]
    filter_throughput = filter_throughput[order]
    valid_filter = np.isfinite(filter_wave_nm) & np.isfinite(filter_throughput)
    filter_wave_nm = filter_wave_nm[valid_filter]
    filter_throughput = np.clip(filter_throughput[valid_filter], 0.0, None)

    filter_response = filter_throughput * np.gradient(filter_wave_nm)

    # 2. Weights
    if wave_weights is not None:
        if isinstance(wave_weights, u.Quantity):
            sample_wave_nm = np.asarray(wave_weights.to(u.nm).value, dtype=float)
        else:
            sample_wave_nm = np.asarray(wave_weights, dtype=float)

        if magnitude_weights is not None:
            sample_magnitudes = np.asarray(magnitude_weights, dtype=float)
            throughput = np.interp(sample_wave_nm, filter_wave_nm, filter_throughput, left=0, right=0)
            response = throughput * np.abs(np.gradient(sample_wave_nm))
            sample_flux = 10.0**(0.4 * (zp - sample_magnitudes))
            weights = sample_flux * response / response.sum()
        else:
            weights = np.asarray(flux_weights, dtype=float)
    elif SED is not None:
        sample_wave_nm = filter_wave_nm
        spectral_fnu = np.asarray(SED(sample_wave_nm), dtype=float)
        sample_flux = spectral_fnu / 1e-32
        response = filter_response / filter_response.sum()
        weights = sample_flux * response
    else:
        sample_wave_nm = filter_wave_nm
        weights = filter_response / filter_response.sum()

    valid_samples = np.isfinite(sample_wave_nm) & np.isfinite(weights) & (weights >= 0)
    sample_wave_nm = sample_wave_nm[valid_samples]
    weights = weights[valid_samples]

    if shape is None:
        ny = nx = 91
    else:
        ny, nx = map(int, shape)
    fov = max(ny, nx)

    # 3. Strut geometry from detector coordinates
    thetax_rad, thetay_rad = rs.psf.get_detector_position_telescope_frame(
        SCA, idl_x=position[0], idl_y=position[1]
    )
    bard = superback.getbardict(thetax_rad, thetay_rad)

    # 4. Monochromatic calculation across anchor wavelengths
    wave_anchors_nm = np.linspace(filter_wave_nm[0], filter_wave_nm[-1], num_anchors)
    beta = compute_anchor_weights(sample_wave_nm, weights, wave_anchors_nm, device=device, dtype=dtype)

    mem_info = get_device_memory_info(device=device, memory_fraction=memory_fraction)
    eff_chunk_size = chunk_size if chunk_size is not None else mem_info["chunk_size"]
    eff_chunk_threshold = mem_info["chunk_threshold"]

    if fov > eff_chunk_threshold:
        # Spatial chunking for gigantic single stamps to prevent OOM
        calc_dtype = torch.float32
        wave_anchors_m = torch.tensor(wave_anchors_nm * 1e-9, dtype=calc_dtype, device=device)
        beta_t = beta.to(device=device, dtype=calc_dtype)

        # Precompute normalization on central odd-sized core
        S_norm = min(fov, 2049)
        if S_norm % 2 == 0:
            S_norm += 1
        x_norm, y_norm = create_psf_grid(
            S_norm, scale_arcsec=0.108, oversample_factor=oversample_factor,
            device=device, dtype=calc_dtype
        )
        comp_norm = romanpsf_pytorch(x_norm, y_norm, wave_anchors_m, barparam=bard, device=device, dtype=calc_dtype)
        norm_sums = (torch.abs(comp_norm)**2).sum(dim=(-2, -1), keepdim=True)
        norm_sums = torch.where(norm_sums > 0, norm_sums, torch.ones_like(norm_sums))
        del x_norm, y_norm, comp_norm

        half_field = (fov * 0.108 / oversample_factor) / 2.0
        fine_scale = 0.108 / oversample_factor
        coords = torch.linspace(
            -half_field + fine_scale / 2.0,
            half_field - fine_scale / 2.0,
            fov,
            device=device,
            dtype=calc_dtype,
        )

        CHUNK_SIZE = eff_chunk_size
        psf_data = np.zeros((ny, nx), dtype=np.float32)
        y_crop_start = (fov - ny) // 2
        x_crop_start = (fov - nx) // 2

        chunks_y = list(range(0, ny, CHUNK_SIZE))
        chunks_x = list(range(0, nx, CHUNK_SIZE))

        for cy0 in chunks_y:
            cy1 = min(ny, cy0 + CHUNK_SIZE)
            sy0 = y_crop_start + cy0
            sy1 = y_crop_start + cy1
            for cx0 in chunks_x:
                cx1 = min(nx, cx0 + CHUNK_SIZE)
                sx0 = x_crop_start + cx0
                sx1 = x_crop_start + cx1

                x_grid, y_grid = torch.meshgrid(coords[sx0:sx1], coords[sy0:sy1], indexing="xy")
                comp_chunk = romanpsf_pytorch(
                    x_grid, y_grid, wave_anchors_m, barparam=bard,
                    device=device, dtype=calc_dtype
                )
                del x_grid, y_grid
                mono_chunk = torch.abs(comp_chunk)**2
                del comp_chunk
                mono_norm = mono_chunk / norm_sums
                del mono_chunk
                chrom_chunk = torch.sum(mono_norm * beta_t.view(-1, 1, 1), dim=0)
                del mono_norm

                psf_data[cy0:cy1, cx0:cx1] = chrom_chunk.detach().cpu().numpy()
                del chrom_chunk

        del coords, norm_sums, beta_t

        requested_flux = float(np.sum(weights))
        image_flux = float(np.sum(psf_data))
        if image_flux <= 0 or not np.isfinite(image_flux):
            raise RuntimeError("PyTorch SUPERBACK produced non-positive flux.")
        psf_data *= (requested_flux / image_flux)

        if return_tensor:
            return torch.from_numpy(psf_data).to(device=device, dtype=dtype)

        hdu = fits.PrimaryHDU(data=psf_data)
        hdu.header["TELESCOP"] = "ROMAN"
        hdu.header["INSTRUME"] = "WFI"
        hdu.header["DETECTOR"] = f"WFI{SCA:02d}_FULL"
        hdu.header["FILTER"] = bandpass
        hdu.header["X_POS"] = position[0]
        hdu.header["Y_POS"] = position[1]
        hdu.header["BUNIT"] = "nJy/pixel"
        hdu.header["MAGZP"] = (float(zp), "AB zero point for flux values in nJy")
        hdu.header["PIXSCALE"] = (float(pixel_scale), "Pixel scale in arcsec/pixel")
        hdu.header["TOTFLUX"] = (requested_flux, "Integrated PSF flux in nJy")
        hdu.header["ABMAG"] = (
            float(zp - 2.5 * np.log10(requested_flux)),
            "Integrated AB magnitude per square arcsec",
        )
        return hdu

    # Coordinate grid for standard-sized stamps (fov <= eff_chunk_threshold)
    x_src, y_src = create_psf_grid(
        fov, scale_arcsec=0.108, oversample_factor=oversample_factor,
        device=device, dtype=dtype
    )

    wave_anchors_m = torch.tensor(wave_anchors_nm * 1e-9, dtype=dtype, device=device)
    complex_psf = romanpsf_pytorch(
        x_src, y_src, wave_anchors_m, barparam=bard, device=device, dtype=dtype
    )
    mono_psf = torch.abs(complex_psf)**2
    sums = mono_psf.sum(dim=(-2, -1), keepdim=True)
    psf_cube = mono_psf / sums
    chromatic_psf = torch.sum(psf_cube * beta.view(-1, 1, 1), dim=0)
    del complex_psf, mono_psf, psf_cube

    # 6. Normalize to requested integrated flux
    requested_flux = float(np.sum(weights))
    image_flux = float(chromatic_psf.sum().item())
    if image_flux <= 0 or not np.isfinite(image_flux):
        raise RuntimeError("PyTorch SUPERBACK produced non-positive flux.")
    chromatic_psf *= (requested_flux / image_flux)

    # 7. Crop if rectangular shape requested
    if ny != fov or nx != fov:
        y_start = (fov - ny) // 2
        x_start = (fov - nx) // 2
        chromatic_psf = chromatic_psf[y_start:y_start + ny, x_start:x_start + nx]

    if return_tensor:
        return chromatic_psf

    psf_data = chromatic_psf.detach().cpu().numpy()
    hdu = fits.PrimaryHDU(data=psf_data)
    hdu.header["TELESCOP"] = "ROMAN"
    hdu.header["INSTRUME"] = "WFI"
    hdu.header["DETECTOR"] = f"WFI{SCA:02d}_FULL"
    hdu.header["FILTER"] = bandpass
    hdu.header["X_POS"] = position[0]
    hdu.header["Y_POS"] = position[1]
    hdu.header["BUNIT"] = "nJy/pixel"
    hdu.header["MAGZP"] = (float(zp), "AB zero point for flux values in nJy")
    hdu.header["PIXSCALE"] = (float(pixel_scale), "Pixel scale in arcsec/pixel")
    hdu.header["TOTFLUX"] = (requested_flux, "Integrated PSF flux in nJy")
    hdu.header["ABMAG"] = (
        float(zp - 2.5 * np.log10(requested_flux)),
        "Integrated AB magnitude per square arcsec",
    )
    return hdu


def partition_stars_by_size(star_items, n_workers=None, budget_scale=1.0, chunk_threshold=2048, is_cpu=False):
    """
    Partition stars into batches of identical stamp sizes.

    Guarantee:
    Every star in a batch is evaluated with stamp size S_batch = max(raw_sizes in batch).
    Thus, star stamp sizes are ALWAYS bigger than necessary, NEVER less.
    """
    # Sort stars by required raw size ascending
    star_items.sort(key=lambda s: s["raw_stamp_size"])

    batches = []
    i = 0
    n_stars = len(star_items)

    while i < n_stars:
        current_size = star_items[i]["raw_stamp_size"]

        # Dynamic batch sizing adapted to hardware memory budget
        if current_size <= 129:
            max_b = min(256, max(16, int(64 * budget_scale)))
        elif current_size <= 257:
            max_b = min(128, max(8, int(32 * budget_scale)))
        elif current_size <= 513:
            max_b = min(64, max(4, int(16 * budget_scale)))
        elif current_size <= 1025:
            max_b = min(32, max(2, int(8 * budget_scale)))
        elif current_size <= chunk_threshold:
            max_b = min(8, max(1, int(2 * budget_scale)))
        else:
            max_b = 1  # Gigantic stars spatially chunked individually

        if is_cpu and n_workers is not None and n_workers > 0 and current_size <= 1025:
            max_b = min(max_b, int(n_workers * 4))

        # Collect stars for this batch
        batch_stars = [star_items[i]]
        j = i + 1
        while j < n_stars and len(batch_stars) < max_b:
            next_size = star_items[j]["raw_stamp_size"]
            # Don't group if next star is more than 35% larger (keeps padding minimal)
            if next_size > current_size * 1.35 and len(batch_stars) >= 2:
                break
            batch_stars.append(star_items[j])
            j += 1

        # The batch stamp size is the maximum required size in this batch (odd integer)
        batch_size = max(s["raw_stamp_size"] for s in batch_stars)
        if batch_size % 2 == 0:
            batch_size += 1

        batches.append((batch_size, batch_stars))
        i = j

    return batches


@torch.inference_mode()
def render_batch_into_canvas(batch_size, stars, wave_anchors_nm, canvas,
                             device=None, dtype=torch.float64, pixel_scale=0.11,
                             verbose=False, pbar=None, batch_info=None,
                             chunk_size=2048, chunk_threshold=2048):
    """
    Render a batch of stars sharing the exact same stamp size and inject in-place into canvas.

    For gigantic stamps (S > chunk_threshold), stars are spatially chunked into tiles of size `chunk_size`
    intersecting the Roman canvas bounding box. This prevents any multi-gigabyte memory allocations
    and eliminates OOM / process kills while preserving exact flux conservation.
    """
    S = batch_size
    n_batch = len(stars)
    num_anchors = len(wave_anchors_nm)
    wave_anchors_m = torch.tensor(wave_anchors_nm * 1e-9, dtype=dtype, device=device)

    if S > chunk_threshold:
        # Spatial chunking for gigantic stamps (e.g. m <= 2 bright stars Alpha Cen, Sirius, etc.)
        calc_dtype = torch.float32
        wave_anchors_m = wave_anchors_m.to(calc_dtype)
        grid_scale = 0.108

        for star in stars:
            center_x = int(np.rint(star["x"]))
            center_y = int(np.rint(star["y"]))

            x_start = center_x - S // 2
            y_start = center_y - S // 2
            x_end = x_start + S
            y_end = y_start + S

            canvas_x0 = max(0, x_start)
            canvas_x1 = min(canvas.shape[1], x_end)
            canvas_y0 = max(0, y_start)
            canvas_y1 = min(canvas.shape[0], y_end)

            if canvas_x0 >= canvas_x1 or canvas_y0 >= canvas_y1:
                continue

            # 1. Normalization sum on central odd-sized core
            S_norm = min(S, 2049)
            if S_norm % 2 == 0:
                S_norm += 1
            x_norm, y_norm = create_psf_grid(S_norm, scale_arcsec=grid_scale, device=device, dtype=calc_dtype)
            comp_norm = romanpsf_pytorch(x_norm, y_norm, wave_anchors_m, barparam=star["bard"], device=device, dtype=calc_dtype)
            norm_sums = (torch.abs(comp_norm)**2).sum(dim=(-2, -1), keepdim=True)
            norm_sums = torch.where(norm_sums > 0, norm_sums, torch.ones_like(norm_sums))
            del x_norm, y_norm, comp_norm

            # 2. 1D coordinate tensor (only ~200 KB for S = 50,000)
            half_field = (S * grid_scale) / 2.0
            coords = torch.linspace(
                -half_field + grid_scale / 2.0,
                half_field - grid_scale / 2.0,
                S,
                device=device,
                dtype=calc_dtype,
            )

            # 3. Canvas-bounded tiles
            chunks_y = list(range(canvas_y0, canvas_y1, chunk_size))
            chunks_x = list(range(canvas_x0, canvas_x1, chunk_size))
            total_chunks = len(chunks_y) * len(chunks_x)

            beta = star["beta"].to(device=device, dtype=calc_dtype)
            chunk_idx = 0

            for cy0 in chunks_y:
                cy1 = min(canvas_y1, cy0 + chunk_size)
                sy0 = cy0 - y_start
                sy1 = cy1 - y_start
                for cx0 in chunks_x:
                    cx1 = min(canvas_x1, cx0 + chunk_size)
                    sx0 = cx0 - x_start
                    sx1 = cx1 - x_start
                    chunk_idx += 1

                    if verbose:
                        chunk_str = f"star {star['source_id']} ({S}x{S}) chunk {chunk_idx}/{total_chunks}"
                        if pbar is not None:
                            pbar.set_postfix_str(chunk_str)
                        elif chunk_idx % 20 == 1 or chunk_idx == total_chunks:
                            print(f"      ... evaluating chunk {chunk_idx}/{total_chunks} for stamp {S}x{S}...", flush=True)

                    x_grid, y_grid = torch.meshgrid(coords[sx0:sx1], coords[sy0:sy1], indexing="xy")
                    comp_chunk = romanpsf_pytorch(
                        x_grid, y_grid, wave_anchors_m, barparam=star["bard"],
                        device=device, dtype=calc_dtype
                    )
                    del x_grid, y_grid

                    mono_chunk = torch.abs(comp_chunk)**2
                    del comp_chunk

                    mono_norm = mono_chunk / norm_sums
                    del mono_chunk

                    chrom_chunk = torch.sum(mono_norm * beta.view(-1, 1, 1), dim=0)
                    del mono_norm

                    canvas[cy0:cy1, cx0:cx1] += chrom_chunk.detach().cpu().numpy()
                    del chrom_chunk

            del coords, norm_sums, beta

    else:
        # Precompute coordinate grid ONCE for the whole batch of matching stamp sizes
        calc_dtype = dtype
        x_src, y_src = create_psf_grid(S, scale_arcsec=0.108, device=device, dtype=calc_dtype)
        x_meters = (x_src * ARCSEC_TO_RAD * DEFAULT_FLENGTH).unsqueeze(0)
        y_meters = (y_src * ARCSEC_TO_RAD * DEFAULT_FLENGTH).unsqueeze(0)
        del x_src, y_src

        theta = torch.sqrt(x_meters**2 + y_meters**2)
        k = (2.0 * np.pi) / wave_anchors_m.to(calc_dtype).view(-1, 1, 1) / DEFAULT_FLENGTH  # (A, 1, 1)

        complex_dtype = torch.complex64 if calc_dtype == torch.float32 else torch.complex128

        # Precompute the common primary mirror Bessel pattern ONCE for all anchors
        z_out = k * DEFAULT_AOUTER * theta  # (A, S, S)
        z_in = k * DEFAULT_AINNER * theta

        j1_out = bessel_j1_safe(z_out, device)
        del z_out
        j1_in = bessel_j1_safe(z_in, device)
        del z_in

        m_out = j1_out * (2.0 * np.pi * DEFAULT_AOUTER**2)
        del j1_out
        m_in = j1_in * (2.0 * np.pi * DEFAULT_AINNER**2)
        del j1_in
        base_mirror = (m_out - m_in).to(complex_dtype)  # (A, S, S)
        del m_out, m_in

        pixscale_arcsec = DEFAULT_PIXSIZE / DEFAULT_FLENGTH / ARCSEC_TO_RAD

        for star in stars:
            final_star = base_mirror.clone()
            bars = normalize_barparam(star["bard"])

            bararea = []
            for bar in bars:
                amp = bar.get("amp", 1.0)
                bx0 = bar.get("x0", 0.0)
                by0 = bar.get("y0", 0.0)
                rot = bar.get("theta", 0.0)
                b_val = barpsf_pytorch(x_meters, y_meters, k, bar["w"], bar["l"], bx0, by0, rotation=rot)
                final_star -= amp * b_val
                bararea.append(bar["w"] * bar["l"])

            open_pupil_area = np.pi * (DEFAULT_AOUTER**2 - DEFAULT_AINNER**2) - np.sum(bararea)
            norm_factor = (k / 2.0 / np.pi / np.sqrt(open_pupil_area)) * (DEFAULT_PIXSIZE / pixscale_arcsec)
            final_star *= norm_factor

            mono_cube = torch.abs(final_star)**2  # (A, S, S)
            del final_star
            sums = mono_cube.sum(dim=(-2, -1), keepdim=True)
            mono_cube /= torch.where(sums > 0, sums, torch.ones_like(sums))

            chromatic_star = torch.sum(mono_cube * star["beta"].to(calc_dtype).view(-1, 1, 1), dim=0)
            del mono_cube

            req_flux = star["requested_flux"]
            c_sum = chromatic_star.sum()
            if c_sum > 0:
                chromatic_star *= (req_flux / c_sum)

            _inject_psf_patch(canvas, chromatic_star, star["x"], star["y"], S)
            del chromatic_star

        del base_mirror, x_meters, y_meters, theta, k

    if device.type == "cuda":
        torch.cuda.empty_cache()
    elif device.type == "mps":
        torch.mps.empty_cache()


def _inject_psf_patch(canvas, psf_tensor, star_x, star_y, S):
    """
    In-place injection of PSF patch into the NumPy canvas.
    """
    center_x = int(np.rint(star_x))
    center_y = int(np.rint(star_y))

    x_start = center_x - S // 2
    y_start = center_y - S // 2
    x_end = x_start + S
    y_end = y_start + S

    canvas_x_start = max(0, x_start)
    canvas_y_start = max(0, y_start)
    canvas_x_end = min(canvas.shape[1], x_end)
    canvas_y_end = min(canvas.shape[0], y_end)

    if canvas_x_start < canvas_x_end and canvas_y_start < canvas_y_end:
        psf_x_start = canvas_x_start - x_start
        psf_y_start = canvas_y_start - y_start
        psf_x_end = psf_x_start + (canvas_x_end - canvas_x_start)
        psf_y_end = psf_y_start + (canvas_y_end - canvas_y_start)

        patch = psf_tensor[psf_y_start:psf_y_end, psf_x_start:psf_x_end].detach().cpu().numpy()
        canvas[canvas_y_start:canvas_y_end, canvas_x_start:canvas_x_end] += patch
        del patch


@torch.inference_mode()
def build_psf_canvas_pytorch(stars_catalog, shape_out, bandpass, depth=31,
                             n_workers=None, device=None, dtype=torch.float64,
                             pixel_scale=0.11, num_anchors=10, verbose=True,
                             memory_fraction=0.70, chunk_size=None):
    """
    Build a full Roman WFI PSF canvas in PyTorch using size-bucketed multiprocessing.

    Parameters
    ----------
    stars_catalog : pandas.DataFrame
        Star catalog containing required astronomical & pixel coordinates.
    shape_out : tuple of int
        (ny, nx) shape of output canvas.
    bandpass : str
        Roman WFI filter name (e.g. 'F129', 'F158').
    depth : float, default 31
        Surface brightness depth for PSF extension in mag arcsec^-2.
    n_workers : int, optional
        Parallel CPU workers hint.
    device : str or torch.device, optional
        Target device ('auto', 'cuda', 'mps', 'cpu').
    dtype : torch.dtype, default torch.float64
        Calculation precision.
    pixel_scale : float, default 0.11
        Detector pixel scale in arcsec/pixel.
    num_anchors : int, default 10
        Anchor wavelengths for chromatic integration.
    verbose : bool, default True
        Whether to display progress.
    memory_fraction : float, default 0.70
        Fraction of available RAM/VRAM to allocate as the dynamic operational budget (typically 0.60 to 0.80).
    chunk_size : int, optional
        Override spatial chunk size (if None, adaptively computed based on available memory).

    Returns
    -------
    np.ndarray
        Accumulated PSF canvas of shape `shape_out`.
    """
    device, dtype = resolve_device_and_dtype(device, dtype)
    dtype = dtype or torch.float64

    # Configure CPU parallelism if executing on CPU
    cpu_threads = None
    if device.type == "cpu":
        if n_workers is not None and n_workers > 0:
            torch.set_num_threads(int(n_workers))
        cpu_threads = torch.get_num_threads()

    mem_info = get_device_memory_info(device=device, memory_fraction=memory_fraction)
    eff_chunk_size = chunk_size if chunk_size is not None else mem_info["chunk_size"]
    eff_chunk_threshold = mem_info["chunk_threshold"]
    budget_scale = mem_info["budget_scale"]

    canvas = np.zeros(shape_out, dtype=float)
    if stars_catalog.empty:
        return canvas

    # 1. Filter curves
    bp = rs.telescopes.Roman.get_filter(instrument="WFI", filter_name=bandpass)
    filter_wave_nm = np.asarray(bp["wavelength_bins"].to(u.nm).value, dtype=float)
    filter_throughput = np.asarray(bp["transmission_bins"], dtype=float)
    order = np.argsort(filter_wave_nm)
    filter_wave_nm = filter_wave_nm[order]
    filter_throughput = filter_throughput[order]
    valid_filter = np.isfinite(filter_wave_nm) & np.isfinite(filter_throughput)
    filter_wave_nm = filter_wave_nm[valid_filter]
    filter_throughput = np.clip(filter_throughput[valid_filter], 0.0, None)

    wave_anchors_nm = np.linspace(filter_wave_nm[0], filter_wave_nm[-1], num_anchors)

    # 2. Precompute star sizes and chromatic weights
    star_items = []
    for idx, row in stars_catalog.iterrows():
        x, y = float(row["x"]), float(row["y"])
        x_det, y_det = float(row["x_det"]), float(row["y_det"])
        detector_id = int(row["detector_id"])

        star = rs.psf.Star(row)
        magnitude = star.get_magnitudes(bandpass)
        extension = rs.psf.get_psf_extension_moffat(magnitude, depth=depth)
        radius_pixels = max(1, int(np.ceil(extension / pixel_scale)))
        raw_stamp_size = 2 * radius_pixels + 1

        wave_weights, flux_weights = star.get_linear_weights(bandpass)
        thetax_rad, thetay_rad = rs.psf.get_detector_position_telescope_frame(
            detector_id, idl_x=x_det, idl_y=y_det
        )
        bard = superback.getbardict(thetax_rad, thetay_rad)

        beta = compute_anchor_weights(
            wave_weights, flux_weights, wave_anchors_nm, device=device, dtype=dtype
        )

        star_items.append({
            "source_id": row.get("source_id", idx),
            "x": x,
            "y": y,
            "detector_id": detector_id,
            "raw_stamp_size": raw_stamp_size,
            "bard": bard,
            "beta": beta,
            "requested_flux": float(np.sum(flux_weights)),
        })

    # 3. Partition into batches of matching stamp sizes (always >= required)
    batches = partition_stars_by_size(
        star_items, n_workers=n_workers, budget_scale=budget_scale,
        chunk_threshold=eff_chunk_threshold, is_cpu=(device.type == "cpu")
    )
    num_batches = len(batches)
    total_stars = len(star_items)

    # 4. Dispatch batches across device with progress tracking
    if verbose:
        free_gb = mem_info["free_bytes"] / 1024**3
        total_gb = mem_info["total_bytes"] / 1024**3
        budget_gb = mem_info["budget_bytes"] / 1024**3
        pct = int(mem_info["memory_fraction"] * 100)
        avg_stars = total_stars / num_batches if num_batches > 0 else 0.0

        # Batch with most stars (if tied, pick the one with highest stamp size)
        max_stars_batch = max(batches, key=lambda b: (len(b[1]), b[0]))
        max_count = len(max_stars_batch[1])
        max_size = max_stars_batch[0]

        # Batch with smallest stamp size (if tied, pick the one with most stars)
        min_stamp_batch = min(batches, key=lambda b: (b[0], -len(b[1])))
        min_stamp_count = len(min_stamp_batch[1])
        min_stamp_size = min_stamp_batch[0]

        # Batch with fewest stars and largest stamp size
        min_star_count = min(len(b[1]) for b in batches)
        min_stars_batches = [b for b in batches if len(b[1]) == min_star_count]
        min_stars_biggest_batch = max(min_stars_batches, key=lambda b: b[0])
        min_count = len(min_stars_biggest_batch[1])
        min_size = min_stars_biggest_batch[0]
        min_note = " (spatially-chunked)" if min_size > eff_chunk_threshold else ""

        if device.type == "cpu":
            dev_str = f"CPU ({cpu_threads} OpenMP threads)"
        else:
            dev_str = f"{device.type.upper()}"

        print(f"\n[PyTorch SUPERBACK] Hardware device: {dev_str} ({dtype}) | Canvas: {shape_out[0]}x{shape_out[1]}", flush=True)
        print(f"[PyTorch SUPERBACK] Memory: {free_gb:.1f} GB free of {total_gb:.1f} GB | Adaptive budget ({pct}%): {budget_gb:.1f} GB (tile size: {eff_chunk_size}x{eff_chunk_size})", flush=True)
        print(f"[PyTorch SUPERBACK] Star catalog: {total_stars} stars | Batches: {num_batches} | Avg stars/batch: {avg_stars:.1f}", flush=True)
        if num_batches == 1:
            note = " (spatially-chunked)" if max_size > eff_chunk_threshold else ""
            print(f"  • Single batch: {max_count} star(s) | stamp size {max_size}x{max_size}{note}", flush=True)
        else:
            print(f"  • Batch with most stars: {max_count} star(s) | stamp size {max_size}x{max_size}", flush=True)
            print(f"  • Batch with smallest stamp: {min_stamp_count} star(s) | stamp size {min_stamp_size}x{min_stamp_size}", flush=True)
            print(f"  • Batch with fewest stars (largest stamp): {min_count} star(s) | stamp size {min_size}x{min_size}{min_note}", flush=True)

    pbar = None
    if verbose:
        try:
            from tqdm import tqdm
            pbar = tqdm(
                total=total_stars,
                desc=f"Injecting star PSFs ({device.type.upper()})",
                unit="star",
                leave=True,
            )
        except Exception:
            pbar = None

    try:
        for b_idx, (batch_size, stars) in enumerate(batches, 1):
            if pbar is not None:
                pbar.set_postfix_str(f"batch {b_idx}/{num_batches} ({batch_size}x{batch_size})")
            elif verbose:
                print(f"  --> Processing batch {b_idx}/{num_batches} ({len(stars)} star(s), stamp {batch_size}x{batch_size})...", flush=True)

            render_batch_into_canvas(
                batch_size, stars, wave_anchors_nm, canvas,
                device=device, dtype=dtype, pixel_scale=pixel_scale,
                verbose=verbose, pbar=pbar, batch_info=(b_idx, num_batches),
                chunk_size=eff_chunk_size, chunk_threshold=eff_chunk_threshold
            )
            if pbar is not None:
                pbar.update(len(stars))
    finally:
        if pbar is not None:
            pbar.close()

    if verbose:
        print(f"[PyTorch SUPERBACK] Completed injection of {total_stars} stars into canvas.\n", flush=True)

    gc.collect()
    return canvas
