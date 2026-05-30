# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build

**Local Mac** (conda env `pdd-build`):
```bash
make gpdd_monthly_inout LIB=/Users/heig/miniforge3/envs/pdd-build/lib INC=/Users/heig/miniforge3/envs/pdd-build/include
```
The conda env was created with `conda create -n pdd-build netcdf-fortran -c conda-forge`. Debug build: append `debug=1`.

**Cluster (Betzy / Olivia)** — both use the same module as of 2026-05:
```bash
source setenv.sh && make gpdd_monthly_inout
make clean   # removes executables, obj/*.o, obj/*.mod
```
`setenv.sh` is the single source of truth for the module environment (`netCDF-Fortran/4.6.1-gompi-2024a`). Source it before compiling — `runPDD` sources it automatically via `SLURM_SUBMIT_DIR`.

**Compiler quirk**: gfortran 13.3.0 (gompi-2024a on Betzy) treats lines >132 chars as a hard error (`-Werror=line-truncation`). Fixed with `-ffree-line-length-none` in `DFLAGS` (already in Makefile). The Mac conda gfortran doesn't enforce this.

**Makefile `LIBC`**: on clusters, netCDF-C and netCDF-Fortran live in separate module paths. `LIBC = $(if ${EBROOTNETCDF},${EBROOTNETCDF}/lib,$(LIB))` captures the C library path for `-L` and `-Wl,-rpath`. On Mac (conda) both libraries are co-located in `LIB`.

## Running

**Local Mac:**
```bash
DYLD_LIBRARY_PATH=/Users/heig/miniforge3/envs/pdd-build/lib ./gpdd_monthly_inout.x [params.nml]
```

**Cluster (SLURM):**
```bash
sbatch runPDD params_CESM2_historical_anom_16km.nml   # namelist as positional arg
```
Do NOT use `sbatch --export=NML=...` — restricting the exported environment breaks srun's ability to see module-set `LD_LIBRARY_PATH`. Pass the namelist as a positional argument instead; `runPDD` defaults to `params.nml` if none given.

**First-time setup on a new cluster directory** — replaces `Forcing/` and `output/` with symlinks to work storage:
```bash
./setup_work.sh   # creates /cluster/work/users/$USER/pdd/<setup>/{Forcing,output}
```
`setup_work.sh` targets `/cluster/work/users/$USER/...` (Betzy layout). On Olivia the work storage is project-based (`/cluster/work/projects/nn5011k/heig/pdd/`), so create symlinks manually there.

**Olivia**: use `runPDD_oli` instead of `runPDD` (`--partition=small`, no `--qos`). Olivia requires `--account=` and `--time=` explicitly in sbatch calls. No `devel` partition — use `small` for testing. Results bit-identical to Betzy (verified 2026-05-30).

The active driver is `gpdd_monthly_inout.x`. To switch scenarios pass a different `.nml` file — no recompile needed.

**Available namelist files:**

| File | Grid | Scenario |
|------|------|----------|
| `params.nml` | 1 km | MRI-ESM2-0 SSP585 anomaly (default) |
| `params_MRI-ESM2-0_historical_anom.nml` | 1 km | MRI-ESM2-0 historical |
| `params_CESM2_ssp370_anom.nml` | 1 km | CESM2 SSP370 anomaly |
| `params_CESM2_historical_anom.nml` | 1 km | CESM2 historical |
| `params_CESM2_ssp370_direct.nml` | 1 km | CESM2 SSP370 direct (fmode=0) |
| `params_CESM2_historical_anom_16km.nml` | 16 km | CESM2 historical, acabf only (outputvars=1) |
| `params_MRI-ESM2-0_historical_anom_16km.nml` | 16 km | MRI-ESM2-0 historical, acabf only |

**Namelist parameters** (all in `&config` group):
- `fmode` — 0 = direct tas/pr forcing; 1 = anomaly/ratio forcing (needs reference files)
- `outputmode` — 0 = SMBMIP units (kg m-2 s-1, K); 1 = human units (mm/month, C)
- `outputvars` — 0 = full output (all 9 variables, default); 1 = acabf only (faster for development)
- `outputfreq` — 0 = monthly output, 12 time steps per year file (default); 1 = annual mean output, 1 time step per year file (~12× smaller files; use for calibration ensembles)
- `year0`, `nt` — start year and number of years
- `inpathname_pr/tas`, `fileroot_pr/tas`, `res_suffix` — forcing file location and naming pattern
- `nx`, `ny`, `res` — grid dimensions and resolution (1681×2881 at 1 km; 211×361 at 8 km; 106×181 at 16 km)
- `filename_prref`, `filename_tasref` — reference fields for fmode=1 (not opened in fmode=0)
- `filename_defmask` — Greenland ice mask (`sftgif` variable); required for both fmodes
- `outpathname` — output directory (must exist before run)
- `institution`, `contact_name`, `contact_email` — written as global NetCDF attributes
- `ddfactorsnow`, `ddfactorice` — degree-day factors for snow and ice melt (m/yr/PDD); defaults 0.00297 and 0.00791
- `sigma` — temperature variability for PDD calculation (°C); default 4.5
- `rainlimit` — temperature threshold for rain/snow partitioning (°C); default 1.0
- `deflate_level` — zlib compression level for output NetCDF4 files (0 = off, 1–9); default 5 (~3× compression at 16 km, ~4–6× at 1 km due to large off-ice zero regions)

**Forcing files** live under `./Forcing/` (mirrored from cluster). The `output/` directory must exist before running.

## Architecture

Dependency chain: `ncio.f90` → `massbalance_module.f90` → driver programs.

**`ncio.f90`** — third-party NetCDF I/O wrapper (Alexander Robinson). Provides `nc_read`, `nc_write`, `nc_create`, `nc_write_dim`, `nc_write_attr`, `nc_ndims`, `nc_dims`, `nc_open`, `nc_close`. Two local modifications have been made: (1) `nc_create` now correctly passes `nf90_netcdf4` as the creation mode when `netcdf4=.TRUE.` — the original code had this line commented out since 2015, silently creating NetCDF3 files; (2) `ncio_deflate_level` module variable controls zlib compression (0 = off, 1–9 = deflate level); set by the driver after reading the namelist.

**`massbalance_module.f90`** — all PDD physics. Six subroutines in three tiers:

| Subroutine | Input | Output | Status |
|---|---|---|---|
| `pdd_model_greenland_total_monthly_inout` | monthly tp, t2m (3D) | monthly smb + components (3D) | **Active** — called by `gpdd_monthly_inout.f90` |
| `pdd_model_greenland_total_monthly` | monthly tp, t2m (3D) | annual smb + components (2D) | Legacy |
| `pdd_model_greenland_total_yearly` | annual acc, t2m, t2j (2D) | annual smb + components (2D) | Legacy |
| `massbalance_pdd_model_greenland` | anomaly-based, lat+Hs (2D) | annual smb (2D) | Legacy |
| `calculate_pdd_monthly_inout_taj` | monthly t2m (3D) | monthly pdd, rfr (3D) | Called by active path (uses seasonal parameterisation) |
| `calculate_pdd_monthly_inout` | monthly t2m (3D) | monthly pdd, rfr (3D) | Direct monthly variant (not currently called) |

All subroutines use `(nx, ny, ...)` array ordering. PDD lookup tables (`taberf`, `tabepdd`) are module-level `SAVE` arrays initialised once via `init_pdd_lut()` (private, called at the start of each `calculate_pdd_*` subroutine).

**Module-level constants** (`dp`, `pi`, `valmax`, `nintx`) are declared once above `CONTAINS` and shared by all subroutines.

**Key physics parameters** — all passed as `INTENT(IN)` arguments; set via namelist:
- `ddfactorsnow`, `ddfactorice`, `sigma`, `rainlimit` — exposed to all `pdd_model_*` and `calculate_pdd_*` subroutines; defaults in `gpdd_monthly_inout.f90` are 0.00297, 0.00791, 4.5, 1.0
- `pmax = 0.3` (refreezing cap) — still a `PARAMETER` constant in each `pdd_model_*` subroutine, not yet exposed

**Private helpers** (`init_pdd_lut`, `melt_cascade_2d`, `melt_cascade_3d`) — the snow/ice melt + refreezing cascade is extracted into these subroutines to eliminate duplication.

**`gpdd_monthly_inout.f90`** — active driver. Year loop reads one forcing file per year, calls `pdd_model_greenland_total_monthly_inout`, applies the ice mask, converts units, and writes one NetCDF file per variable per year via the `write_nc_file` contained subroutine. Output is on the ISMIP6 Greenland grid (epsg:3413, origin x=−720000, y=−3450000 m).

**`gpdd_monthly_inout_mpi.f90`** — MPI parallel driver. ny-axis domain decomposition; rank 0 handles all I/O. Uses transposed (nx,12,ny) buffers (`global_t`, `local_t`) to scatter/gather all 12 months in a single `MPI_Scatterv`/`MPI_Gatherv` call per variable (11 MPI collectives/year, down from 132). Build with `make gpdd_monthly_inout_mpi`; submit with `runPDD_mpi`.

**`gpdd_monthly.f90`**, **`gpdd.f90`** — older drivers, kept for reference. Not used in the current workflow.

## MPI parallelisation — conclusions

MPI within a single run does **not** pay off at any tested resolution. The bottleneck is serial I/O (rank 0 reads/writes all forcing and output) combined with memory transposition overhead in scatter/gather. Measured 1 km scaling on Betzy with the optimised code (11 MPI collectives/year):

| Tasks | Nodes | Queue | Elapsed | vs serial (1:01) |
|-------|-------|-------|---------|-----------------|
| 128 | 1 | preproc | 1:02 | ~same |
| 256 | 4 | normal | 1:08 | slightly slower |
| 512 | 4 | normal | 1:09 | slightly slower |

**Root causes of no speedup:** (1) Disk I/O is fully serial on rank 0 — reading ~900 MB/year of forcing dominates wall time. (2) The pack/unpack transposition loops in `scatter_3d`/`gather_3d` access `global(nx,ny,12)` with a 38 MB stride between months, thrashing L3 cache and adding overhead that offsets compute savings. (3) Compute is only ~15% of total run time even at 1 km.

**Right approach for throughput:** SLURM job arrays with serial runs — one job per model/scenario/parameter set. No code changes required; fully exploits the embarrassingly parallel nature of the ensemble.

**MPI submission quirks on Betzy:**
- Preproc queue: use `--mem=32G` (total node) not `--mem-per-cpu=8G` (128 tasks × 8G = 1024G overflows the 256G preproc node)
- Normal queue: minimum 4 nodes; memory is allocated automatically (no `--mem` needed)
- Preproc allows multiple jobs on the same node simultaneously — avoid submitting concurrent jobs that write to the same output directory
- Override `runPDD_mpi` defaults: `sbatch --partition=preproc --nodes=1 --ntasks=N --mem=32G runPDD_mpi params.nml`

## Resolution convergence

CESM2 historical 1950–1954, all resolutions with the same forcing and period:

| Res | Mean Gt/yr | Annual values |
|-----|-----------|---------------|
| 1 km | 274.9 | 200.0  224.9  335.3  239.8  374.5 |
| 2 km | 274.8 | 199.9  224.8  335.2  239.8  374.5 |
| 4 km | 275.4 | 200.3  225.5  335.8  240.4  375.1 |
| 8 km | 276.3 | 201.2  226.1  336.1  242.1  375.8 |
| 16 km | 273.4 | 198.1  223.7  333.1  239.3  373.1 |

Spread across all resolutions: **2.9 Gt/yr (~1%)**. Model is well-converged at 16 km for domain-integrated SMB. Plot script: `diag/plot_resolution_comparison.py` (run with `plotting` conda env).

## Refactoring status

**Done:**
- Configuration extracted from `gpdd_monthly_inout.f90` into runtime-read Fortran namelist files (`.nml`)
- `write_nc_file` now receives institution/contact as arguments instead of hardcoded strings
- `def_mask` read unconditionally (was only read inside `fmode==1`, leaving it uninitialised for `fmode==0`)
- `ddfactorsnow`, `ddfactorice`, `sigma`, `rainlimit` exposed as namelist parameters — removed from all `PARAMETER` declarations in `massbalance_module.f90` and threaded through all subroutine signatures
- Legacy drivers `gpdd.f90` and `gpdd_monthly.f90` call sites updated to pass the four physics parameters
- `dp`, `pi`, `valmax`, `nintx` promoted to module-level constants (were redeclared in every subroutine); `pi` corrected from `REAL` (single precision) to `REAL(dp)` (double precision)
- PDD lookup tables (`taberf`, `tabepdd`) promoted to module-level `SAVE` arrays; duplicated 20-line initialisation blocks replaced by a single private `init_pdd_lut()` subroutine
- Snow/ice melt + refreezing cascade extracted into private `melt_cascade_2d` and `melt_cascade_3d`, eliminating the four copies in `pdd_model_*` subroutines

**Remaining known issues in `massbalance_module.f90`**:
- `pmax = 0.3` still a local `PARAMETER` in each `pdd_model_*` subroutine; could be promoted to module level or namelist
- `calculate_pdd_monthly_inout_taj` exists "for testing backward compatibility" — unclear if it should become permanent or be replaced by `calculate_pdd_monthly_inout`
- `calculate_pdd_monthly_inout` (the non-taj variant) is dead code — defined but never called

## Numerical reproducibility

The refactoring was verified by comparing against the last committed baseline (pre Tasks 1–4). Domain-integrated totals agree to < 4×10⁻⁷ (machine precision) for all variables. Pixel-level breakdown:

| Category | Pixels | Magnitude | Cause |
|---|---|---|---|
| Floating-point noise | ~35% of grid | < 10⁻¹⁵ | Arithmetic reordering from subroutine extraction (Tasks 3 & 4) |
| LUT boundary rounding | ~200 pixels (0.0003%) | ~0.5 PDD, 0.002 rfr | `π` REAL→REAL(dp) shifts `nint()` at a handful of cold/warm boundaries |
| Negative pdd fix | 6.7 M pixels corrected | 0 → +O(10⁻⁸) | Single-precision `π` produced pdd ∈ [−2.5×10⁻⁸, 0) at the cold LUT edge; now correctly ≥ 0 |

The `output_ref/` directory holds output from before any refactoring; expect the same small differences there. The pixel-level changes are scientifically negligible — no systematic physics change.

## Testing

Reference output for the default scenario (MRI-ESM2-0 SSP585, 2015–2019) lives in `output_ref/`. After a run, compare numerically using the `nc` conda env:

```bash
/Users/heig/miniforge3/envs/nc/bin/python - <<'EOF'
import os, netCDF4 as nc, numpy as np
ref_dir, new_dir = "output_ref", "output"
for fname in sorted(f for f in os.listdir(ref_dir) if f.endswith(".nc")):
    r = nc.Dataset(f"{ref_dir}/{fname}")
    n = nc.Dataset(f"{new_dir}/{fname}")
    for vname in r.variables:
        if vname in ("x","y","time"): continue
        diff = np.max(np.abs(r.variables[vname][:] - n.variables[vname][:]))
        print(f"{fname}  {vname:20s}  {'IDENTICAL' if diff==0 else f'DIFFERS max={diff:.3e}'}")
    r.close(); n.close()
EOF
```

`output_ref/` was created before the `pi` REAL→REAL(dp) fix, so exact identity is no longer expected. Small differences in `pdd` (~0.5) and `rfr` (~0.002) at a handful of LUT boundary pixels are normal — see Numerical reproducibility below. To check for regressions, compare domain-integrated totals rather than pointwise equality.

The `nc` env (Python 3.13, netCDF4, scipy, numpy, cdo 2.4.4, nco) is the standard analysis environment. CLI tools like `ncdiff`, `ncdump`, `cdo` are also available in that env.
