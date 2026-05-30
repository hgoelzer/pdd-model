# PDD model

Implements the Positive Degree Day mass balance model for the Greenland ice sheet following Huybrechts and De Wolde (1999).

Uses [ncio](git@github.com:alex-robinson/ncio.git) by Alexander Robinson for NetCDF I/O.

## Setup

Forcing and grid data are typically linked or copied from the cluster storage. See `Forcing/README.txt` and `data/README.txt` for rsync instructions.

The `output/` directory must exist before running:

```bash
mkdir -p output
```

## Build

```bash
make gpdd_monthly_inout
```

On the cluster (SAGA/Betzy/Olivia) the NetCDF-Fortran paths are picked up automatically via `${EBROOTNETCDFMINFORTRAN}`. For local builds, override:

```bash
make gpdd_monthly_inout LIB=/path/to/netcdf-fortran/lib INC=/path/to/netcdf-fortran/include
```

## Configure

All run settings are in a Fortran namelist file. Copy and edit one of the provided templates:

| File | Scenario |
|------|----------|
| `params.nml` | MRI-ESM2-0 SSP585 anomaly mode (default) |
| `params_MRI-ESM2-0_historical_anom.nml` | MRI-ESM2-0 historical |
| `params_CESM2_ssp370_anom.nml` | CESM2 SSP370 anomaly |
| `params_CESM2_historical_anom.nml` | CESM2 historical |
| `params_CESM2_ssp370_direct.nml` | CESM2 SSP370 direct forcing |

Key settings in each `.nml` file:

```fortran
&config
  fmode        = 1        ! 0=direct tas/pr  1=anomaly/ratio forcing
  outputmode   = 0        ! 0=SMBMIP units (kg m-2 s-1)  1=human units (mm/month)
  year0        = 2015     ! start year
  nt           = 86       ! number of years
  nx           = 1681     ! grid x size (1 km resolution)
  ny           = 2881     ! grid y size

  ! PDD physics — tune without recompiling
  ddfactorsnow = 0.00297  ! degree-day factor for snow melt (m/yr/PDD)
  ddfactorice  = 0.00791  ! degree-day factor for ice melt  (m/yr/PDD)
  sigma        = 4.5      ! temperature variability (°C)
  rainlimit    = 1.0      ! rain/snow partitioning threshold (°C)
  ...
/
```

## Run

```bash
./gpdd_monthly_inout.x [params.nml]
```

If no argument is given, `params.nml` in the working directory is used. Pass a different file to switch scenarios without recompiling:

```bash
./gpdd_monthly_inout.x params_CESM2_ssp370_anom.nml
```

Output is written to `output/` — one NetCDF file per variable per year (tas, acabf, pr, mrro, prra, prsn, snicefreez, pdd, rfr).

## Diagnostics

Basin-integrated scalars and mass change diagnostics are in `diag/`:

```bash
./diag/makeDiag.sh      # generate diagnostic files
./diag/calc_scalars.sh  # basin-integrated scalar output
./diag/get_dM.sh        # mass change time series
./diag/get_dMap.sh      # mass change maps
```

See `diag/README.txt` for details.

## Testing

Reference output for the default scenario (MRI-ESM2-0 SSP585, 2015–2019) lives in `output_ref/`. After a run, compare numerically using the `nc` conda env (Python 3.13, netCDF4, numpy):

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

The `output_ref/` files were generated before the `pi` precision fix (see Code structure below); small differences in `pdd` and `rfr` at a handful of LUT boundary pixels are expected and physically negligible — domain-integrated totals agree to < 4×10⁻⁷.

## Code structure

All PDD physics lives in `massbalance_module.f90`. Module-level constants (`dp`, `pi`, `valmax`, `nintx`) and the PDD lookup tables (`taberf`, `tabepdd`) are declared once at module scope. Private helpers keep the public subroutines concise:

- `init_pdd_lut()` — initialises the PDD lookup tables on first call; shared by all `calculate_pdd_*` subroutines.
- `melt_cascade_2d` / `melt_cascade_3d` — snow/ice melt and refreezing loop; shared by all `pdd_model_*` subroutines.

Key physics parameters (`ddfactorsnow`, `ddfactorice`, `sigma`, `rainlimit`) are all namelist-configurable without recompiling. The legacy drivers `gpdd.f90` and `gpdd_monthly.f90` are kept for reference; the active driver is `gpdd_monthly_inout.f90`.
