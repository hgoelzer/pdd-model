program gpdd_mpi
  ! MPI-parallel PDD model: spatial domain decomposition along ny axis.
  ! Rank 0 handles all I/O; each rank owns a contiguous block of ny rows.
  ! massbalance_module is unchanged — called with local_ny instead of ny.

  use mpi
  use ncio
  use massbalance_module

  implicit none

  ! MPI
  integer :: nproc, rank, ierr, p
  integer, allocatable :: sendcounts(:), displs(:)  ! per-rank row counts * nx
  integer :: local_ny, jstart, jend

  REAL, PARAMETER :: rhoi = 910.0
  REAL, PARAMETER :: rhof = 1000.0

  character(len=256) :: filename
  character(len=256) :: fileroot_pr, fileroot_tas
  character(len=256) :: inpathname_pr, inpathname_tas
  character(len=256) :: outpathname, fileroot_out
  character(len=256) :: filename_prref, filename_tasref, filename_defmask
  character(len=7)   :: res_suffix
  character(len=4)   :: cyear
  character(len=256) :: nml_file
  character(len=256) :: institution
  character(len=256) :: contact_name
  character(len=256) :: contact_email

  integer :: nx, ny, nt
  double precision :: ddfactorsnow, ddfactorice, sigma, rainlimit
  integer :: nml_unit, nml_ios, narg
  integer :: i, j, t, year, year0, m, k
  integer :: ndims
  integer :: ncid
  integer :: fmode
  integer :: outputmode
  integer :: outputvars
  integer :: deflate_level
  double precision :: res

  ! dimensions (rank 0 only, for write_nc_file)
  double precision, allocatable :: x(:), y(:)

  ! Global arrays: allocated on rank 0 only, used for I/O
  double precision, allocatable :: tas_g(:,:,:)
  double precision, allocatable :: pr_g(:,:,:)
  double precision, allocatable :: smb_g(:,:,:)
  double precision, allocatable :: snow_g(:,:,:)
  double precision, allocatable :: rain_g(:,:,:)
  double precision, allocatable :: sir_g(:,:,:)
  double precision, allocatable :: abl_g(:,:,:)
  double precision, allocatable :: pdd_g(:,:,:)
  double precision, allocatable :: rfr_g(:,:,:)
  double precision, allocatable :: tas_ref_g(:,:,:)
  double precision, allocatable :: pr_ref_g(:,:,:)
  double precision, allocatable :: tas_anom_g(:,:,:)
  double precision, allocatable :: pr_ratio_g(:,:,:)
  double precision, allocatable :: def_mask_g(:,:)

  ! Transposed scatter/gather buffers — allocated once, reused every year.
  ! Layout (nx, 12, ny/local_ny) makes each rank's full data contiguous for a single MPI call.
  double precision, allocatable :: global_t(:,:,:)  ! rank 0: (nx, 12, ny)
  double precision, allocatable :: local_t(:,:,:)   ! all ranks: (nx, 12, local_ny)

  ! Local arrays: allocated on all ranks, size (nx, local_ny, 12)
  double precision, allocatable :: tas(:,:,:)
  double precision, allocatable :: pr(:,:,:)
  double precision, allocatable :: smb(:,:,:)
  double precision, allocatable :: snow(:,:,:)
  double precision, allocatable :: rain(:,:,:)
  double precision, allocatable :: sir(:,:,:)
  double precision, allocatable :: abl(:,:,:)
  double precision, allocatable :: pdd(:,:,:)
  double precision, allocatable :: rfr(:,:,:)
  double precision, allocatable :: tas_ref(:,:,:)
  double precision, allocatable :: pr_ref(:,:,:)
  double precision, allocatable :: tas_anom(:,:,:)
  double precision, allocatable :: pr_ratio(:,:,:)
  double precision, allocatable :: def_mask(:,:)

  character(len=32),  allocatable :: dimnames(:)
  integer,            allocatable :: dimlens(:)

  NAMELIST /config/ &
       outputmode, outputvars, fmode, deflate_level, &
       inpathname_pr, inpathname_tas, &
       fileroot_pr, fileroot_tas, &
       year0, nt, fileroot_out, &
       nx, ny, res, res_suffix, &
       filename_prref, filename_tasref, filename_defmask, &
       outpathname, &
       institution, contact_name, contact_email, &
       ddfactorsnow, ddfactorice, sigma, rainlimit

  ! ── MPI init ──────────────────────────────────────────────────────────────
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  ! ── Namelist (all ranks read; avoids broadcast of many scalars) ───────────
  outputmode    = 0
  outputvars    = 0
  fmode         = 1
  deflate_level = 0
  inpathname_pr = ""
  inpathname_tas= ""
  fileroot_pr   = ""
  fileroot_tas  = ""
  year0         = 2015
  nt            = 5
  fileroot_out  = ""
  nx            = 1681
  ny            = 2881
  res           = 1000.0d0
  res_suffix    = "i01000m"
  filename_prref   = ""
  filename_tasref  = ""
  filename_defmask = ""
  outpathname   = "./output"
  institution   = "NORCE Research"
  contact_name  = "Heiko Goelzer"
  contact_email = "heig@norceresearch.no"
  ddfactorsnow  = 0.00297d0
  ddfactorice   = 0.00791d0
  sigma         = 4.5d0
  rainlimit     = 1.0d0

  narg = command_argument_count()
  if (narg >= 1) then
     call get_command_argument(1, nml_file)
  else
     nml_file = "params.nml"
  end if

  nml_unit = 10
  open(unit=nml_unit, file=trim(nml_file), status="old", action="read", iostat=nml_ios)
  if (nml_ios /= 0) then
     if (rank == 0) write(*,*) "ERROR: Cannot open namelist file: ", trim(nml_file)
     call MPI_Finalize(ierr); stop 1
  end if
  read(nml_unit, nml=config, iostat=nml_ios)
  if (nml_ios /= 0) then
     if (rank == 0) write(*,*) "ERROR: Cannot read &config from: ", trim(nml_file)
     close(nml_unit); call MPI_Finalize(ierr); stop 1
  end if
  close(nml_unit)

  if (rank == 0) ncio_deflate_level = deflate_level

  if (rank == 0) then
     write(*,*) "## MPI tasks:     ", nproc
     write(*,*) "## fmode=", fmode, "  outputmode=", outputmode, "  outputvars=", outputvars
     write(*,*) "## year0=", year0, "  nt=", nt
     write(*,*) "## nx=", nx, "  ny=", ny, "  res=", res
     write(*,*) "## res_suffix=    ", trim(res_suffix)
     write(*,*) "## inpathname_pr= ", trim(inpathname_pr)
     write(*,*) "## inpathname_tas=", trim(inpathname_tas)
     write(*,*) "## outpathname=   ", trim(outpathname)
  end if

  ! ── Domain decomposition ──────────────────────────────────────────────────
  ! Block decomposition of ny along y-axis.
  ! The first mod(ny, nproc) ranks get one extra row.
  allocate(sendcounts(0:nproc-1), displs(0:nproc-1))
  displs(0) = 0
  do p = 0, nproc-1
     sendcounts(p) = ny / nproc
     if (p < mod(ny, nproc)) sendcounts(p) = sendcounts(p) + 1
     if (p > 0) displs(p) = displs(p-1) + sendcounts(p-1)
  end do
  local_ny = sendcounts(rank)
  jstart   = displs(rank) + 1
  jend     = jstart + local_ny - 1

  if (rank == 0) write(*,*) "## ny=", ny, "  local_ny per rank: min/max=", &
       minval(sendcounts), "/", maxval(sendcounts)

  ! ── Allocate scatter/gather buffers ─────────────────────────────────────
  allocate(local_t(nx, 12, local_ny))
  if (rank == 0) allocate(global_t(nx, 12, ny))

  ! ── Allocate global arrays (rank 0 only) ─────────────────────────────────
  if (rank == 0) then
     allocate(x(nx), y(ny))
     allocate(tas_g(nx,ny,12), pr_g(nx,ny,12))
     allocate(smb_g(nx,ny,12), snow_g(nx,ny,12), rain_g(nx,ny,12))
     allocate(sir_g(nx,ny,12), abl_g(nx,ny,12))
     allocate(pdd_g(nx,ny,12), rfr_g(nx,ny,12))
     allocate(def_mask_g(nx,ny))
     if (fmode == 1) then
        allocate(tas_ref_g(nx,ny,12), pr_ref_g(nx,ny,12))
        allocate(tas_anom_g(nx,ny,12), pr_ratio_g(nx,ny,12))
     end if
  end if

  ! ── Allocate local arrays (all ranks) ────────────────────────────────────
  allocate(tas(nx,local_ny,12), pr(nx,local_ny,12))
  allocate(smb(nx,local_ny,12), snow(nx,local_ny,12), rain(nx,local_ny,12))
  allocate(sir(nx,local_ny,12), abl(nx,local_ny,12))
  allocate(pdd(nx,local_ny,12), rfr(nx,local_ny,12))
  allocate(def_mask(nx,local_ny))
  if (fmode == 1) then
     allocate(tas_ref(nx,local_ny,12), pr_ref(nx,local_ny,12))
     allocate(tas_anom(nx,local_ny,12), pr_ratio(nx,local_ny,12))
  end if

  ! ── Read static fields (rank 0) and scatter ───────────────────────────────
  if (fmode == 1) then
     if (rank == 0) then
        write(*,*) "## Running in anomaly mode"
        ndims = nc_ndims(filename_tasref,"tas")
        call nc_dims(filename_tasref,"tas",dimnames,dimlens)
        write(*,*) "ndims= ", ndims, "  dimlens= ", dimlens
        call nc_open(filename_tasref, ncid, .FALSE.)
        call nc_read(filename_tasref, "tas", tas_ref_g)
        call nc_close(ncid)
        ndims = nc_ndims(filename_prref,"pr")
        call nc_open(filename_prref, ncid, .FALSE.)
        call nc_read(filename_prref, "pr", pr_ref_g)
        call nc_close(ncid)
        pr_ref_g = pr_ref_g / 31556926 * 12.
     end if
     ! Scatter refs to all ranks
     call scatter_3d(tas_ref_g, tas_ref)
     call scatter_3d(pr_ref_g,  pr_ref)
  else
     if (rank == 0) write(*,*) "## Running with tas, pr"
  end if

  ! Mask
  if (rank == 0) then
     ndims = nc_ndims(filename_defmask,"sftgif")
     call nc_dims(filename_defmask,"sftgif",dimnames,dimlens)
     write(*,*) "ndims= ", ndims, "  dimlens= ", dimlens
     call nc_open(filename_defmask, ncid, .FALSE.)
     call nc_read(filename_defmask, "sftgif", def_mask_g)
     call nc_close(ncid)
  end if
  call scatter_2d(def_mask_g, def_mask)

  ! ── Main year loop ────────────────────────────────────────────────────────
  t = 1
  DO WHILE(t <= nt)

     year = t + year0 - 1
     if (rank == 0) write(*,*) "Time counter ", t, year, nt
     write (cyear, "(I0.4)") year

     ! Reading tas
     filename = trim(inpathname_tas) // "/" // trim(fileroot_tas) // "_" // trim(cyear) // "_" // trim(res_suffix) // ".nc"
     if (rank == 0) then
        write(*,*) "### File tas: ", trim(filename)
        call nc_open(filename, ncid, .FALSE.)
        if (fmode == 0) then
           call nc_read(filename, "tas", tas_g)
        else
           call nc_read(filename, "tas", tas_anom_g)
        end if
        call nc_close(ncid)
     end if

     ! Reading pr
     filename = trim(inpathname_pr) // "/" // trim(fileroot_pr) // "_" // trim(cyear) // "_" // trim(res_suffix) // ".nc"
     if (rank == 0) then
        write(*,*) "### File pr: ", trim(filename)
        call nc_open(filename, ncid, .FALSE.)
        if (fmode == 0) then
           call nc_read(filename, "pr", pr_g)
        else
           call nc_read(filename, "pr", pr_ratio_g)
        end if
        call nc_close(ncid)
     end if

     ! Scatter forcing to all tasks
     if (fmode == 0) then
        call scatter_3d(tas_g, tas)
        call scatter_3d(pr_g,  pr)
     else
        call scatter_3d(tas_anom_g,  tas_anom)
        call scatter_3d(pr_ratio_g,  pr_ratio)
        ! Construct full local tas, pr from anomalies + references
        tas = tas_ref + tas_anom
        pr  = pr_ref  * pr_ratio
     end if

     ! Convert pr units: kg m-2 s-1 → m i.e. per month
     pr = pr * 31556926 / rhoi / 12.

     ! ── Physics (all ranks, on local arrays) ────────────────────────────────
     call pdd_model_greenland_total_monthly_inout(nx, local_ny, ddfactorsnow, ddfactorice, &
          sigma, rainlimit, pr, tas, smb, snow, rain, sir, abl, pdd, rfr)

     ! ── Gather output to rank 0 ─────────────────────────────────────────────
     call gather_3d(smb,  smb_g)
     call gather_3d(snow, snow_g)
     call gather_3d(rain, rain_g)
     call gather_3d(sir,  sir_g)
     call gather_3d(abl,  abl_g)
     call gather_3d(pdd,  pdd_g)
     call gather_3d(rfr,  rfr_g)
     call gather_3d(tas,  tas_g)
     call gather_3d(pr,   pr_g)

     ! ── Post-processing and write (rank 0 only) ─────────────────────────────
     if (rank == 0) then

        ! Unit conversion
        if (outputmode == 0) then
           smb_g  = smb_g  * rhoi / 31556926. * 12.
           pr_g   = pr_g   * rhoi / 31556926. * 12.
           rain_g = rain_g * rhoi / 31556926. * 12.
           snow_g = snow_g * rhoi / 31556926. * 12.
           abl_g  = abl_g  * rhoi / 31556926. * 12.
           sir_g  = sir_g  * rhoi / 31556926. * 12.
        else if (outputmode == 1) then
           tas_g  = tas_g  - 273.15
           smb_g  = smb_g  * rhoi
           pr_g   = pr_g   * rhoi
           rain_g = rain_g * rhoi
           snow_g = snow_g * rhoi
           abl_g  = abl_g  * rhoi
           sir_g  = sir_g  * rhoi
        end if

        ! Apply ice mask
        DO j = 1, ny
           DO i = 1, nx
              DO k = 1, 12
                 tas_g(i,j,k)  = tas_g(i,j,k)  * def_mask_g(i,j)
                 pdd_g(i,j,k)  = pdd_g(i,j,k)  * def_mask_g(i,j)
                 pr_g(i,j,k)   = pr_g(i,j,k)   * def_mask_g(i,j)
                 rfr_g(i,j,k)  = rfr_g(i,j,k)  * def_mask_g(i,j)
                 smb_g(i,j,k)  = smb_g(i,j,k)  * def_mask_g(i,j)
                 rain_g(i,j,k) = rain_g(i,j,k) * def_mask_g(i,j)
                 snow_g(i,j,k) = snow_g(i,j,k) * def_mask_g(i,j)
                 abl_g(i,j,k)  = abl_g(i,j,k)  * def_mask_g(i,j)
                 sir_g(i,j,k)  = sir_g(i,j,k)  * def_mask_g(i,j)
              END DO
           END DO
        END DO

        ! Write output
        if (outputmode == 0) then
           if (outputvars == 0) &
           call write_nc_file(tas_g, res, year, outpathname, &
                fileroot_out, "tas", "K", "air_temperature", institution, contact_name, contact_email)
           call write_nc_file(smb_g, res, year, outpathname, &
                fileroot_out, "acabf", "kg m-2 s-1", "land_ice_surface_specific_mass_balance_flux", &
                institution, contact_name, contact_email)
           if (outputvars == 0) then
              call write_nc_file(pr_g,   res, year, outpathname, fileroot_out, "pr", "kg m-2 s-1", "precipitation_flux", institution, contact_name, contact_email)
              call write_nc_file(abl_g,  res, year, outpathname, fileroot_out, "mrro", "kg m-2 s-1", "runoff_flux", institution, contact_name, contact_email)
              call write_nc_file(rain_g, res, year, outpathname, fileroot_out, "prra", "kg m-2 s-1", "rainfall_flux", institution, contact_name, contact_email)
              call write_nc_file(snow_g, res, year, outpathname, fileroot_out, "prsn", "kg m-2 s-1", "snowfall_flux", institution, contact_name, contact_email)
              call write_nc_file(sir_g,  res, year, outpathname, fileroot_out, "snicefreez", "kg m-2 s-1", "surface_snow_and_ice_refreezing_flux", institution, contact_name, contact_email)
           end if
        else if (outputmode == 1) then
           if (outputvars == 0) &
           call write_nc_file(tas_g, res, year, outpathname, &
                fileroot_out, "tas", "C", "air_temperature", institution, contact_name, contact_email)
           call write_nc_file(smb_g, res, year, outpathname, &
                fileroot_out, "acabf", "kg m-2 yr-1", "land_ice_surface_specific_mass_balance_flux", &
                institution, contact_name, contact_email)
           if (outputvars == 0) then
              call write_nc_file(pr_g,   res, year, outpathname, fileroot_out, "pr", "kg m-2 yr-1", "precipitation_flux", institution, contact_name, contact_email)
              call write_nc_file(abl_g,  res, year, outpathname, fileroot_out, "mrro", "kg m-2 yr-1", "runoff_flux", institution, contact_name, contact_email)
              call write_nc_file(rain_g, res, year, outpathname, fileroot_out, "prra", "kg m-2 yr-1", "rainfall_flux", institution, contact_name, contact_email)
              call write_nc_file(snow_g, res, year, outpathname, fileroot_out, "prsn", "kg m-2 yr-1", "snowfall_flux", institution, contact_name, contact_email)
              call write_nc_file(sir_g,  res, year, outpathname, fileroot_out, "snicefreez", "kg m-2 yr-1", "surface_snow_and_ice_refreezing_flux", institution, contact_name, contact_email)
           end if
        end if
        if (outputvars == 0) then
           call write_nc_file(pdd_g, res, year, outpathname, fileroot_out, "pdd", "1", "positive_degree_days", institution, contact_name, contact_email)
           call write_nc_file(rfr_g, res, year, outpathname, fileroot_out, "rfr", "1", "rain_fraction", institution, contact_name, contact_email)
        end if

     end if ! rank == 0

     t = t + 1
  END DO

  ! ── Clean up ───────────────────────────────────────────────────────────────
  deallocate(sendcounts, displs)
  deallocate(local_t)
  if (rank == 0) deallocate(global_t)
  deallocate(tas, pr, smb, snow, rain, sir, abl, pdd, rfr, def_mask)
  if (fmode == 1) deallocate(tas_ref, pr_ref, tas_anom, pr_ratio)
  if (rank == 0) then
     deallocate(x, y)
     deallocate(tas_g, pr_g, smb_g, snow_g, rain_g, sir_g, abl_g, pdd_g, rfr_g, def_mask_g)
     if (fmode == 1) deallocate(tas_ref_g, pr_ref_g, tas_anom_g, pr_ratio_g)
  end if

  call MPI_Finalize(ierr)

contains

  ! ── Internal scatter/gather helpers ─────────────────────────────────────────
  ! Scatter a global (nx, ny, 12) array from rank 0 to local (nx, local_ny, 12).
  ! Uses a transposed (nx, 12, ny) buffer so all 12 months for a rank's row-block
  ! are contiguous, enabling a single MPI_Scatterv instead of 12.
  subroutine scatter_3d(global, local)
    double precision, intent(in)  :: global(:,:,:)   ! (nx, ny, 12) on rank 0; ignored elsewhere
    double precision, intent(out) :: local(:,:,:)    ! (nx, local_ny, 12)
    integer :: j, sc(0:nproc-1), dp(0:nproc-1)
    sc = sendcounts * nx * 12
    dp = displs    * nx * 12
    if (rank == 0) then
       do j = 1, ny
          global_t(:, :, j) = global(:, j, :)
       end do
    end if
    call MPI_Scatterv(global_t, sc, dp, MPI_DOUBLE_PRECISION, &
                      local_t, local_ny*nx*12, MPI_DOUBLE_PRECISION, &
                      0, MPI_COMM_WORLD, ierr)
    do j = 1, local_ny
       local(:, j, :) = local_t(:, :, j)
    end do
  end subroutine scatter_3d

  ! Scatter a global (nx, ny) array from rank 0 to local (nx, local_ny).
  subroutine scatter_2d(global, local)
    double precision, intent(in)  :: global(:,:)
    double precision, intent(out) :: local(:,:)
    integer :: sc(0:nproc-1), dp(0:nproc-1)
    sc = sendcounts * nx
    dp = displs    * nx
    call MPI_Scatterv(global, sc, dp, MPI_DOUBLE_PRECISION, &
                      local, local_ny*nx, MPI_DOUBLE_PRECISION, &
                      0, MPI_COMM_WORLD, ierr)
  end subroutine scatter_2d

  ! Gather local (nx, local_ny, 12) from all ranks to global (nx, ny, 12) on rank 0.
  ! Uses the same transposed buffer trick as scatter_3d: single MPI_Gatherv per call.
  subroutine gather_3d(local, global)
    double precision, intent(in)  :: local(:,:,:)    ! (nx, local_ny, 12)
    double precision, intent(out) :: global(:,:,:)   ! (nx, ny, 12) on rank 0
    integer :: j, sc(0:nproc-1), dp(0:nproc-1)
    sc = sendcounts * nx * 12
    dp = displs    * nx * 12
    do j = 1, local_ny
       local_t(:, :, j) = local(:, j, :)
    end do
    call MPI_Gatherv(local_t, local_ny*nx*12, MPI_DOUBLE_PRECISION, &
                     global_t, sc, dp, MPI_DOUBLE_PRECISION, &
                     0, MPI_COMM_WORLD, ierr)
    if (rank == 0) then
       do j = 1, ny
          global(:, j, :) = global_t(:, :, j)
       end do
    end if
  end subroutine gather_3d

  ! ── write_nc_file (rank 0 only; identical to serial driver) ─────────────────
  subroutine write_nc_file(avar, res, year, outpathname, fileroot, varname, units, longname, &
                            institution, contact_name, contact_email)
    use ncio
    implicit none
    integer, parameter :: dp = KIND(1.0D0)
    real(dp), intent(in) :: avar(:,:,:), res
    integer,  intent(in) :: year
    character(len=*), intent(in) :: outpathname, fileroot, varname, units, longname
    character(len=*), intent(in) :: institution, contact_name, contact_email
    character(len=256) :: cyear, filename
    double precision :: days
    write(cyear, "(I0.4)") year
    days = (year-1850.0)*360.0+15
    filename = trim(outpathname) // "/" // trim(varname) // "_" // trim(fileroot) // "_" // trim(cyear) // ".nc"
    write(*,*) "### output file: ", trim(filename)
    call nc_create(filename, overwrite=.TRUE., netcdf4=.TRUE.)
    call nc_write_attr(filename, "title",            "NORCE-PDD output")
    call nc_write_attr(filename, "institution",      trim(institution))
    call nc_write_attr(filename, "contact_name",     trim(contact_name))
    call nc_write_attr(filename, "contact_email",    trim(contact_email))
    call nc_write_attr(filename, "grid_information", "ISMIP6 grid, epsg:3413")
    call nc_write_dim(filename, "x",    x=-720000.d0, dx=res,    nx=nx,    units="m")
    call nc_write_dim(filename, "y",    x=-3450000.d0, dx=res,   nx=ny,    units="m")
    call nc_write_dim(filename, "time", x=days, dx=30.d0,        nx=12, &
         units="days since 1850-01-01 00:00:00", calendar="360_day", unlimited=.TRUE.)
    call nc_write(filename, varname, avar(:,:,:), dim1="x", dim2="y", dim3="time", &
         units=units, long_name=longname)
  end subroutine write_nc_file

end program gpdd_mpi
