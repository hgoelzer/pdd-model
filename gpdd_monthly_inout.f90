program gpdd
  ! pdd model for greenland with monthly input and output

  use ncio 
  use massbalance_module

  implicit none

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
  integer :: time_ref_year
  double precision :: ddfactorsnow, ddfactorice, sigma, rainlimit
  integer :: nml_unit, nml_ios, narg
  integer :: i, j, t, year, year0, m, k
  integer :: ndims
  integer :: ncid
  integer :: fmode
  integer :: outputmode
  integer :: outputvars
  integer :: outputfreq
  integer :: deflate_level

  ! Annual mean arrays (nx, ny, 1) — allocated when outputfreq == 1
  double precision, allocatable :: smb_ann(:,:,:)
  double precision, allocatable :: tas_ann(:,:,:)
  double precision, allocatable :: pr_ann(:,:,:)
  double precision, allocatable :: abl_ann(:,:,:)
  double precision, allocatable :: snow_ann(:,:,:)
  double precision, allocatable :: rain_ann(:,:,:)
  double precision, allocatable :: sir_ann(:,:,:)
  double precision, allocatable :: pdd_ann(:,:,:)
  double precision, allocatable :: rfr_ann(:,:,:)
  double precision :: res
  
  ! dimensions
  double precision,   allocatable :: x(:), y(:)

  ! variables, all monthly
  double precision,   allocatable :: tas(:,:,:)
  double precision,   allocatable :: pr(:,:,:)
  double precision,   allocatable :: smb(:,:,:)
  double precision,   allocatable :: snow(:,:,:)
  double precision,   allocatable :: rain(:,:,:)
  double precision,   allocatable :: sir(:,:,:)
  double precision,   allocatable :: abl(:,:,:)
  double precision,   allocatable :: pdd(:,:,:)
  double precision,   allocatable :: rfr(:,:,:)
  ! references, anomalies, rations
  double precision,   allocatable :: tas_ref(:,:,:)
  double precision,   allocatable :: pr_ref(:,:,:)
  double precision,   allocatable :: tas_anom(:,:,:)
  double precision,   allocatable :: pr_ratio(:,:,:)
  double precision,   allocatable :: def_mask(:,:)
  
  character(len=32),  allocatable :: dimnames(:)
  integer,            allocatable :: dimlens(:)

  NAMELIST /config/ &
       outputmode, outputvars, outputfreq, fmode, deflate_level, &
       inpathname_pr, inpathname_tas, &
       fileroot_pr, fileroot_tas, &
       year0, nt, time_ref_year, fileroot_out, &
       nx, ny, res, res_suffix, &
       filename_prref, filename_tasref, filename_defmask, &
       outpathname, &
       institution, contact_name, contact_email, &
       ddfactorsnow, ddfactorice, sigma, rainlimit

  ! Default configuration (overridden by namelist file)
  outputmode     = 0
  outputvars     = 0
  outputfreq     = 0
  fmode          = 1
  deflate_level  = 5
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
  time_ref_year = 1850

  ! Read namelist file: use first CLI argument, else default "params.nml"
  narg = command_argument_count()
  if (narg >= 1) then
     call get_command_argument(1, nml_file)
  else
     nml_file = "params.nml"
  end if
  write(*,*) "## Reading config from: ", trim(nml_file)
  nml_unit = 10
  open(unit=nml_unit, file=trim(nml_file), status="old", action="read", iostat=nml_ios)
  if (nml_ios /= 0) then
     write(*,*) "ERROR: Cannot open namelist file: ", trim(nml_file)
     stop 1
  end if
  read(nml_unit, nml=config, iostat=nml_ios)
  if (nml_ios /= 0) then
     write(*,*) "ERROR: Cannot read &config from: ", trim(nml_file)
     close(nml_unit)
     stop 1
  end if
  close(nml_unit)

  ncio_deflate_level = deflate_level

  write(*,*) "## fmode=", fmode, "  outputmode=", outputmode, "  outputvars=", outputvars, "  outputfreq=", outputfreq
  write(*,*) "## year0=", year0, "  nt=", nt, "  time_ref_year=", time_ref_year
  write(*,*) "## nx=", nx, "  ny=", ny, "  res=", res
  write(*,*) "## res_suffix=    ", trim(res_suffix)
  write(*,*) "## inpathname_pr= ", trim(inpathname_pr)
  write(*,*) "## inpathname_tas=", trim(inpathname_tas)
  write(*,*) "## fileroot_out=  ", trim(fileroot_out)
  write(*,*) "## outpathname=   ", trim(outpathname)
  write(*,*) "## institution=   ", trim(institution)

  ! #########################################################################
  
  ! Allocate dimensions
  allocate(x(nx),y(ny))

  ! Allocate arrays
  allocate(tas(nx,ny,12))
  allocate(pr(nx,ny,12))

  allocate(smb(nx,ny,12))
  allocate(snow(nx,ny,12))
  allocate(rain(nx,ny,12))
  allocate(sir(nx,ny,12))
  allocate(abl(nx,ny,12))
  allocate(pdd(nx,ny,12))
  allocate(rfr(nx,ny,12))

  allocate(def_mask(nx,ny))

  if (outputfreq == 1) then
     allocate(smb_ann(nx,ny,1))
     allocate(tas_ann(nx,ny,1))
     allocate(pr_ann(nx,ny,1))
     allocate(abl_ann(nx,ny,1))
     allocate(snow_ann(nx,ny,1))
     allocate(rain_ann(nx,ny,1))
     allocate(sir_ann(nx,ny,1))
     allocate(pdd_ann(nx,ny,1))
     allocate(rfr_ann(nx,ny,1))
  end if

  ! Mode specific
  if (fmode == 0) then
     ! Forcing with tas and pr
     write(*,*) "## Running with tas, pr"

  elseif (fmode == 1) then
     ! Run with anomalies/ratios; need references
     write(*,*) "## Running in anomaly mode"
     
     allocate(tas_ref(nx,ny,12))
     allocate(pr_ref(nx,ny,12))
     allocate(tas_anom(nx,ny,12))
     allocate(pr_ratio(nx,ny,12))
     
     ! read references
     !ndims = nc_ndims(filename_tasref,"tas")
     !call nc_dims(filename_tasref,"tas",dimnames,dimlens)
     !!write(*,*) "ndims= ", ndims
     !!write(*,*) "dimnames= ", dimnames
     !!write(*,*) "dimlens=  ", dimlens 
     !call nc_open(filename_tasref, ncid, .FALSE.)
     !call nc_read(filename_tasref, "tas",tas_ref)
     !call nc_close(ncid)

     ndims = nc_ndims(filename_tasref,"tas")
     call nc_dims(filename_tasref,"tas",dimnames,dimlens)
     write(*,*) "ndims= ", ndims
     write(*,*) "dimnames= ", dimnames
     write(*,*) "dimlens=  ", dimlens 
     call nc_open(filename_tasref, ncid, .FALSE.)
     call nc_read(filename_tasref, "tas",tas_ref) ! deg K
     call nc_close(ncid)
     ! alread in K
     !tas_ref(:,:,:) = tas_ref(:,:,:) + 273.15

     ndims = nc_ndims(filename_prref,"pr")
     call nc_dims(filename_prref,"pr",dimnames,dimlens)
     !write(*,*) "ndims= ", ndims
     !write(*,*) "dimnames= ", dimnames
     !write(*,*) "dimlens=  ", dimlens 
     call nc_open(filename_prref, ncid, .FALSE.)
     call nc_read(filename_prref, "pr",pr_ref) ! mm w.e. month-1
     call nc_close(ncid)
     ! convert to rate: mm w.e. s-1
     pr_ref(:,:,:) = pr_ref(:,:,:)/31556926*12.

  end if

  ! Read Greenland mask (required for both fmodes)
  ndims = nc_ndims(filename_defmask,"sftgif")
  call nc_dims(filename_defmask,"sftgif",dimnames,dimlens)
  write(*,*) "ndims= ", ndims
  write(*,*) "dimnames= ", dimnames
  write(*,*) "dimlens=  ", dimlens
  call nc_open(filename_defmask, ncid, .FALSE.)
  call nc_read(filename_defmask, "sftgif",def_mask) ! 1
  call nc_close(ncid)


  ! Main year loop
  t = 1
  DO WHILE(t <= nt)

     year = t + year0-1
     write(*,*) "Time counter ", t, year, nt

     ! Reading yearly files with monthly data
     write (cyear, "(I0.4)") year

     ! Reading tas [K]
     filename = trim(inpathname_tas) // "/" // trim(fileroot_tas) // "_" // trim(cyear) // "_" // trim(res_suffix) // ".nc"
     write(*,*) "### File tas: ", trim(filename)
     ndims = nc_ndims(filename,"tas")
     call nc_dims(filename,"tas",dimnames,dimlens)
     !write(*,*) "ndims= ", ndims
     !write(*,*) "dimnames= ", dimnames
     !write(*,*) "dimlens=  ", dimlens 
     call nc_open(filename, ncid, .FALSE.)
     if (fmode == 0) then
        call nc_read(filename, "tas",tas)
     else if (fmode == 1) then
        call nc_read(filename, "tas",tas_anom)
     end if
     call nc_close(ncid)

     ! Reading pr [kg m-2 s-1]
     filename = trim(inpathname_pr) // "/" // trim(fileroot_pr) // "_" // trim(cyear) // "_" // trim(res_suffix) // ".nc"
     write(*,*) "### File pr: ", trim(filename)
     ndims = nc_ndims(filename,"pr")
     call nc_dims(filename,"pr",dimnames,dimlens)
     !write(*,*) "ndims= ", ndims
     !write(*,*) "dimnames= ", dimnames
     !write(*,*) "dimlens=  ", dimlens 
     call nc_open(filename, ncid, .FALSE.)
     if (fmode == 0) then
        call nc_read(filename, "pr",pr)
     else if (fmode == 1) then
        call nc_read(filename, "pr",pr_ratio)
     end if
     call nc_close(ncid)
     ! End reading files   

     if (fmode == 1) then
        ! Forcing from anomalies/ratios
        ! Construct full forcing
        tas = tas_ref + tas_anom
        pr = pr_ref * pr_ratio
     end if
     
     ! Convert to units of PDD model: monthly total precip in m i.e.
     ! convert from kg m-2 s-1 = mm w.e. s-1  to  m i.e. per month
     pr(:,:,:) = pr(:,:,:)*31556926/rhoi/12.

     ! Model call
     call pdd_model_greenland_total_monthly_inout(nx, ny, ddfactorsnow, ddfactorice, sigma, rainlimit, pr, tas, smb, snow, rain, sir, abl, pdd, rfr)

     ! Convert for output depending on outputmode
     if (outputmode == 0) then
     ! Convert for output: temp in K
        !tas(:,:,:) = tas(:,:,:)
     ! Convert for output: SMB terms from monthly totals to rates kg m-2 s-1 = mm w.e. s-1
        smb(:,:,:) = smb(:,:,:) * rhoi / 31556926. * 12.
        pr(:,:,:) = pr(:,:,:) * rhoi / 31556926. * 12.
        rain(:,:,:) = rain(:,:,:) * rhoi / 31556926. * 12.
        snow(:,:,:) = snow(:,:,:) * rhoi / 31556926. * 12.
        abl(:,:,:) = abl(:,:,:) * rhoi / 31556926. * 12.
        sir(:,:,:) = sir(:,:,:) * rhoi / 31556926. * 12.
     else if (outputmode == 1) then
     ! Convert for output: temp in C
        tas(:,:,:) = tas(:,:,:)-273.15
     ! Convert for output: SMB terms to kg m-2 per month = mm w.e. per month
        smb(:,:,:) = smb(:,:,:) * rhoi 
        pr(:,:,:) = pr(:,:,:) * rhoi 
        rain(:,:,:) = rain(:,:,:) * rhoi 
        snow(:,:,:) = snow(:,:,:) * rhoi 
        abl(:,:,:) = abl(:,:,:) * rhoi 
        sir(:,:,:) = sir(:,:,:) * rhoi 
     end if

     ! Mask output to def_mask
    DO j=1,ny
       DO i=1,nx
          DO k=1,12
             tas(i,j,k) = tas(i,j,k) * def_mask(i,j)
             pdd(i,j,k) = pdd(i,j,k) * def_mask(i,j)
             pr(i,j,k) = pr(i,j,k) * def_mask(i,j)
             rfr(i,j,k) = rfr(i,j,k) * def_mask(i,j)
             smb(i,j,k) = smb(i,j,k) * def_mask(i,j)
             rain(i,j,k) = rain(i,j,k) * def_mask(i,j)
             snow(i,j,k) = snow(i,j,k) * def_mask(i,j)
             abl(i,j,k) = abl(i,j,k) * def_mask(i,j)
             sir(i,j,k) = sir(i,j,k) * def_mask(i,j)
          END DO
       END DO
    END DO
          
     ! Write files
     if (outputfreq == 1) then
        ! Annual mean output: average 12 months → 1 time step per year
        smb_ann(:,:,1)  = sum(smb(:,:,:),  dim=3) / 12.d0
        if (outputvars == 0) then
           tas_ann(:,:,1)  = sum(tas(:,:,:),  dim=3) / 12.d0
           pr_ann(:,:,1)   = sum(pr(:,:,:),   dim=3) / 12.d0
           abl_ann(:,:,1)  = sum(abl(:,:,:),  dim=3) / 12.d0
           snow_ann(:,:,1) = sum(snow(:,:,:), dim=3) / 12.d0
           rain_ann(:,:,1) = sum(rain(:,:,:), dim=3) / 12.d0
           sir_ann(:,:,1)  = sum(sir(:,:,:),  dim=3) / 12.d0
           pdd_ann(:,:,1)  = sum(pdd(:,:,:),  dim=3) / 12.d0
           rfr_ann(:,:,1)  = sum(rfr(:,:,:),  dim=3) / 12.d0
        end if
        if (outputmode == 0) then
           if (outputvars == 0) &
           call write_nc_file_annual(tas_ann, res, year, outpathname, &
                fileroot_out, "tas", "K", "air_temperature", &
                institution, contact_name, contact_email)
           call write_nc_file_annual(smb_ann, res, year, outpathname, &
                fileroot_out, "acabf", "kg m-2 s-1", "land_ice_surface_specific_mass_balance_flux", &
                institution, contact_name, contact_email)
           if (outputvars == 0) then
           call write_nc_file_annual(pr_ann, res, year, outpathname, &
                fileroot_out, "pr", "kg m-2 s-1", "precipitation_flux", &
                institution, contact_name, contact_email)
           call write_nc_file_annual(abl_ann, res, year, outpathname, &
                fileroot_out, "mrro", "kg m-2 s-1", "runoff_flux", &
                institution, contact_name, contact_email)
           call write_nc_file_annual(rain_ann, res, year, outpathname, &
                fileroot_out, "prra", "kg m-2 s-1", "rainfall_flux", &
                institution, contact_name, contact_email)
           call write_nc_file_annual(snow_ann, res, year, outpathname, &
                fileroot_out, "prsn", "kg m-2 s-1", "snowfall_flux", &
                institution, contact_name, contact_email)
           call write_nc_file_annual(sir_ann, res, year, outpathname, &
                fileroot_out, "snicefreez", "kg m-2 s-1", "surface_snow_and_ice_refreezing_flux", &
                institution, contact_name, contact_email)
           end if
        else if (outputmode == 1) then
           if (outputvars == 0) &
           call write_nc_file_annual(tas_ann, res, year, outpathname, &
                fileroot_out, "tas", "C", "air_temperature", &
                institution, contact_name, contact_email)
           call write_nc_file_annual(smb_ann, res, year, outpathname, &
                fileroot_out, "acabf", "kg m-2 yr-1", "land_ice_surface_specific_mass_balance_flux", &
                institution, contact_name, contact_email)
           if (outputvars == 0) then
           call write_nc_file_annual(pr_ann, res, year, outpathname, &
                fileroot_out, "pr", "kg m-2 yr-1", "precipitation_flux", &
                institution, contact_name, contact_email)
           call write_nc_file_annual(abl_ann, res, year, outpathname, &
                fileroot_out, "mrro", "kg m-2 yr-1", "runoff_flux", &
                institution, contact_name, contact_email)
           call write_nc_file_annual(rain_ann, res, year, outpathname, &
                fileroot_out, "prra", "kg m-2 yr-1", "rainfall_flux", &
                institution, contact_name, contact_email)
           call write_nc_file_annual(snow_ann, res, year, outpathname, &
                fileroot_out, "prsn", "kg m-2 yr-1", "snowfall_flux", &
                institution, contact_name, contact_email)
           call write_nc_file_annual(sir_ann, res, year, outpathname, &
                fileroot_out, "snicefreez", "kg m-2 yr-1", "surface_snow_and_ice_refreezing_flux", &
                institution, contact_name, contact_email)
           end if
        end if
        if (outputvars == 0) then
        call write_nc_file_annual(pdd_ann, res, year, outpathname, &
             fileroot_out, "pdd", "1", "positive_degree_days", &
             institution, contact_name, contact_email)
        call write_nc_file_annual(rfr_ann, res, year, outpathname, &
             fileroot_out, "rfr", "1", "rain_fraction", &
             institution, contact_name, contact_email)
        end if

     else
        ! Monthly output (default)
        if (outputmode == 0) then
           if (outputvars == 0) &
           call write_nc_file(tas, res, year, outpathname, &
                fileroot_out, "tas", "K", "air_temperature", &
                institution, contact_name, contact_email)
           call write_nc_file(smb, res, year, outpathname, &
                fileroot_out, "acabf", "kg m-2 s-1", "land_ice_surface_specific_mass_balance_flux", &
                institution, contact_name, contact_email)
           if (outputvars == 0) then
           call write_nc_file(pr, res, year, outpathname, &
                fileroot_out, "pr", "kg m-2 s-1", "precipitation_flux", &
                institution, contact_name, contact_email)
           call write_nc_file(abl, res, year, outpathname, &
                fileroot_out, "mrro", "kg m-2 s-1", "runoff_flux", &
                institution, contact_name, contact_email)
           call write_nc_file(rain, res, year, outpathname, &
                fileroot_out, "prra", "kg m-2 s-1", "rainfall_flux", &
                institution, contact_name, contact_email)
           call write_nc_file(snow, res, year, outpathname, &
                fileroot_out, "prsn", "kg m-2 s-1", "snowfall_flux", &
                institution, contact_name, contact_email)
           call write_nc_file(sir, res, year, outpathname, &
                fileroot_out, "snicefreez", "kg m-2 s-1", "surface_snow_and_ice_refreezing_flux", &
                institution, contact_name, contact_email)
           end if
        else if (outputmode == 1) then
           if (outputvars == 0) &
           call write_nc_file(tas, res, year, outpathname, &
                fileroot_out, "tas", "C", "air_temperature", &
                institution, contact_name, contact_email)
           call write_nc_file(smb, res, year, outpathname, &
                fileroot_out, "acabf", "kg m-2 yr-1", "land_ice_surface_specific_mass_balance_flux", &
                institution, contact_name, contact_email)
           if (outputvars == 0) then
           call write_nc_file(pr, res, year, outpathname, &
                fileroot_out, "pr", "kg m-2 yr-1", "precipitation_flux", &
                institution, contact_name, contact_email)
           call write_nc_file(abl, res, year, outpathname, &
                fileroot_out, "mrro", "kg m-2 yr-1", "runoff_flux", &
                institution, contact_name, contact_email)
           call write_nc_file(rain, res, year, outpathname, &
                fileroot_out, "prra", "kg m-2 yr-1", "rainfall_flux", &
                institution, contact_name, contact_email)
           call write_nc_file(snow, res, year, outpathname, &
                fileroot_out, "prsn", "kg m-2 yr-1", "snowfall_flux", &
                institution, contact_name, contact_email)
           call write_nc_file(sir, res, year, outpathname, &
                fileroot_out, "snicefreez", "kg m-2 yr-1", "surface_snow_and_ice_refreezing_flux", &
                institution, contact_name, contact_email)
           end if
        end if
        if (outputvars == 0) then
        call write_nc_file(pdd, res, year, outpathname, &
             fileroot_out, "pdd", "1", "positive_degree_days", &
             institution, contact_name, contact_email)
        call write_nc_file(rfr, res, year, outpathname, &
             fileroot_out, "rfr", "1", "rain_fraction", &
             institution, contact_name, contact_email)
        end if
     end if
     
     ! update timer
     t = t+1

  END DO ! end main loop


  ! Clean up
  deallocate(x,y)

  deallocate(tas)
  deallocate(pr)

  deallocate(smb)
  deallocate(snow)
  deallocate(rain)
  deallocate(sir)
  deallocate(abl)
  deallocate(pdd)
  deallocate(rfr)

  deallocate(def_mask)

  if (fmode == 1) then
     deallocate(tas_ref)
     deallocate(pr_ref)
     deallocate(tas_anom)
     deallocate(pr_ratio)
  end if

  if (outputfreq == 1) then
     deallocate(smb_ann, tas_ann, pr_ann, abl_ann)
     deallocate(snow_ann, rain_ann, sir_ann, pdd_ann, rfr_ann)
  end if

  
contains

  subroutine write_nc_file(avar,res,year,outpathname,fileroot,varname,units,longname, &
                            institution,contact_name,contact_email)
    ! Write variable to netcdf file on ISMIP6 grid

    use ncio

    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(1.0D0)  ! Kind of double precision numbers.

    REAL(dp), INTENT(IN)  :: avar(:,:,:)
    REAL(dp), INTENT(IN)  :: res ! 1000., 8000.
    INTEGER, INTENT(IN)  :: year
    character(len=*), INTENT(IN) :: outpathname
    character(len=*), INTENT(IN) :: fileroot
    character(len=*), INTENT(IN) :: varname
    character(len=*), INTENT(IN) :: units
    character(len=*), INTENT(IN) :: longname
    character(len=*), INTENT(IN) :: institution
    character(len=*), INTENT(IN) :: contact_name
    character(len=*), INTENT(IN) :: contact_email
    ! Internal variables
    character(len=256) :: cyear
    character(len=256) :: filename
    character(len=40)  :: time_units
    double precision :: days

    ! Construct output filename
    write (cyear, "(I0.4)") year
    days = (year - dble(time_ref_year)) * 360.0d0 + 15.0d0
    write (time_units, "(A,I0.4,A)") "days since ", time_ref_year, "-01-01 00:00:00"
    filename = trim(outpathname) // "/" // trim(varname) &
         // "_" // trim(fileroot) // "_" // trim(cyear) // ".nc"
    write(*,*) "### output file: ", trim(filename)

    ! Create the netcdf file, write global attributes
    call nc_create(filename,overwrite=.TRUE.,netcdf4=.TRUE.)
    call nc_write_attr(filename,"title","NORCE-PDD output")
    call nc_write_attr(filename,"institution", trim(institution))
    call nc_write_attr(filename,"contact_name", trim(contact_name))
    call nc_write_attr(filename,"contact_email", trim(contact_email))
    call nc_write_attr(filename,"grid_information", "ISMIP6 grid, epsg:3413")

    ! Write the dimensions (x, y, time), defined inline
    call nc_write_dim(filename,"x",x=-720000.d0,dx=res,nx=nx,units="m")
    call nc_write_dim(filename,"y",x=-3450000.d0,dx=res,nx=ny,units="m")
    call nc_write_dim(filename,"time",x=days,dx=30.d0,nx=12, &
         units=trim(time_units),calendar="360_day", unlimited=.TRUE.)

    ! Write output
    call nc_write(filename,varname,avar(:,:,:),dim1="x",dim2="y",dim3="time",units=units, long_name=longname)

  end subroutine write_nc_file

  subroutine write_nc_file_annual(avar,res,year,outpathname,fileroot,varname,units,longname, &
                                   institution,contact_name,contact_email)
    ! Write annual-mean variable to netcdf file — one time step per year.

    use ncio

    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(1.0D0)

    REAL(dp), INTENT(IN)  :: avar(:,:,:)    ! shape (nx, ny, 1)
    REAL(dp), INTENT(IN)  :: res
    INTEGER,  INTENT(IN)  :: year
    character(len=*), INTENT(IN) :: outpathname, fileroot
    character(len=*), INTENT(IN) :: varname, units, longname
    character(len=*), INTENT(IN) :: institution, contact_name, contact_email

    character(len=256) :: cyear, filename
    character(len=40)  :: time_units
    double precision :: days

    write (cyear, "(I0.4)") year
    ! Register annual mean at mid-year (day 180 on 360-day calendar)
    days = (year - dble(time_ref_year)) * 360.0d0 + 180.0d0
    write (time_units, "(A,I0.4,A)") "days since ", time_ref_year, "-01-01 00:00:00"
    filename = trim(outpathname) // "/" // trim(varname) &
         // "_" // trim(fileroot) // "_" // trim(cyear) // ".nc"
    write(*,*) "### output file: ", trim(filename)

    call nc_create(filename, overwrite=.TRUE., netcdf4=.TRUE.)
    call nc_write_attr(filename, "title",            "NORCE-PDD output")
    call nc_write_attr(filename, "institution",       trim(institution))
    call nc_write_attr(filename, "contact_name",      trim(contact_name))
    call nc_write_attr(filename, "contact_email",     trim(contact_email))
    call nc_write_attr(filename, "grid_information",  "ISMIP6 grid, epsg:3413")
    call nc_write_attr(filename, "output_frequency",  "annual")

    call nc_write_dim(filename, "x", x=-720000.d0,   dx=res, nx=size(avar,1), units="m")
    call nc_write_dim(filename, "y", x=-3450000.d0,  dx=res, nx=size(avar,2), units="m")
    call nc_write_dim(filename, "time", x=days, dx=360.d0, nx=1, &
         units=trim(time_units), calendar="360_day", unlimited=.TRUE.)

    call nc_write(filename, varname, avar(:,:,:), dim1="x", dim2="y", dim3="time", &
         units=units, long_name=longname)

  end subroutine write_nc_file_annual

end program gpdd
