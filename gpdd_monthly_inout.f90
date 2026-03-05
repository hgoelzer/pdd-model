program gpdd
  ! pdd model for greenland with monthly input and output

  use ncio 
  use massbalance_module

  implicit none

  REAL, PARAMETER :: rhoi = 910.0
  REAL, PARAMETER :: rhof = 1000.0

  character(len=256) :: filename
  character(len=256) :: inpathname, outpathname

  integer :: nx, ny, nt 
  integer :: i, j, t, year, m
  integer :: ndims
  integer :: ncid

  ! dimensions
  double precision,   allocatable :: x(:), y(:)

  ! forcing
  double precision,   allocatable :: t2m_in(:,:,:)
  double precision,   allocatable :: tp_in(:,:,:)

  ! variables, all monthly
  double precision,   allocatable :: t2m(:,:,:)
  double precision,   allocatable :: tp(:,:,:)
  double precision,   allocatable :: smb(:,:,:)
  double precision,   allocatable :: snow(:,:,:)
  double precision,   allocatable :: rain(:,:,:)
  double precision,   allocatable :: sir(:,:,:)
  double precision,   allocatable :: abl(:,:,:)
  double precision,   allocatable :: pdd(:,:,:)
  double precision,   allocatable :: rfr(:,:,:)

  character(len=32),  allocatable :: dimnames(:)
  integer,            allocatable :: dimlens(:)

  inpathname = "./data"
  outpathname = "./output"

  ! Define array sizes 
  nx = 1681
  ny = 2881
  nt = 2 ! 1991 - 

  ! Allocate dimensions
  allocate(x(nx),y(ny))

  ! Allocate arrays
  allocate(t2m_in(nx,ny,1))
  allocate(tp_in(nx,ny,1))

  allocate(t2m(nx,ny,12))
  allocate(tp(nx,ny,12))

  allocate(smb(nx,ny,12))
  allocate(snow(nx,ny,12))
  allocate(rain(nx,ny,12))
  allocate(sir(nx,ny,12))
  allocate(abl(nx,ny,12))
  allocate(pdd(nx,ny,12))
  allocate(rfr(nx,ny,12))

  ! Main year loop
  t = 1
  DO WHILE(t <= nt)

     year = t+1990
     write(*,*) "Time counter ", t, year, nt

     ! month loop
     m = 1
     DO WHILE(m <= 12)

        ! construct t2m filename
        write (filename, "(A18,I0.2,A1,I0.4,A10)") "t2m_CARRA-monthly-", m, "-", year, "_e01000.nc"
        filename = trim(inpathname) // "/" // trim(filename)
        write(*,*) 
        write(*,*) "### File t2m: ", filename
        ndims = nc_ndims(filename,"t2m")
        call nc_dims(filename,"t2m",dimnames,dimlens)
        write(*,*) "ndims= ", ndims
        !write(*,*) "dimnames= ", dimnames
        write(*,*) "dimlens=  ", dimlens 
        call nc_open(filename, ncid, .FALSE.)
        call nc_read(filename, "t2m",t2m_in)
        call nc_close(ncid)
        t2m(:,:,m) = t2m_in(:,:,1)

        ! construct tp filename
        write (filename, "(A17,I0.2,A1,I0.4,A10)") "tp_CARRA-monthly-", m, "-", year, "_e01000.nc"
        filename = trim(inpathname) // "/" // trim(filename)
        write(*,*) 
        write(*,*) "### File tp: ", filename
        ndims = nc_ndims(filename,"tp")
        call nc_dims(filename,"tp",dimnames,dimlens)
        write(*,*) "ndims= ", ndims
        !write(*,*) "dimnames= ", dimnames
        write(*,*) "dimlens=  ", dimlens 
        call nc_open(filename, ncid, .FALSE.)
        call nc_read(filename, "tp",tp_in)
        call nc_close(ncid)
        tp(:,:,m) = tp_in(:,:,1)

        ! update timer
        m = m+1

     END DO ! month loop


     ! Precip forcing; convert from kg/m2/yr = mm/yr w.e. to m/yr i.e.
     tp(:,:,:) = tp(:,:,:)/1000.*rhof/rhoi

     ! Model call
     call pdd_model_greenland_total_monthly_inout(nx, ny, tp, t2m, smb, snow, rain, sir, abl, pdd, rfr)

     ! convert SMB terms to kg m-2 s-1 = mm w.e. / s
     smb(:,:,:) = smb(:,:,:) * rhoi / 31556926
     tp(:,:,:) = tp(:,:,:) * rhoi / 31556926
     rain(:,:,:) = rain(:,:,:) * rhoi / 31556926
     abl(:,:,:) = abl(:,:,:) * rhoi / 31556926
     sir(:,:,:) = sir(:,:,:) * rhoi / 31556926

     ! Write files
     call write_nc_file(t2m, year, outpathname, "tas", "K", "air_temperature")
     call write_nc_file(tp, year, outpathname, "pr", "kg m-2 s-1", "precipitation_flux")
     call write_nc_file(smb, year, outpathname, "acabf", "kg m-2 s-1", "land_ice_surface_specific_mass_balance_flux")
     call write_nc_file(abl, year, outpathname, "mrro", "kg m-2 s-1", "runoff_flux")
     call write_nc_file(rain, year, outpathname, "prra", "kg m-2 s-1", "rainfall_flux")
     call write_nc_file(sir, year, outpathname, "snicefreez", "kg m-2 s-1", "surface_snow_and_ice_refreezing_flux")
     
     ! update timer
     t = t+1

  END DO ! end main loop


  ! Clean up
  deallocate(x,y)

  deallocate(t2m_in)
  deallocate(tp_in)

  deallocate(t2m)
  deallocate(tp)

  deallocate(smb)
  deallocate(snow)
  deallocate(rain)
  deallocate(sir)
  deallocate(abl)
  deallocate(pdd)
  deallocate(rfr)


contains

  subroutine write_nc_file(avar,year,outpathname,varname,units,longname)
    ! Write variable to netcdf file

    use ncio

    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(1.0D0)  ! Kind of double precision numbers.

    REAL(dp), INTENT(IN)  :: avar(:,:,:)
    INTEGER, INTENT(IN)  :: year
    character(len=*), INTENT(IN) :: outpathname
    character(len=*), INTENT(IN) :: varname
    character(len=*), INTENT(IN) :: units
    character(len=*), INTENT(IN) :: longname

    ! Internal variables
    character(len=256) :: fileroot, cyear
    character(len=256) :: filename
    double precision :: days

    ! Construct output filename
    write (fileroot, "(A41)") "GIS_NORCE-CISM-PDD_CESM2_ssp370_r11i1p1f1"
    write (cyear, "(I0.4)") year
    !write(*,*) "### root: ", fileroot
    !write(*,*) "### cyear: ", cyear
    filename = trim(outpathname) // "/" // trim(varname) &
         // "_" // trim(fileroot) // "_" // trim(cyear) // ".nc" 
    write(*,*) "### output file: ", filename

    ! Create the netcdf file, write global attributes
    call nc_create(filename,overwrite=.TRUE.,netcdf4=.TRUE.)
    call nc_write_attr(filename,"title","NORCE-CISM-PDD output")
    call nc_write_attr(filename,"institution", "NORCE Research")
    call nc_write_attr(filename,"contact_name", "Heiko Goelzer")
    call nc_write_attr(filename,"contact_email", "heig@norceresearch.no")
    call nc_write_attr(filename,"grid_information", "ISMIP6 grid 1km, epsg:3413")

    ! Write the dimensions (x, y, time), defined inline
    call nc_write_dim(filename,"x",x=-720000.d0,dx=1000.d0,nx=nx,units="m")
    call nc_write_dim(filename,"y",x=-3450000.d0,dx=1000.d0,nx=ny,units="m")
    call nc_write_dim(filename,"time",x=days,dx=30.d0,nx=12, &
         units="days since 1900-01-01 00:00:00",calendar="360_day", unlimited=.TRUE.)
    days = (year-1900.0)*360.0+15 ! register monthly means to mid month on 360-day calendar

    ! Write output
    call nc_write(filename,varname,avar(:,:,:),dim1="x",dim2="y",dim3="time",units=units, long_name=longname)

  end subroutine write_nc_file

end program gpdd
