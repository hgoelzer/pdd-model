program gpdd
! pdd model for greenland with monthly input and annual output

    use ncio 
    use massbalance_module

    implicit none

    REAL, PARAMETER :: rhoi = 910.0
    REAL, PARAMETER :: rhof = 1000.0

    character(len=256) :: testchar

    character(len=256) :: filename
    character(len=256) :: inpathname, outpathname

    integer :: nx, ny, nt 
    integer :: i, j, t, year, m
    integer :: ndims
    integer :: ncid
    double precision :: dyear

    ! dimensions
    double precision,   allocatable :: x(:), y(:), time(:) 

    ! forcing
    double precision,   allocatable :: t2m_in(:,:,:)
    double precision,   allocatable :: tp_in(:,:,:)

    double precision,   allocatable :: t2m(:,:,:) ! monthly
    double precision,   allocatable :: tp(:,:,:) ! monthly

    ! variables 
    double precision,   allocatable :: tma(:,:) ! annual
    double precision,   allocatable :: tpa(:,:) ! annual
    double precision,   allocatable :: smb(:,:) ! annual

    double precision,   allocatable :: snow(:,:)
    double precision,   allocatable :: rain(:,:)
    double precision,   allocatable :: sir(:,:)
    double precision,   allocatable :: abl(:,:)
    double precision,   allocatable :: pdd(:,:)
    double precision,   allocatable :: rfr(:,:)

    ! output
    double precision,   allocatable :: tma3(:,:,:) ! annual time series
    double precision,   allocatable :: tpa3(:,:,:) ! annual time series
    double precision,   allocatable :: smb3(:,:,:) ! annual time series

    double precision,   allocatable :: snow3(:,:,:)
    double precision,   allocatable :: rain3(:,:,:)
    double precision,   allocatable :: sir3(:,:,:)
    double precision,   allocatable :: abl3(:,:,:)
    double precision,   allocatable :: pdd3(:,:,:)
    double precision,   allocatable :: rfr3(:,:,:)
    
    character(len=32),  allocatable :: dimnames(:)
    integer,            allocatable :: dimlens(:)

    inpathname = "./data"
    outpathname = "./output"
    
    ! Define array sizes 
    nx = 1681
    ny = 2881
    nt = 1 ! 1991 - 

    ! Allocate dimensions
    allocate(x(nx),y(ny),time(nt))

    ! Allocate arrays
    allocate(t2m_in(nx,ny,1))
    allocate(tp_in(nx,ny,1))

    allocate(t2m(nx,ny,12))
    allocate(tp(nx,ny,12))

    allocate(tma(nx,ny))
    allocate(tpa(nx,ny))
    allocate(smb(nx,ny))

    allocate(snow(nx,ny))
    allocate(rain(nx,ny))
    allocate(sir(nx,ny))
    allocate(abl(nx,ny))
    allocate(pdd(nx,ny))
    allocate(rfr(nx,ny))

    allocate(tma3(nx,ny,nt))
    allocate(tpa3(nx,ny,nt))
    allocate(smb3(nx,ny,nt))

    allocate(snow3(nx,ny,nt))
    allocate(rain3(nx,ny,nt))
    allocate(sir3(nx,ny,nt))
    allocate(abl3(nx,ny,nt))
    allocate(pdd3(nx,ny,nt))
    allocate(rfr3(nx,ny,nt))

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

       ! Annual diagnostics
       tma = SUM(t2m(:,:,:),DIM=3)/12
       tpa = SUM(tp(:,:,:),DIM=3)

       ! Precip forcing; convert from kg/m2/yr = mm/yr w.e. to m/yr i.e.
       tp(:,:,:) = tp(:,:,:)/1000.*rhof/rhoi

       ! Model call
       call pdd_model_greenland_total_monthly(nx, ny, 0.003d0, 0.008d0, 4.5d0, 1.0d0, tp, t2m, smb, snow, rain, sir, abl, pdd, rfr)

       !! Update output container
       tma3(:,:,1) = tma(:,:)
       tpa3(:,:,1) = tpa(:,:)
       smb3(:,:,1) = smb(:,:)
       
       snow3(:,:,1) = snow(:,:)
       rain3(:,:,1) = rain(:,:)
       sir3(:,:,1) = sir(:,:)
       abl3(:,:,1) = abl(:,:)
       pdd3(:,:,1) = pdd(:,:)
       rfr3(:,:,1) = rfr(:,:)

       ! Write out annual files
       ! convert SMB terms back to mm w.e/yr
       smb3(:,:,:) = smb3(:,:,:) * rhoi 
       tpa3(:,:,:) = tpa3(:,:,:) * rhoi 

       ! Writing output file
       write (filename, "(A42,I0.4,A3)") "GIS_NORCE-CISM-PDD_CESM2_ssp370_r11i1p1f1_", year, ".nc"
       filename = trim(outpathname) // "/" // trim(filename)
       dyear = year
       
       ! Create the netcdf file, write global attributes
       call nc_create(filename,overwrite=.TRUE.,netcdf4=.TRUE.)
       call nc_write_attr(filename,"title","NORCE-CISM-PDD")
       call nc_write_attr(filename,"institution", "NORCE")

       ! Write the dimensions (x, y, time), defined inline
       call nc_write_dim(filename,"x",x=0.d0,dx=2500.d0,nx=nx,units="m")
       call nc_write_dim(filename,"y",x=0.d0,dx=2500.d0,nx=ny,units="m")
       call nc_write_dim(filename,"time",x=dyear,dx=1.d0,nx=1, &
            units="years",calendar="360_day", unlimited=.TRUE.)

       ! Write time dependent output
       call nc_write(filename,"smb",smb3(:,:,1),dim1="x",dim2="y",dim3="time")
       ! forcing
       call nc_write(filename,"tma",tma3(:,:,1),dim1="x",dim2="y",dim3="time")
       call nc_write(filename,"tpa",tpa3(:,:,1),dim1="x",dim2="y",dim3="time")
       ! components
       call nc_write(filename,"rfr",rfr3(:,:,1),dim1="x",dim2="y",dim3="time")
       call nc_write(filename,"snow",snow3(:,:,1),dim1="x",dim2="y",dim3="time")
       call nc_write(filename,"rain",rain3(:,:,1),dim1="x",dim2="y",dim3="time")
       call nc_write(filename,"abl",abl3(:,:,1),dim1="x",dim2="y",dim3="time")
       call nc_write(filename,"pdd",pdd3(:,:,1),dim1="x",dim2="y",dim3="time")
       call nc_write(filename,"sir",sir3(:,:,1),dim1="x",dim2="y",dim3="time")

       
       ! update timer
       t = t+1

    END DO ! end main loop


    ! Clean up
    deallocate(x,y,time)

    deallocate(t2m_in)
    deallocate(tp_in)

    deallocate(t2m)
    deallocate(tp)

    deallocate(tma3)
    deallocate(tpa3)
    deallocate(smb3)

    deallocate(snow3)
    deallocate(rain3)
    deallocate(abl3)
    deallocate(pdd3)
    deallocate(sir3)
    deallocate(rfr3)


end program 
