program gpdd
! pdd model for greenland

    use ncio 
    use massbalance_module

    implicit none

    REAL, PARAMETER :: rhoi = 910.0
    REAL, PARAMETER :: rhof = 1000.0

    character(len=256) :: testchar

    character(len=:), allocatable :: filename, pathname

    integer :: nx, ny, nt 
    integer :: i, j, t
    integer :: ndims

    ! dimensions
    double precision,   allocatable :: x(:), y(:), time(:) 

    ! input 
    double precision,   allocatable :: lat(:,:)
    double precision,   allocatable :: orog(:,:)

    ! forcing
    real,               allocatable :: t2m(:,:,:)
    real,               allocatable :: tp(:,:,:)

    ! variables 
    double precision,   allocatable :: acc(:,:)
    double precision,   allocatable :: T_anomaly(:,:)
    double precision,   allocatable :: smb(:,:)
    double precision,   allocatable :: T_avg(:,:)

    ! output
    double precision,   allocatable :: acc3(:,:,:)
    double precision,   allocatable :: T_anomaly3(:,:,:)
    double precision,   allocatable :: smb3(:,:,:)
    
    character(len=32),  allocatable :: dimnames(:)
    integer,            allocatable :: dimlens(:)

    pathname = "./data"
    
    ! Reading orog as base for dimensions 
    filename = pathname // "/" // "orog_e16000.nc"
    write(*,*)
    write(*,*) "### File orog: ", filename
    ndims = nc_ndims(filename,"orog")
    call nc_dims(filename,"orog",dimnames,dimlens)
    write(*,*) "ndims= ", ndims
    write(*,*) "dimnames= ", dimnames
    write(*,*) "dimlens=  ", dimlens 
    write(*,*) dimlens(1) 
    write(*,*) dimlens(2) 

    ! Define array sizes 
    nx = dimlens(1)
    ny = dimlens(2)
    nt = 10

    ! Allocate dimensions
    allocate(x(nx),y(ny),time(nt))

    ! Allocate arrays
    allocate(lat(nx,ny))
    allocate(orog(nx,ny))

    allocate(t2m(nx,ny,nt))
    !allocate(tp(nx,ny,nt))
    allocate(tp(nx,ny,1))

    allocate(T_anomaly(nx,ny))
    allocate(smb(nx,ny))
    allocate(acc(nx,ny))
    allocate(T_avg(nx,ny))

    allocate(T_anomaly3(nx,ny,nt))
    allocate(smb3(nx,ny,nt))
    allocate(acc3(nx,ny,nt))


    !! Initialisation
    call nc_read(filename, "orog",orog)


    ! Reading lat
    filename = pathname // "/" // "lat_lon_16km.nc"
    write(*,*)
    write(*,*) "### File lat: ", filename
    ndims = nc_ndims(filename,"lat")
    call nc_dims(filename,"lat",dimnames,dimlens)
    write(*,*) "ndims= ", ndims
    write(*,*) "dimnames= ", dimnames
    write(*,*) "dimlens=  ", dimlens 
    call nc_read(filename, "lat",lat)

    ! Read forcing files 
    filename = pathname // "/" // "t2m_CARRA-yearly_e16000.nc"
    write(*,*)
    write(*,*) "### File t2m: ", filename
    !call nc_read_attr(filename, "institution", testchar)
    !write(*,*) "Institution: ", trim(testchar)
    call nc_read(filename, "time",time)
    write(*,"(a10,100f12.1)") "time: ", time 
    ndims = nc_ndims(filename,"t2m")
    call nc_dims(filename,"t2m",dimnames,dimlens)
    write(*,*) "ndims= ", ndims
    write(*,*) "dimnames= ", dimnames
    write(*,*) "dimlens=  ", dimlens 
    call nc_read(filename, "t2m",t2m)
    !write(*,"(a10,100f12.1)") "t2m: ", t2m 

    ! Reading precip
    filename = pathname // "/" // "tp_CARRA-yearly-1991.nc_e16000.nc"
    write(*,*)
    write(*,*) "### File tp: ", filename
    ndims = nc_ndims(filename,"tp")
    call nc_dims(filename,"tp",dimnames,dimlens)
    write(*,*) "ndims= ", ndims
    write(*,*) "dimnames= ", dimnames
    write(*,*) "dimlens=  ", dimlens 
    call nc_read(filename, "tp",tp)


    ! Calculate long-term average t2m 
    T_avg = SUM(t2m(:,:,:),DIM=3)/nt


    ! Main loop
    t = 1
    DO WHILE(t <= nt)

       write(*,*) "Time counter ", t, nt

       ! construct filename
       write (filename, "(A8,I0.4,A3)") "forcing-", t, ".nc"
       write(*,*) trim(filename)

       ! temperature forcing
       T_anomaly = t2m(:,:,t)-T_avg(:,:)
       
       ! Precip forcing; convert from kg/m2/yr = mm/yr w.e. to m/yr i.e.
       acc = tp(:,:,1)/1000.*rhof/rhoi
       !acc = tp(:,:,t)/1000.*rhof/rhoi

       ! Model call
       call massbalance_pdd_model_greenland(nx,ny,lat,orog,acc,T_anomaly,smb)

       ! Update output container
       T_anomaly3(:,:,t) = T_anomaly(:,:)
       smb3(:,:,t) = smb(:,:)
       acc3(:,:,t) = acc(:,:)

       ! update timer
       t = t+1

    END DO ! end main loop


    ! Writing output file 
    filename = "smb_gpdd.nc"

    ! Create the netcdf file, write global attributes
    call nc_create(filename,overwrite=.TRUE.,netcdf4=.TRUE.)
    call nc_write_attr(filename,"title","CARRA_PDD")
    call nc_write_attr(filename,"institution", "NORCE")

    ! Write the dimensions (x, y, time), defined inline
    call nc_write_dim(filename,"x",x=0.d0,dx=16000.d0,nx=nx,units="m")
    call nc_write_dim(filename,"y",x=0.d0,dx=16000.d0,nx=ny,units="m")
    call nc_write_dim(filename,"time",x=1991.d0,dx=1.d0,nx=nt, &
                      units="years",calendar="360_day", unlimited=.TRUE.)
    
    ! Write output fixed
    call nc_write(filename,"lat",lat(:,:),dim1="x",dim2="y")
    call nc_write(filename,"orog",orog(:,:),dim1="x",dim2="y")
    call nc_write(filename,"tavg",T_avg(:,:),dim1="x",dim2="y")

    ! Write time dependent output
    call nc_write(filename,"smb",smb3(:,:,:),dim1="x",dim2="y",dim3="time")
    call nc_write(filename,"acc",acc3(:,:,:),dim1="x",dim2="y",dim3="time")
    call nc_write(filename,"dT",T_anomaly3(:,:,:),dim1="x",dim2="y",dim3="time")


    ! Clean up
    deallocate(x,y,time)

    ! 2D
    deallocate(lat)
    deallocate(orog)
    deallocate(acc)
    deallocate(T_anomaly)
    deallocate(smb)
    deallocate(t2m)
    deallocate(T_avg)

    ! 3D
    deallocate(T_anomaly3)
    deallocate(smb3)
    deallocate(acc3)

end program 
