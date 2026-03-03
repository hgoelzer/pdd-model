program gpdd
! pdd model for greenland

    use ncio 
    use massbalance_module

    implicit none

    REAL, PARAMETER :: rhoi = 910.0
    REAL, PARAMETER :: rhof = 1000.0

    character(len=256) :: testchar

    character(len=256) :: filename, pathname

    integer :: nx, ny, nt 
    integer :: i, j, t, year
    integer :: ndims
    integer :: ncid

    ! dimensions
    double precision,   allocatable :: x(:), y(:), time(:) 

    ! forcing
    double precision,   allocatable :: t2m(:,:,:)
    double precision,   allocatable :: t2j(:,:,:)
    double precision,   allocatable :: tp(:,:,:)

    ! variables 
    double precision,   allocatable :: acc(:,:)
    double precision,   allocatable :: smb(:,:)

    ! output
    double precision,   allocatable :: acc3(:,:,:)
    double precision,   allocatable :: t2m3(:,:,:)
    double precision,   allocatable :: t2j3(:,:,:)
    double precision,   allocatable :: smb3(:,:,:)
    
    character(len=32),  allocatable :: dimnames(:)
    integer,            allocatable :: dimlens(:)

    pathname = "./data"
    
    ! Define array sizes 
    nx = 1681
    ny = 2881
    nt = 30

    ! Allocate dimensions
    allocate(x(nx),y(ny),time(nt))

    ! Allocate arrays
    allocate(t2m(nx,ny,1))
    allocate(t2j(nx,ny,1))
    allocate(tp(nx,ny,1))

    allocate(smb(nx,ny))
    allocate(acc(nx,ny))

    allocate(smb3(nx,ny,nt))
    allocate(t2m3(nx,ny,nt))
    allocate(t2j3(nx,ny,nt))
    allocate(acc3(nx,ny,nt))


    ! Main loop
    t = 1
    DO WHILE(t <= nt)

       year = t+1990
       write(*,*) "Time counter ", t, year, nt

       ! construct t2m filename
       write(*,*) "construct t2m filename"
       write (filename, "(A17,I0.4,A10)") "t2m_CARRA-yearly-", year, "_e01000.nc"
       write(*,*) "construct t2m filename end"
       filename = trim(pathname) // "/" // trim(filename)
       write(*,*) "### File t2m: ", filename
       !call nc_read(filename, "time",time)
       !write(*,"(a10,100f12.1)") "time: ", time 
       ndims = nc_ndims(filename,"t2m")
       call nc_dims(filename,"t2m",dimnames,dimlens)
       write(*,*) "ndims= ", ndims
       write(*,*) "dimnames= ", dimnames
       write(*,*) "dimlens=  ", dimlens 
       call nc_open(filename, ncid, .FALSE.)
       call nc_read(filename, "t2m",t2m, ncid=ncid)
       call nc_close(ncid)

       ! construct t2j filename
       write (filename, "(A21,I0.4,A10)") "t2m_CARRA-monthly-07-", year, "_e01000.nc"
       filename = trim(pathname) // "/" // trim(filename)
       write(*,*) "### File t2j: ", filename
       !call nc_read(filename, "time",time)
       !write(*,"(a10,100f12.1)") "time: ", time 
       ndims = nc_ndims(filename,"t2m")
       call nc_dims(filename,"t2m",dimnames,dimlens)
       write(*,*) "ndims= ", ndims
       write(*,*) "dimnames= ", dimnames
       write(*,*) "dimlens=  ", dimlens 
       call nc_open(filename, ncid, .FALSE.)
       call nc_read(filename, "t2m",t2j)
       call nc_close(ncid)

       ! construct t2m filename
       write (filename, "(A16,I0.4,A10)") "tp_CARRA-yearly-", year, "_e01000.nc"
       filename = trim(pathname) // "/" // trim(filename)
       write(*,*) "### File tp: ", filename
       !call nc_read(filename, "time",time)
       !write(*,"(a10,100f12.1)") "time: ", time 
       ndims = nc_ndims(filename,"tp")
       call nc_dims(filename,"tp",dimnames,dimlens)
       write(*,*) "ndims= ", ndims
       write(*,*) "dimnames= ", dimnames
       write(*,*) "dimlens=  ", dimlens 
       call nc_open(filename, ncid, .FALSE.)
       call nc_read(filename, "tp",tp)
       call nc_close(ncid)

       ! Precip forcing; convert from kg/m2/yr = mm/yr w.e. to m/yr i.e.
       acc = tp(:,:,1)/1000.*rhof/rhoi

       ! Model call
       call pdd_model_greenland_total(nx, ny, acc, t2m, t2j, smb)

       ! Update output container
       smb3(:,:,t) = smb(:,:)
       acc3(:,:,t) = acc(:,:)
       t2m3(:,:,t) = t2m(:,:,1)
       t2j3(:,:,t) = t2j(:,:,1)

       ! update timer
       t = t+1

    END DO ! end main loop


    ! convert SMB terms back to mm w.e/yr
    smb3(:,:,:) = smb3(:,:,:) * rhoi 
    acc3(:,:,:) = acc3(:,:,:) * rhoi 

    ! Writing output file 
    filename = "smb_gpdd.nc"

    ! Create the netcdf file, write global attributes
    call nc_create(filename,overwrite=.TRUE.,netcdf4=.TRUE.)
    call nc_write_attr(filename,"title","CARRA_PDD")
    call nc_write_attr(filename,"institution", "NORCE")

    ! Write the dimensions (x, y, time), defined inline
    call nc_write_dim(filename,"x",x=0.d0,dx=2500.d0,nx=nx,units="m")
    call nc_write_dim(filename,"y",x=0.d0,dx=2500.d0,nx=ny,units="m")
    call nc_write_dim(filename,"time",x=1991.d0,dx=1.d0,nx=nt, &
                      units="years",calendar="360_day", unlimited=.TRUE.)
    
    ! Write time dependent output
    call nc_write(filename,"smb",smb3(:,:,:),dim1="x",dim2="y",dim3="time")
    call nc_write(filename,"acc",acc3(:,:,:),dim1="x",dim2="y",dim3="time")
    !call nc_write(filename,"t2m",t2m3(:,:,:),dim1="x",dim2="y",dim3="time")
    !call nc_write(filename,"t2j",t2j3(:,:,:),dim1="x",dim2="y",dim3="time")


    ! Clean up
    deallocate(x,y,time)

    ! 2D
    deallocate(acc)
    deallocate(smb)
    deallocate(t2m)
    deallocate(t2j)

    ! 3D
    deallocate(smb3)
    deallocate(acc3)

end program 
