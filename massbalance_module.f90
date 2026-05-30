! File name: massbalance_module.f90

MODULE massbalance_module

  IMPLICIT NONE
  INTEGER,  PARAMETER :: dp     = KIND(1.0D0)
  REAL(dp), PARAMETER :: pi     = 2.0_dp * ACOS(0.0_dp)
  REAL(dp), PARAMETER :: valmax = 6.0_dp
  INTEGER,  PARAMETER :: nintx  = 1200

  REAL(dp), SAVE :: taberf(-nintx:nintx)
  REAL(dp), SAVE :: tabepdd(-nintx:nintx)
  LOGICAL,  SAVE :: lut_initialized = .FALSE.

CONTAINS

  SUBROUTINE init_pdd_lut()
    REAL(dp) :: deltax, sq2pi, fac1, fdx, xi, xj, yi, yj, help
    INTEGER  :: i
    IF (lut_initialized) RETURN
    taberf(0)  = 0.0
    deltax = valmax / nintx
    sq2pi  = (2*pi)**(0.5)
    fac1   = deltax / (2*sq2pi)
    tabepdd(0) = 1./sq2pi
    xj = 0.
    yj = 1.
    DO i = 1, nintx
      xi  = xj
      yi  = yj
      xj  = xj + deltax
      yj  = exp(-0.5*xj*xj)
      fdx = (yi + yj) * fac1
      taberf(i)   = taberf(i-1) + fdx
      taberf(-i)  = -taberf(i)
      help        = yj / sq2pi + xj * taberf(i)
      tabepdd(i)  = help + xj * 0.5
      tabepdd(-i) = help - xj * 0.5
    END DO
    lut_initialized = .TRUE.
  END SUBROUTINE init_pdd_lut


  SUBROUTINE pdd_model_greenland_total_monthly_inout(nx, ny, ddfactorsnow, ddfactorice, sigma, rainlimit, tp, t2m, smb, snow, rain, sir, abl, pdd, rfr)
    ! Positive degree day model, uddated by Heiko Goelzer, Mar 2026
    ! Implements Greenland pdd model of Huybrechts and De Wolde 1999
    ! Forcing with total monthly fields, output monthly data

    IMPLICIT NONE


    ! --------------------------------------------------------------------------
    ! Declaration of global variables
    ! --------------------------------------------------------------------------

    ! Input variables: 
    INTEGER, INTENT(IN)  :: nx, ny ! grid size

    REAL(dp), INTENT(IN)  :: tp(nx,ny,12)  ! total monthly precip (m/yr)
    REAL(dp), INTENT(IN)  :: t2m(nx,ny,12)  ! monthly mean 2m temperature (K)

    ! Output variables: 
    REAL(dp), INTENT(OUT) :: smb(nx,ny,12)     ! surface mass balance (m/yr)
    REAL(dp), INTENT(OUT) :: snow(nx,ny,12)
    REAL(dp), INTENT(OUT) :: rain(nx,ny,12)
    REAL(dp), INTENT(OUT) :: sir(nx,ny,12)
    REAL(dp), INTENT(OUT) :: abl(nx,ny,12)           ! runoff (m/yr)
    REAL(dp), INTENT(OUT) :: pdd(nx,ny,12)
    REAL(dp), INTENT(OUT) :: rfr(nx,ny,12)

    REAL(dp), INTENT(IN)                :: ddfactorsnow
    REAL(dp), INTENT(IN)                :: ddfactorice
    REAL(dp), INTENT(IN)                :: sigma
    REAL(dp), INTENT(IN)                :: rainlimit

    ! Local variables
    REAL(dp), allocatable               :: tm(:,:,:)

    REAL(dp), PARAMETER                 :: pmax = 0.3 ! See update in Janssens and Huybrechts 2000
 

    ! Allocate arrays
    allocate(tm(nx,ny,12))

    ! Monthly temperature (C)
    tm = t2m - 273.15

    ! Determine number of positive degree days per year and rain fraction
    !call calculate_pdd_monthly_inout(nx, ny, tm, pdd, rfr)
    ! With parameterised seasonal cycle
    call calculate_pdd_monthly_inout_taj(nx, ny, sigma, rainlimit, tm, pdd, rfr)

    ! Distinguish rain and snow according to rain fraction
    rain = tp * rfr
    snow = tp - rain

    call melt_cascade_3d(nx, ny, ddfactorsnow, ddfactorice, pmax, snow, pdd, tp, sir, abl)

    ! no melt where insufficient energy
    WHERE (pdd.LE.0)
     sir = 0
     abl = 0
    END WHERE

    smb = tp - abl - rain ! = snow - abl

  END SUBROUTINE pdd_model_greenland_total_monthly_inout



  SUBROUTINE pdd_model_greenland_total_monthly(nx, ny, ddfactorsnow, ddfactorice, sigma, rainlimit, tp, t2m, smb, snow, rain, sir, abl, pdd, rfr)
    ! Positive degree day model, uddated by Heiko Goelzer, Feb 2022
    ! Implements Greenland pdd model of Huybrechts and De Wolde 1999
    ! Forcing with total monthly fields

    IMPLICIT NONE


    ! --------------------------------------------------------------------------
    ! Declaration of global variables
    ! --------------------------------------------------------------------------

    ! Input variables: 
    INTEGER, INTENT(IN)  :: nx, ny ! grid size

    REAL(dp), INTENT(IN)  :: tp(nx,ny,12)  ! total monthly precip (m/yr)
    REAL(dp), INTENT(IN)  :: t2m(nx,ny,12)  ! monthly mean 2m temperature (deg)

    ! Output variables: 
    REAL(dp), INTENT(OUT) :: smb(nx,ny)     ! surface mass balance (m/yr)
    REAL(dp), INTENT(OUT) :: snow(nx,ny)
    REAL(dp), INTENT(OUT) :: rain(nx,ny)
    REAL(dp), INTENT(OUT) :: sir(nx,ny)
    REAL(dp), INTENT(OUT) :: abl(nx,ny)           ! runoff (m/yr)
    REAL(dp), INTENT(OUT) :: pdd(nx,ny)
    REAL(dp), INTENT(OUT) :: rfr(nx,ny)

    REAL(dp), INTENT(IN)                :: ddfactorsnow
    REAL(dp), INTENT(IN)                :: ddfactorice
    REAL(dp), INTENT(IN)                :: sigma
    REAL(dp), INTENT(IN)                :: rainlimit

    ! Local variables
    REAL(dp), allocatable               :: tm(:,:,:)
    REAL(dp), allocatable               :: tpa(:,:)

    REAL(dp), PARAMETER                 :: pmax = 0.3 ! See update in Janssens and Huybrechts 2000
 

    ! Allocate arrays
    allocate(tm(nx,ny,12))
    allocate(tpa(nx,ny))

    ! Monthly temperature
    tm = t2m - 273.15

    ! Determine number of positive degree days per year and rain fraction
    call calculate_pdd_monthly(nx, ny, sigma, rainlimit, tm, pdd, rfr)

    ! total annual precip
    tpa(:,:) = SUM(tp(:,:,:),dim=3)

    ! Distinguish rain and snow according to rain fraction
    rain = tpa * rfr
    snow = tpa - rain

    call melt_cascade_2d(nx, ny, ddfactorsnow, ddfactorice, pmax, snow, pdd, tpa, sir, abl)

    ! no melt where insufficient energy
    WHERE (pdd.LE.0)
     sir = 0
     abl = 0
    END WHERE

    ! smb = snow - abl
    smb = tpa - abl - rain

  END SUBROUTINE pdd_model_greenland_total_monthly



  SUBROUTINE pdd_model_greenland_total_yearly(nx, ny, ddfactorsnow, ddfactorice, sigma, rainlimit, acc, t2m, t2j, smb, snow, rain, sir, abl, pdd, rfr)
    ! Positive degree day model, uddated by Heiko Goelzer, Feb 2022
    ! Implements Greenland pdd model of Huybrechts and De Wolde 1999
    ! Forcing with total fields, not anomalies

    IMPLICIT NONE



    ! --------------------------------------------------------------------------
    ! Declaration of global variables
    ! --------------------------------------------------------------------------

    ! Input variables: 
    INTEGER, INTENT(IN)  :: nx, ny ! grid size

    REAL(dp), INTENT(IN)  :: acc(nx,ny)    ! total yearly accumulation (m/yr)
    REAL(dp), INTENT(IN)  :: t2m(nx,ny)    ! annual mean 2m temperature (deg)
    REAL(dp), INTENT(IN)  :: t2j(nx,ny)    ! july 2m temperature (deg)

    ! Output variables: 
    REAL(dp), INTENT(OUT) :: smb(nx,ny)        ! surface mass balance (m/yr)
    REAL(dp), INTENT(OUT) :: snow(nx,ny)
    REAL(dp), INTENT(OUT) :: rain(nx,ny)
    REAL(dp), INTENT(OUT) :: sir(nx,ny)
    REAL(dp), INTENT(OUT) :: abl(nx,ny)           ! runoff (m/yr)
    REAL(dp), INTENT(OUT) :: pdd(nx,ny)
    REAL(dp), INTENT(OUT) :: rfr(nx,ny)

    REAL(dp), INTENT(IN)                :: ddfactorsnow
    REAL(dp), INTENT(IN)                :: ddfactorice
    REAL(dp), INTENT(IN)                :: sigma
    REAL(dp), INTENT(IN)                :: rainlimit

    ! Local variables
    REAL(dp), allocatable               :: tma(:,:)
    REAL(dp), allocatable               :: tmj(:,:)

    REAL(dp), PARAMETER                 :: pmax = 0.3 ! See update in Janssens and Huybrechts 2000


    ! Allocate arrays
    allocate(tma(nx,ny))
    allocate(tmj(nx,ny))

    ! Global temperature perturbation
    tma = t2m - 273.15

    ! Summer temperature
    tmj = t2j - 273.15

    ! Determine number of positive degree days per year and rain fraction
    call calculate_pdd_yearly(nx, ny, sigma, rainlimit, tma, tmj, pdd, rfr)

    ! Distinguish rain and snow according to rain fraction
    rain = acc * rfr
    snow = acc - rain

    call melt_cascade_2d(nx, ny, ddfactorsnow, ddfactorice, pmax, snow, pdd, acc, sir, abl)

    ! no melt where insufficient energy
    WHERE (pdd.LE.0)
     sir = 0
     abl = 0
    END WHERE

    ! smb = snow - abl
    smb = acc - abl - rain

    deallocate(tma)
    deallocate(tmj)

  END SUBROUTINE pdd_model_greenland_total_yearly



  SUBROUTINE massbalance_pdd_model_greenland(nx, ny, ddfactorsnow, ddfactorice, sigma, rainlimit, lat, Hs, acc_PD, T_anomaly, smb)
    ! Positive degree day model, added by Heiko Goelzer, Jan 2017
    ! Implements Greenland pdd model of Huybrechts and De Wolde 1999

    IMPLICIT NONE



    ! --------------------------------------------------------------------------
    ! Declaration of global variables
    ! --------------------------------------------------------------------------

    ! Input variables: 
    INTEGER, INTENT(IN)  :: nx, ny ! grid size

    REAL(dp), INTENT(IN)  :: lat(nx,ny)            ! latitude (deg+)
    REAL(dp), INTENT(IN)  :: Hs(nx,ny)             ! surface elevation (m)
    ! Hs_PD not in use; compare with Antarctic model
!    REAL(dp), DIMENSION(   nx,ny), INTENT(IN)  :: Hs_PD         ! reference surface elevation (m)

    REAL(dp), INTENT(IN)  :: acc_PD(nx,ny)    ! reference accumulation (m/yr)
    REAL(dp), INTENT(IN)  :: T_anomaly(nx,ny) ! temperature anomaly w.r.t. PD (deg)

    ! Output variables: 
    REAL(dp), INTENT(OUT) :: smb(nx,ny)        ! surface mass balance (m/yr)

    ! Local variables
    REAL(dp), allocatable               :: tma(:,:)
    REAL(dp), allocatable               :: snow(:,:)
    REAL(dp), allocatable               :: rain(:,:)
    REAL(dp), allocatable               :: sir(:,:)
    REAL(dp), allocatable               :: acc(:,:)
    REAL(dp), allocatable               :: abl(:,:)           ! runoff (m/yr)
    REAL(dp), allocatable               :: tmj(:,:)
    REAL(dp), allocatable               :: pdd(:,:)
    REAL(dp), allocatable               :: rfr(:,:)
    REAL(dp), allocatable               :: s_prec(:,:)
    REAL(dp), allocatable               :: h_inv(:,:)

    REAL(dp), INTENT(IN)                :: ddfactorsnow
    REAL(dp), INTENT(IN)                :: ddfactorice
    REAL(dp), INTENT(IN)                :: sigma
    REAL(dp), INTENT(IN)                :: rainlimit

    REAL(dp), PARAMETER                 :: pmax = 0.3 ! See update in Janssens and Huybrechts 2000


    ! Allocate arrays
    allocate(tma(nx,ny))
    allocate(snow(nx,ny))
    allocate(rain(nx,ny))
    allocate(sir(nx,ny))
    allocate(acc(nx,ny))
    allocate(abl(nx,ny))
    allocate(tmj(nx,ny))
    allocate(pdd(nx,ny))
    allocate(rfr(nx,ny))
    allocate(s_prec(nx,ny))
    allocate(h_inv(nx,ny))

    ! Determination of annual temperature 
    h_inv=20.*(abs(lat)-65.) 
    WHERE (Hs<h_inv)
     tma=49.13-0.007992*h_inv-0.7576*abs(lat)
    ELSEWHERE
     tma=49.13-0.007992*Hs-0.7576*abs(lat)
    END WHERE

    ! Global temperature perturbation
    tma = tma + T_anomaly

    ! Calculate summer temperatures
    tmj=30.78-0.006277*hs-0.3262*abs(lat)+T_anomaly  

    ! Calculate precipitation
    ! See equation (C13) and (C14) in Huybrechts and De Wolde 1999
    WHERE (T_anomaly>=0.) 
     s_prec=0.05
    ELSEWHERE (T_anomaly<=0. .AND. T_anomaly>=-10.)
     s_prec=0.05 - 0.005*T_anomaly
    ELSEWHERE
     s_prec=0.1
    END WHERE
    acc = acc_PD * (1.+s_prec)**T_anomaly

    ! Determine number of positive degree days per year and rain fraction
    call calculate_pdd_yearly(nx, ny, sigma, rainlimit, tma, tmj, pdd, rfr)

    ! Distinguish rain and snow according to rain fraction
    rain = acc * rfr
    snow = acc - rain

    call melt_cascade_2d(nx, ny, ddfactorsnow, ddfactorice, pmax, snow, pdd, acc, sir, abl)

    ! no melt where insufficient energy
    WHERE (pdd.LE.0)
     sir = 0
     abl = 0
    END WHERE

    ! smb = snow - abl
    smb = acc - abl - rain

  END SUBROUTINE massbalance_pdd_model_greenland



  SUBROUTINE calculate_pdd_yearly(nx, ny, sigma, rainlimit, tma, tmj, pdd, rfr)
    ! Positive degree day model, added by Heiko Goelzer, Jan 2017
    ! PDD model from Huybrechts and De Wolde 1999
    ! Seasonal cycle is parametised using annual mean and july temperature

    IMPLICIT NONE



    ! Input variables: 
    INTEGER, INTENT(IN)  :: nx, ny ! grid size

    REAL(dp), INTENT(IN)  :: tma(nx,ny)
    REAL(dp), INTENT(IN)  :: tmj(nx,ny)

    ! Output variables: 
    REAL(dp), INTENT(OUT)  :: pdd(nx,ny)
    REAL(dp), INTENT(OUT)  :: rfr(nx,ny)

    ! Local variables:
    REAL(dp), INTENT(IN)      :: sigma
    REAL(dp), INTENT(IN)      :: rainlimit

    INTEGER                   :: i, j, k
    REAL(dp), allocatable     :: pdd12(:,:,:)
    REAL(dp), allocatable     :: rfr12(:,:,:)
    REAL(dp)                  :: help1, help2, help3, ampl, tempnorm, fac2, ntemp12

    ! Allocate arrays
    allocate(pdd12(nx,ny,12))
    allocate(rfr12(nx,ny,12))

    call init_pdd_lut()

    ! Calculate rain fraction and number of PDDs
    ! Huybrechts and De Wolde 1999 (C10), (C15)
    help1=sigma*360./12.
    fac2=pi/6.
    DO j=1,ny
      DO i=1,nx
        pdd(i,j)=0.0
        rfr(i,j)=0.0
        ! Seosonal amplitude
        ampl=abs(tmj(i,j)-tma(i,j))
        DO k=1,12
          ! Current temperature
          ntemp12=tma(i,j)+ampl*cos(fac2*(k-1))
          ntemp12=ntemp12/sigma
          help2=nintx*ntemp12/amax1(abs(ntemp12),valmax)
          ! Monthly PDDs
          pdd12(i,j,k)=help1*amax1(tabepdd(nint(help2)),ntemp12)
          ! Annual PDDS
          pdd(i,j)=pdd(i,j)+pdd12(i,j,k)  
          tempnorm=ntemp12-rainlimit/sigma
          help3=nintx*tempnorm/amax1(abs(tempnorm),valmax)  
          ! Monthly rain fraction
          rfr12(i,j,k)=taberf(nint(help3))+0.5
          ! Annual rain fraction
          rfr(i,j)=rfr(i,j)+rfr12(i,j,k)/12.  
        END DO
      END DO
    END DO


    deallocate(pdd12)
    deallocate(rfr12)


  END SUBROUTINE calculate_pdd_yearly



  SUBROUTINE calculate_pdd_monthly(nx, ny, sigma, rainlimit, tm, pdd, rfr)
    ! Positive degree day model, added by Heiko Goelzer, Feb 2022
    ! PDD model from Huybrechts and De Wolde 1999
    ! Monthly input fields

    IMPLICIT NONE



    ! Input variables: 
    INTEGER, INTENT(IN)  :: nx, ny ! grid size

    REAL(dp), INTENT(IN)  :: tm(nx,ny,12) ! monthly temperature 

    ! Output variables: 
    REAL(dp), INTENT(OUT)  :: pdd(nx,ny)
    REAL(dp), INTENT(OUT)  :: rfr(nx,ny)

    ! Local variables:
    REAL(dp), INTENT(IN)      :: sigma
    REAL(dp), INTENT(IN)      :: rainlimit

    INTEGER                   :: i, j, k
    REAL(dp), allocatable     :: pdd12(:,:,:)
    REAL(dp), allocatable     :: rfr12(:,:,:)
    REAL(dp)                  :: help1, help2, help3, ampl, tempnorm, ntemp12

    ! Allocate arrays
    allocate(pdd12(nx,ny,12))
    allocate(rfr12(nx,ny,12))

    call init_pdd_lut()

    ! Calculate rain fraction and number of PDDs
    ! Huybrechts and De Wolde 1999 (C10), (C15)
    help1=sigma*360./12.
    DO j=1,ny
      DO i=1,nx
        pdd(i,j)=0.0
        rfr(i,j)=0.0
        ! Seosonal cycle
        DO k=1,12
          ! Current temperature
          ntemp12=tm(i,j,k)
          ntemp12=ntemp12/sigma
          help2=nintx*ntemp12/amax1(abs(ntemp12),valmax)
          ! Monthly PDDs
          pdd12(i,j,k)=help1*amax1(tabepdd(nint(help2)),ntemp12)
          ! Annual PDDS
          pdd(i,j)=pdd(i,j)+pdd12(i,j,k)  
          tempnorm=ntemp12-rainlimit/sigma
          help3=nintx*tempnorm/amax1(abs(tempnorm),valmax)  
          ! Monthly rain fraction
          rfr12(i,j,k)=taberf(nint(help3))+0.5
          ! Annual rain fraction
          rfr(i,j)=rfr(i,j)+rfr12(i,j,k)/12.  
        END DO
      END DO
    END DO


    deallocate(pdd12)
    deallocate(rfr12)


  END SUBROUTINE calculate_pdd_monthly

  SUBROUTINE calculate_pdd_monthly_inout(nx, ny, tm12, pdd12, rfr12)
    ! Positive degree day model, added by Heiko Goelzer, Mar 2026
    ! PDD model from Huybrechts and De Wolde 1999
    ! Monthly input and output fields

    IMPLICIT NONE



    ! Input variables: 
    INTEGER, INTENT(IN)  :: nx, ny ! grid size

    REAL(dp), INTENT(IN)  :: tm12(nx,ny,12) ! monthly temperature (C)

    ! Output variables: 
    REAL(dp), INTENT(OUT)  :: pdd12(nx,ny,12)
    REAL(dp), INTENT(OUT)  :: rfr12(nx,ny,12)

    ! Local variables:
    INTEGER                   :: i, j, k
    REAL(dp), PARAMETER       :: sigma = 4.5
    REAL(dp), PARAMETER       :: rainlimit = 1.0
    REAL(dp)                  :: help1, help2, help3, ampl, tempnorm, ntemp

    call init_pdd_lut()

    ! Calculate rain fraction and number of PDDs
    ! Huybrechts and De Wolde 1999 (C10), (C15)
    help1=sigma*360./12.
    DO j=1,ny
      DO i=1,nx
        DO k=1,12
          ! Current temperature
          ntemp=tm12(i,j,k)
          ntemp=ntemp/sigma
          help2=nintx*ntemp/amax1(abs(ntemp),valmax)
          ! Monthly PDDs
          pdd12(i,j,k)=help1*amax1(tabepdd(nint(help2)),ntemp)
          tempnorm=ntemp-rainlimit/sigma
          help3=nintx*tempnorm/amax1(abs(tempnorm),valmax)  
          ! Monthly rain fraction
          rfr12(i,j,k)=taberf(nint(help3))+0.5
        END DO
      END DO
    END DO


  END SUBROUTINE calculate_pdd_monthly_inout


  SUBROUTINE calculate_pdd_monthly_inout_taj(nx, ny, sigma, rainlimit, tm12, pdd12, rfr12)
    ! Positive degree day model, added by Heiko Goelzer, Mar 2026
    ! PDD model from Huybrechts and De Wolde 1999
    ! Monthly input and output fields
    !   Seasonal cycle is parametised using annual mean and july temperature
    !   This is only for testing backwards compatability ! 

    IMPLICIT NONE



    ! Input variables: 
    INTEGER, INTENT(IN)  :: nx, ny ! grid size

    REAL(dp), INTENT(IN)  :: tm12(nx,ny,12) ! monthly temperature (C)

    ! Output variables: 
    REAL(dp), INTENT(OUT)  :: pdd12(nx,ny,12)
    REAL(dp), INTENT(OUT)  :: rfr12(nx,ny,12)

    ! Local variables:
    REAL(dp), INTENT(IN)      :: sigma
    REAL(dp), INTENT(IN)      :: rainlimit

    INTEGER                   :: i, j, k
    REAL(dp)                  :: help1, help2, help3, ampl, tempnorm, ntemp, fac2, tma

    call init_pdd_lut()

    ! Calculate rain fraction and number of PDDs
    ! Huybrechts and De Wolde 1999 (C10), (C15)
    help1=sigma*360./12.
    fac2=pi/6. 
    DO j=1,ny
      DO i=1,nx
        ! Annual mean
        tma = 0.0
        DO k=1,12
          tma = tma + tm12(i,j,k)/12.
        END DO
        ! Seosonal amplitude
        ampl=abs(tm12(i,j,7)-tma)
        DO k=1,12
          ! Current temperature, parameterised seasonal cyclce
          ! was monthly temp: ntemp=tm12(i,j,k)
          ! note -cos to aling seasonal cycle for monthly output 
          ntemp=tma-ampl*cos(fac2*(k-1)) 
          ntemp=ntemp/sigma
          help2=nintx*ntemp/amax1(abs(ntemp),valmax)
          ! Monthly PDDs
          pdd12(i,j,k)=help1*amax1(tabepdd(nint(help2)),ntemp)
          tempnorm=ntemp-rainlimit/sigma
          help3=nintx*tempnorm/amax1(abs(tempnorm),valmax)  
          ! Monthly rain fraction
          rfr12(i,j,k)=taberf(nint(help3))+0.5
        END DO
      END DO
    END DO


  END SUBROUTINE calculate_pdd_monthly_inout_taj



  SUBROUTINE melt_cascade_2d(nx, ny, ddfactorsnow, ddfactorice, pmax, snow, pdd, tpa, sir, abl)
    INTEGER,  INTENT(IN)  :: nx, ny
    REAL(dp), INTENT(IN)  :: ddfactorsnow, ddfactorice, pmax
    REAL(dp), INTENT(IN)  :: snow(nx,ny), pdd(nx,ny), tpa(nx,ny)
    REAL(dp), INTENT(OUT) :: sir(nx,ny), abl(nx,ny)
    INTEGER  :: i, j
    REAL(dp) :: pdds, ablv, sifm
    DO j = 1, ny
      DO i = 1, nx
        pdds = snow(i,j) / ddfactorsnow
        sifm = pmax * snow(i,j)
        IF (sifm > tpa(i,j)) sifm = tpa(i,j)
        IF (pdds <= pdd(i,j)) THEN
          ablv = (pdd(i,j) - pdds) * ddfactorice + snow(i,j)
        ELSE
          ablv = pdd(i,j) * ddfactorsnow
        END IF
        IF (ablv > tpa(i,j) + sifm) THEN
          sir(i,j) = 0.
        ELSE IF (ablv > tpa(i,j)) THEN
          sir(i,j) = tpa(i,j) + sifm - ablv
        ELSE IF (ablv > sifm) THEN
          sir(i,j) = sifm
        ELSE
          sir(i,j) = ablv
        END IF
        abl(i,j) = ablv - sifm
        IF (abl(i,j) < 0.) abl(i,j) = 0.
      END DO
    END DO
  END SUBROUTINE melt_cascade_2d


  SUBROUTINE melt_cascade_3d(nx, ny, ddfactorsnow, ddfactorice, pmax, snow, pdd, tp, sir, abl)
    INTEGER,  INTENT(IN)  :: nx, ny
    REAL(dp), INTENT(IN)  :: ddfactorsnow, ddfactorice, pmax
    REAL(dp), INTENT(IN)  :: snow(nx,ny,12), pdd(nx,ny,12), tp(nx,ny,12)
    REAL(dp), INTENT(OUT) :: sir(nx,ny,12), abl(nx,ny,12)
    INTEGER  :: i, j, k
    REAL(dp) :: pdds, ablv, sifm
    DO j = 1, ny
      DO i = 1, nx
        DO k = 1, 12
          pdds = snow(i,j,k) / ddfactorsnow
          sifm = pmax * snow(i,j,k)
          IF (sifm > tp(i,j,k)) sifm = tp(i,j,k)
          IF (pdds <= pdd(i,j,k)) THEN
            ablv = (pdd(i,j,k) - pdds) * ddfactorice + snow(i,j,k)
          ELSE
            ablv = pdd(i,j,k) * ddfactorsnow
          END IF
          IF (ablv > tp(i,j,k) + sifm) THEN
            sir(i,j,k) = 0.
          ELSE IF (ablv > tp(i,j,k)) THEN
            sir(i,j,k) = tp(i,j,k) + sifm - ablv
          ELSE IF (ablv > sifm) THEN
            sir(i,j,k) = sifm
          ELSE
            sir(i,j,k) = ablv
          END IF
          abl(i,j,k) = ablv - sifm
          IF (abl(i,j,k) < 0.) abl(i,j,k) = 0.
        END DO
      END DO
    END DO
  END SUBROUTINE melt_cascade_3d

END MODULE massbalance_module
