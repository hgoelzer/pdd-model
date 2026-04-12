! File name: massbalance_module.f90

MODULE massbalance_module

CONTAINS

  SUBROUTINE pdd_model_greenland_total_monthly_inout(nx, ny, tp, t2m, smb, snow, rain, sir, abl, pdd, rfr)
    ! Positive degree day model, uddated by Heiko Goelzer, Mar 2026
    ! Implements Greenland pdd model of Huybrechts and De Wolde 1999
    ! Forcing with total monthly fields, output monthly data

    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(1.0D0)  ! Kind of double precision numbers.

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

    ! Local variables
    REAL(dp), allocatable               :: tm(:,:,:)

    REAL(dp), PARAMETER                 :: ddfactorsnow = 0.003
    REAL(dp), PARAMETER                 :: ddfactorice = 0.008
    REAL(dp), PARAMETER                 :: pmax = 0.3 ! See update in Janssens and Huybrechts 2000
 
    INTEGER                             :: i, j, k
    REAL(dp)                            :: pdds, ablv, sifm

    ! Allocate arrays
    allocate(tm(nx,ny,12))

    ! Monthly temperature (C)
    tm = t2m - 273.15

    ! Determine number of positive degree days per year and rain fraction
    !call calculate_pdd_monthly_inout(nx, ny, tm, pdd, rfr)
    ! With parameterised seasonal cycle
    call calculate_pdd_monthly_inout_taj(nx, ny, tm, pdd, rfr)

    ! Distinguish rain and snow according to rain fraction
    rain = tp * rfr
    snow = tp - rain

    ! Melt calculation
    DO j=1,ny
       DO i=1,nx
          DO k=1,12

             ! pdd needed for snow melting
             pdds = snow(i,j,k)/ddfactorsnow
             ! potential for refreezing 
             sifm = pmax*snow(i,j,k)
             ! limit potential by total precipitation (tpa)
             IF(sifm.GT.tp(i,j,k)) sifm=tp(i,j,k)
             
             ! Estimate available melt 
             IF(pdds.LE.pdd(i,j,k)) THEN
                ! Remainig energy (pdd) used for ice melt
                ablv = (pdd(i,j,k)-pdds)*ddfactorice+snow(i,j,k)
             ELSE
                ! All energy (pdd) used for snow melt
                ablv = pdd(i,j,k)*ddfactorsnow
             ENDIF
             
             ! Calculate refreezing
             IF(ablv.GT.tp(i,j,k)+sifm) THEN
                ! entire snowpack melted, no refreezing
                sir(i,j,k) = 0.
             ELSEIF(ablv.GT.tp(i,j,k)) THEN
                sir(i,j,k) = tp(i,j,k)+sifm-ablv
             ELSEIF(ablv.GT.sifm) THEN
                sir(i,j,k) = sifm
             ELSE
                sir(i,j,k) = ablv
             ENDIF
             abl(i,j,k) = ablv - sifm
             ! Sanity check
             IF(abl(i,j,k).lt.0) abl(i,j,k)=0
             
          END DO
       END DO
    END DO

    ! no melt where insuffient energy
    WHERE (pdd.LE.0) 
     sir = 0
     abl = 0
    END WHERE

    smb = tp - abl - rain ! = snow - abl

  END SUBROUTINE pdd_model_greenland_total_monthly_inout



  SUBROUTINE pdd_model_greenland_total_monthly(nx, ny, tp, t2m, smb, snow, rain, sir, abl, pdd, rfr)
    ! Positive degree day model, uddated by Heiko Goelzer, Feb 2022
    ! Implements Greenland pdd model of Huybrechts and De Wolde 1999
    ! Forcing with total monthly fields

    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(1.0D0)  ! Kind of double precision numbers.

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

    ! Local variables
    REAL(dp), allocatable               :: tm(:,:,:)
    REAL(dp), allocatable               :: tpa(:,:)

    REAL(dp), PARAMETER                 :: ddfactorsnow = 0.003
    REAL(dp), PARAMETER                 :: ddfactorice = 0.008
    REAL(dp), PARAMETER                 :: pmax = 0.3 ! See update in Janssens and Huybrechts 2000
 
    INTEGER                             :: i, j
    REAL(dp)                            :: pdds, ablv, sifm

    ! Allocate arrays
    allocate(tm(nx,ny,12))
    allocate(tpa(nx,ny))

    ! Monthly temperature
    tm = t2m - 273.15

    ! Determine number of positive degree days per year and rain fraction
    call calculate_pdd_monthly(nx, ny, tm, pdd, rfr)

    ! total annual precip
    tpa(:,:) = SUM(tp(:,:,:),dim=3)

    ! Distinguish rain and snow according to rain fraction
    rain = tpa * rfr
    snow = tpa - rain

    ! Melt calculation
    DO j=1,ny
      DO i=1,nx

        ! pdd needed for snow melting
        pdds = snow(i,j)/ddfactorsnow
        ! potential for refreezing 
        sifm = pmax*snow(i,j)
        ! limit potential by total precipitation (tpa)
        IF(sifm.GT.tpa(i,j)) sifm=tpa(i,j)

        ! Estimate available melt 
        IF(pdds.LE.pdd(i,j)) THEN
         ! Remainig energy (pdd) used for ice melt
         ablv = (pdd(i,j)-pdds)*ddfactorice+snow(i,j)
        ELSE
         ! All energy (pdd) used for snow melt
         ablv = pdd(i,j)*ddfactorsnow
        ENDIF

        ! Calculate refreezing
        IF(ablv.GT.tpa(i,j)+sifm) THEN
         ! entire snowpack melted, no refreezing
         sir(i,j) = 0.
        ELSEIF(ablv.GT.tpa(i,j)) THEN
         sir(i,j) = tpa(i,j)+sifm-ablv
        ELSEIF(ablv.GT.sifm) THEN
         sir(i,j) = sifm
        ELSE
         sir(i,j) = ablv
        ENDIF
        abl(i,j) = ablv - sifm
        ! Sanity check
        IF(abl(i,j).lt.0) abl(i,j)=0

      END DO
    END DO

    ! no melt where insuffient energy
    WHERE (pdd.LE.0) 
     sir = 0
     abl = 0
    END WHERE

    ! smb = snow - abl
    smb = tpa - abl - rain

  END SUBROUTINE pdd_model_greenland_total_monthly



  SUBROUTINE pdd_model_greenland_total_yearly(nx, ny, acc, t2m, t2j, smb, snow, rain, sir, abl, pdd, rfr)
    ! Positive degree day model, uddated by Heiko Goelzer, Feb 2022
    ! Implements Greenland pdd model of Huybrechts and De Wolde 1999
    ! Forcing with total fields, not anomalies

    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(1.0D0)  ! Kind of double precision numbers.

    REAL, PARAMETER ::  pi = 2.0_dp * ACOS(0.0_dp)

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

    ! Local variables
    REAL(dp), allocatable               :: tma(:,:)
    REAL(dp), allocatable               :: tmj(:,:)

    REAL(dp), PARAMETER                 :: ddfactorsnow = 0.003
    REAL(dp), PARAMETER                 :: ddfactorice = 0.008
    REAL(dp), PARAMETER                 :: pmax = 0.3 ! See update in Janssens and Huybrechts 2000
 
    INTEGER                             :: i, j
    REAL(dp)                            :: pdds, ablv, sifm

    ! Allocate arrays
    allocate(tma(nx,ny))
    allocate(tmj(nx,ny))

    ! Global temperature perturbation
    tma = t2m - 273.15

    ! Summer temperature
    tmj = t2j - 273.15

    ! Determine number of positive degree days per year and rain fraction
    call calculate_pdd_yearly(nx, ny, tma, tmj, pdd, rfr)

    ! Distinguish rain and snow according to rain fraction
    rain = acc * rfr
    snow = acc - rain

    ! Melt calculation
    DO j=1,ny
      DO i=1,nx

        ! pdd needed for snow melting
        pdds = snow(i,j)/ddfactorsnow
        ! potential for refreezing 
        sifm = pmax*snow(i,j)
        ! limit potential by total precipitation (acc)
        IF(sifm.GT.acc(i,j)) sifm=acc(i,j)

        ! Estimate available melt 
        IF(pdds.LE.pdd(i,j)) THEN
         ! Remainig energy (pdd) used for ice melt
         ablv = (pdd(i,j)-pdds)*ddfactorice+snow(i,j)
        ELSE
         ! All energy (pdd) used for snow melt
         ablv = pdd(i,j)*ddfactorsnow
        ENDIF

        ! Calculate refreezing
        IF(ablv.GT.acc(i,j)+sifm) THEN
         ! entire snowpack melted, no refreezing
         sir(i,j) = 0.
        ELSEIF(ablv.GT.acc(i,j)) THEN
         sir(i,j) = acc(i,j)+sifm-ablv
        ELSEIF(ablv.GT.sifm) THEN
         sir(i,j) = sifm
        ELSE
         sir(i,j) = ablv
        ENDIF
        abl(i,j) = ablv - sifm
        ! Sanity check
        IF(abl(i,j).lt.0) abl(i,j)=0

      END DO
    END DO

    ! no melt where insuffient energy
    WHERE (pdd.LE.0) 
     sir = 0
     abl = 0
    END WHERE

    ! smb = snow - abl
    smb = acc - abl - rain

    deallocate(tma)
    deallocate(tmj)

  END SUBROUTINE pdd_model_greenland_total_yearly



  SUBROUTINE massbalance_pdd_model_greenland(nx, ny, lat, Hs, acc_PD, T_anomaly, smb)
    ! Positive degree day model, added by Heiko Goelzer, Jan 2017
    ! Implements Greenland pdd model of Huybrechts and De Wolde 1999

    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(1.0D0)  ! Kind of double precision numbers.

    REAL, PARAMETER ::  pi = 2.0_dp * ACOS(0.0_dp)

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

    REAL(dp), PARAMETER                 :: ddfactorsnow = 0.003
    REAL(dp), PARAMETER                 :: ddfactorice = 0.008
    REAL(dp), PARAMETER                 :: pmax = 0.3 ! See update in Janssens and Huybrechts 2000
 
    INTEGER                             :: i, j
    REAL(dp)                            :: pdds, ablv, sifm

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
    call calculate_pdd_yearly(nx, ny, tma, tmj, pdd, rfr)

    ! Distinguish rain and snow according to rain fraction
    rain = acc * rfr
    snow = acc - rain

    ! Melt calculation
    DO j=1,ny
      DO i=1,nx

        ! pdd needed for snow melting
        pdds = snow(i,j)/ddfactorsnow
        ! potential for refreezing 
        sifm = pmax*snow(i,j)
        ! limit potential by total precipitation (acc)
        IF(sifm.GT.acc(i,j)) sifm=acc(i,j)

        ! Estimate available melt 
        IF(pdds.LE.pdd(i,j)) THEN
         ! Remainig energy (pdd) used for ice melt
         ablv = (pdd(i,j)-pdds)*ddfactorice+snow(i,j)
        ELSE
         ! All energy (pdd) used for snow melt
         ablv = pdd(i,j)*ddfactorsnow
        ENDIF

        ! Calculate refreezing
        IF(ablv.GT.acc(i,j)+sifm) THEN
         ! entire snowpack melted, no refreezing
         sir(i,j) = 0.
        ELSEIF(ablv.GT.acc(i,j)) THEN
         sir(i,j) = acc(i,j)+sifm-ablv
        ELSEIF(ablv.GT.sifm) THEN
         sir(i,j) = sifm
        ELSE
         sir(i,j) = ablv
        ENDIF
        abl(i,j) = ablv - sifm
        ! Sanity check
        IF(abl(i,j).lt.0) abl(i,j)=0

      END DO
    END DO

    ! no melt where insuffient energy
    WHERE (pdd.LE.0) 
     sir = 0
     abl = 0
    END WHERE

    ! smb = snow - abl
    smb = acc - abl - rain

  END SUBROUTINE massbalance_pdd_model_greenland



  SUBROUTINE calculate_pdd_yearly(nx, ny, tma, tmj, pdd, rfr)
    ! Positive degree day model, added by Heiko Goelzer, Jan 2017
    ! PDD model from Huybrechts and De Wolde 1999
    ! Seasonal cycle is parametised using annual mean and july temperature

    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(1.0D0) ! Kind of double precision numbers.

    REAL, PARAMETER ::  pi = 2.0_dp * ACOS(0.0_dp)

    ! Input variables: 
    INTEGER, INTENT(IN)  :: nx, ny ! grid size

    REAL(dp), INTENT(IN)  :: tma(nx,ny)
    REAL(dp), INTENT(IN)  :: tmj(nx,ny)

    ! Output variables: 
    REAL(dp), INTENT(OUT)  :: pdd(nx,ny)
    REAL(dp), INTENT(OUT)  :: rfr(nx,ny)

    ! Local variables:
    LOGICAL, SAVE             :: first_call = .TRUE.
    INTEGER                   :: i, j, k
    REAL(dp), PARAMETER       :: sigma = 4.5 
    REAL(dp), PARAMETER       :: rainlimit = 1.0
    REAL(dp), PARAMETER       :: valmax = 6.0
    INTEGER,  PARAMETER       :: nintx=1200
    REAL(dp), allocatable     :: pdd12(:,:,:)
    REAL(dp), allocatable     :: rfr12(:,:,:)
    REAL(dp)                  :: help1, help2, help3, ampl, tempnorm, fac2, ntemp12

    ! PDD
    REAL(dp), SAVE            :: taberf(-nintx:nintx),tabepdd(-nintx:nintx)
    REAL(dp)                  :: deltax,sq2pi,fac1,fdx,help,xi,xj,yi,yj

    ! Allocate arrays
    allocate(pdd12(nx,ny,12))
    allocate(rfr12(nx,ny,12))

    
    ! ------------------------------------------------------------------------
    ! Calculate lookup tables for error function and expected PDD on first call
    ! Huybrechts and De Wolde 1999 (C10), (C15)

    IF(first_call) THEN
     taberf(0)=0.0
     deltax=valmax/nintx
     sq2pi=(2*pi)**(0.5)
     fac1=deltax/(2*sq2pi)
     tabepdd(0)=1./sq2pi
     xj=0.
     yj=1.
     DO i=1,nintx
       xi=xj
       yi=yj
       xj=xj+deltax
       yj=exp(-0.5*xj*xj)
       fdx=(yi+yj)*fac1
       taberf(i) =taberf(i-1)+fdx
       taberf(-i)=-taberf(i)
       help=yj/sq2pi+xj*taberf(i)
       tabepdd(i) =help+xj*0.5
       tabepdd(-i)=help-xj*0.5
     END DO
     first_call = .FALSE.
    END IF

    ! --------------------------------------------------------------
    
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



  SUBROUTINE calculate_pdd_monthly(nx, ny, tm, pdd, rfr)
    ! Positive degree day model, added by Heiko Goelzer, Feb 2022
    ! PDD model from Huybrechts and De Wolde 1999
    ! Monthly input fields

    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(1.0D0) ! Kind of double precision numbers.

    REAL, PARAMETER ::  pi = 2.0_dp * ACOS(0.0_dp)

    ! Input variables: 
    INTEGER, INTENT(IN)  :: nx, ny ! grid size

    REAL(dp), INTENT(IN)  :: tm(nx,ny,12) ! monthly temperature 

    ! Output variables: 
    REAL(dp), INTENT(OUT)  :: pdd(nx,ny)
    REAL(dp), INTENT(OUT)  :: rfr(nx,ny)

    ! Local variables:
    LOGICAL, SAVE             :: first_call = .TRUE.
    INTEGER                   :: i, j, k
    REAL(dp), PARAMETER       :: sigma = 4.5 
    REAL(dp), PARAMETER       :: rainlimit = 1.0
    REAL(dp), PARAMETER       :: valmax = 6.0
    INTEGER,  PARAMETER       :: nintx=1200
    REAL(dp), allocatable     :: pdd12(:,:,:)
    REAL(dp), allocatable     :: rfr12(:,:,:)
    REAL(dp)                  :: help1, help2, help3, ampl, tempnorm, ntemp12

    ! PDD
    REAL(dp), SAVE            :: taberf(-nintx:nintx),tabepdd(-nintx:nintx)
    REAL(dp)                  :: deltax,sq2pi,fac1,fdx,help,xi,xj,yi,yj

    ! Allocate arrays
    allocate(pdd12(nx,ny,12))
    allocate(rfr12(nx,ny,12))

    
    ! ------------------------------------------------------------------------
    ! Calculate lookup tables for error function and expected PDD on first call
    ! Huybrechts and De Wolde 1999 (C10), (C15)

    IF(first_call) THEN
     taberf(0)=0.0
     deltax=valmax/nintx
     sq2pi=(2*pi)**(0.5)
     fac1=deltax/(2*sq2pi)
     tabepdd(0)=1./sq2pi
     xj=0.
     yj=1.
     DO i=1,nintx
       xi=xj
       yi=yj
       xj=xj+deltax
       yj=exp(-0.5*xj*xj)
       fdx=(yi+yj)*fac1
       taberf(i) =taberf(i-1)+fdx
       taberf(-i)=-taberf(i)
       help=yj/sq2pi+xj*taberf(i)
       tabepdd(i) =help+xj*0.5
       tabepdd(-i)=help-xj*0.5
     END DO
     first_call = .FALSE.
    END IF

    ! --------------------------------------------------------------
    
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

    INTEGER, PARAMETER :: dp = KIND(1.0D0) ! Kind of double precision numbers.

    REAL, PARAMETER ::  pi = 2.0_dp * ACOS(0.0_dp)

    ! Input variables: 
    INTEGER, INTENT(IN)  :: nx, ny ! grid size

    REAL(dp), INTENT(IN)  :: tm12(nx,ny,12) ! monthly temperature (C)

    ! Output variables: 
    REAL(dp), INTENT(OUT)  :: pdd12(nx,ny,12)
    REAL(dp), INTENT(OUT)  :: rfr12(nx,ny,12)

    ! Local variables:
    LOGICAL, SAVE             :: first_call = .TRUE.
    INTEGER                   :: i, j, k
    REAL(dp), PARAMETER       :: sigma = 4.5 
    REAL(dp), PARAMETER       :: rainlimit = 1.0
    REAL(dp), PARAMETER       :: valmax = 6.0
    INTEGER,  PARAMETER       :: nintx=1200
    REAL(dp)                  :: help1, help2, help3, ampl, tempnorm, ntemp

    ! PDD
    REAL(dp), SAVE            :: taberf(-nintx:nintx),tabepdd(-nintx:nintx)
    REAL(dp)                  :: deltax,sq2pi,fac1,fdx,help,xi,xj,yi,yj

    ! ------------------------------------------------------------------------
    ! Calculate lookup tables for error function and expected PDD on first call
    ! Huybrechts and De Wolde 1999 (C10), (C15)

    IF(first_call) THEN
     taberf(0)=0.0
     deltax=valmax/nintx
     sq2pi=(2*pi)**(0.5)
     fac1=deltax/(2*sq2pi)
     tabepdd(0)=1./sq2pi
     xj=0.
     yj=1.
     DO i=1,nintx
       xi=xj
       yi=yj
       xj=xj+deltax
       yj=exp(-0.5*xj*xj)
       fdx=(yi+yj)*fac1
       taberf(i) =taberf(i-1)+fdx
       taberf(-i)=-taberf(i)
       help=yj/sq2pi+xj*taberf(i)
       tabepdd(i) =help+xj*0.5
       tabepdd(-i)=help-xj*0.5
     END DO
     first_call = .FALSE.
    END IF

    ! --------------------------------------------------------------
    
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


  SUBROUTINE calculate_pdd_monthly_inout_taj(nx, ny, tm12, pdd12, rfr12)
    ! Positive degree day model, added by Heiko Goelzer, Mar 2026
    ! PDD model from Huybrechts and De Wolde 1999
    ! Monthly input and output fields
    !   Seasonal cycle is parametised using annual mean and july temperature
    !   This is only for testing backwards compatability ! 

    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(1.0D0) ! Kind of double precision numbers.

    REAL, PARAMETER ::  pi = 2.0_dp * ACOS(0.0_dp)

    ! Input variables: 
    INTEGER, INTENT(IN)  :: nx, ny ! grid size

    REAL(dp), INTENT(IN)  :: tm12(nx,ny,12) ! monthly temperature (C)

    ! Output variables: 
    REAL(dp), INTENT(OUT)  :: pdd12(nx,ny,12)
    REAL(dp), INTENT(OUT)  :: rfr12(nx,ny,12)

    ! Local variables:
    LOGICAL, SAVE             :: first_call = .TRUE.
    INTEGER                   :: i, j, k
    REAL(dp), PARAMETER       :: sigma = 4.5 
    REAL(dp), PARAMETER       :: rainlimit = 1.0
    REAL(dp), PARAMETER       :: valmax = 6.0
    INTEGER,  PARAMETER       :: nintx=1200
    REAL(dp)                  :: help1, help2, help3, ampl, tempnorm, ntemp, fac2, tma

    ! PDD
    REAL(dp), SAVE            :: taberf(-nintx:nintx),tabepdd(-nintx:nintx)
    REAL(dp)                  :: deltax,sq2pi,fac1,fdx,help,xi,xj,yi,yj

    ! ------------------------------------------------------------------------
    ! Calculate lookup tables for error function and expected PDD on first call
    ! Huybrechts and De Wolde 1999 (C10), (C15)

    IF(first_call) THEN
     taberf(0)=0.0
     deltax=valmax/nintx
     sq2pi=(2*pi)**(0.5)
     fac1=deltax/(2*sq2pi)
     tabepdd(0)=1./sq2pi
     xj=0.
     yj=1.
     DO i=1,nintx
       xi=xj
       yi=yj
       xj=xj+deltax
       yj=exp(-0.5*xj*xj)
       fdx=(yi+yj)*fac1
       taberf(i) =taberf(i-1)+fdx
       taberf(-i)=-taberf(i)
       help=yj/sq2pi+xj*taberf(i)
       tabepdd(i) =help+xj*0.5
       tabepdd(-i)=help-xj*0.5
     END DO
     first_call = .FALSE.
    END IF

    ! --------------------------------------------------------------

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

  
END MODULE massbalance_module
