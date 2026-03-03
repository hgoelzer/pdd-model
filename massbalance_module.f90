! File name: massbalance_module.f90

MODULE massbalance_module

CONTAINS

  SUBROUTINE pdd_model_greenland_total(nx, ny, acc, t2m, t2j, smb)
    ! Positive degree day model, uddated by Heiko Goelzer, Feb 2022
    ! Implements Greenland pdd model of Huybrechts and De Wolde 1999
    ! Forcing with total fields, not anomalies

    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(1.0D0)  ! Kind of double precision numbers.
    REAL, PARAMETER :: triple_point_of_water = 273.16_dp

    REAL, PARAMETER ::  pi = 2.0_dp * ACOS(0.0_dp)

    ! --------------------------------------------------------------------------
    ! Declaration of global variables
    ! --------------------------------------------------------------------------

    ! Input variables: 
    INTEGER, INTENT(IN)  :: nx, ny ! grid size

    REAL(dp), INTENT(IN)  :: acc(ny,nx)    ! total yearly accumulation (m/yr)
    REAL(dp), INTENT(IN)  :: t2m(ny,nx)    ! annual mean 2m temperature (deg)
    REAL(dp), INTENT(IN)  :: t2j(ny,nx)    ! july 2m temperature (deg)

    ! Output variables: 
    REAL(dp), INTENT(OUT) :: smb(ny,nx)        ! surface mass balance (m/yr)

    ! Local variables
    REAL(dp), allocatable               :: snow(:,:)
    REAL(dp), allocatable               :: rain(:,:)
    REAL(dp), allocatable               :: sir(:,:)
    REAL(dp), allocatable               :: abl(:,:)           ! runoff (m/yr)
    REAL(dp), allocatable               :: pdd(:,:)
    REAL(dp), allocatable               :: rfr(:,:)
    REAL(dp), allocatable               :: s_prec(:,:)
    REAL(dp), allocatable               :: h_inv(:,:)

    REAL(dp), PARAMETER                 :: ddfactorsnow = 0.0033
    REAL(dp), PARAMETER                 :: ddfactorice = 0.0088
    REAL(dp), PARAMETER                 :: pmax = 0.3 ! See update in Janssens and Huybrechts 2000
 
    INTEGER                             :: i, j
    REAL(dp)                            :: pdds, ablv, sifm

    ! Allocate arrays
    allocate(snow(ny,nx))
    allocate(rain(ny,nx))
    allocate(sir(ny,nx))
    allocate(abl(ny,nx))
    allocate(pdd(ny,nx))
    allocate(rfr(ny,nx))
    allocate(s_prec(ny,nx))
    allocate(h_inv(ny,nx))

    ! Determine number of positive degree days per year and rain fraction
    call calculate_pdd_monthly(nx, ny, t2m, t2j, pdd, rfr)

    ! Distinguish rain and snow according to rain fraction
    rain = acc * rfr
    snow = acc - rain

    ! Melt calculation
    DO j=1,ny
      DO i=1,nx

        ! pdd needed for snow melting
        pdds = snow(j,i)/ddfactorsnow
        ! potential for refreezing 
        sifm = pmax*snow(j,i)
        ! limit potential by total precipitation (acc)
        IF(sifm.GT.acc(j,i)) sifm=acc(j,i)

        ! Estimate available melt 
        IF(pdds.LE.pdd(j,i)) THEN
         ! Remainig energy (pdd) used for ice melt
         ablv = (pdd(j,i)-pdds)*ddfactorice+snow(j,i)
        ELSE
         ! All energy (pdd) used for snow melt
         ablv = pdd(j,i)*ddfactorsnow
        ENDIF

        ! Calculate refreezing
        IF(ablv.GT.acc(j,i)+sifm) THEN
         ! entire snowpack melted, no refreezing
         sir(j,i) = 0.
        ELSEIF(ablv.GT.acc(j,i)) THEN
         sir(j,i) = acc(j,i)+sifm-ablv
        ELSEIF(ablv.GT.sifm) THEN
         sir(j,i) = sifm
        ELSE
         sir(j,i) = ablv
        ENDIF
        abl(j,i) = ablv - sifm
        ! Sanity check
        IF(abl(j,i).lt.0) abl(j,i)=0

      END DO
    END DO

    ! no melt where insuffient energy
    WHERE (pdd.LE.0) 
     sir = 0
     abl = 0
    END WHERE

    ! smb = snow - abl
    smb = acc - abl - rain

  END SUBROUTINE pdd_model_greenland_total



  SUBROUTINE massbalance_pdd_model_greenland(nx, ny, lat, Hs, acc_PD, T_anomaly, smb)
    ! Positive degree day model, added by Heiko Goelzer, Jan 2017
    ! Implements Greenland pdd model of Huybrechts and De Wolde 1999

    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(1.0D0)  ! Kind of double precision numbers.
    REAL, PARAMETER :: triple_point_of_water = 273.16_dp

    REAL, PARAMETER ::  pi = 2.0_dp * ACOS(0.0_dp)

    ! --------------------------------------------------------------------------
    ! Declaration of global variables
    ! --------------------------------------------------------------------------

    ! Input variables: 
    INTEGER, INTENT(IN)  :: nx, ny ! grid size

    REAL(dp), INTENT(IN)  :: lat(ny,nx)            ! latitude (deg+)
    REAL(dp), INTENT(IN)  :: Hs(ny,nx)             ! surface elevation (m)
    ! Hs_PD not in use; compare with Antarctic model
!    REAL(dp), DIMENSION(   ny,nx), INTENT(IN)  :: Hs_PD         ! reference surface elevation (m)

    REAL(dp), INTENT(IN)  :: acc_PD(ny,nx)    ! reference accumulation (m/yr)
    REAL(dp), INTENT(IN)  :: T_anomaly(ny,nx) ! temperature anomaly w.r.t. PD (deg)

    ! Output variables: 
    REAL(dp), INTENT(OUT) :: smb(ny,nx)        ! surface mass balance (m/yr)

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
    allocate(tma(ny,nx))
    allocate(snow(ny,nx))
    allocate(rain(ny,nx))
    allocate(sir(ny,nx))
    allocate(acc(ny,nx))
    allocate(abl(ny,nx))
    allocate(tmj(ny,nx))
    allocate(pdd(ny,nx))
    allocate(rfr(ny,nx))
    allocate(s_prec(ny,nx))
    allocate(h_inv(ny,nx))

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
    call calculate_pdd_monthly(nx, ny, tma, tmj, pdd, rfr)

    ! Distinguish rain and snow according to rain fraction
    rain = acc * rfr
    snow = acc - rain

    ! Melt calculation
    DO j=1,ny
      DO i=1,nx

        ! pdd needed for snow melting
        pdds = snow(j,i)/ddfactorsnow
        ! potential for refreezing 
        sifm = pmax*snow(j,i)
        ! limit potential by total precipitation (acc)
        IF(sifm.GT.acc(j,i)) sifm=acc(j,i)

        ! Estimate available melt 
        IF(pdds.LE.pdd(j,i)) THEN
         ! Remainig energy (pdd) used for ice melt
         ablv = (pdd(j,i)-pdds)*ddfactorice+snow(j,i)
        ELSE
         ! All energy (pdd) used for snow melt
         ablv = pdd(j,i)*ddfactorsnow
        ENDIF

        ! Calculate refreezing
        IF(ablv.GT.acc(j,i)+sifm) THEN
         ! entire snowpack melted, no refreezing
         sir(j,i) = 0.
        ELSEIF(ablv.GT.acc(j,i)) THEN
         sir(j,i) = acc(j,i)+sifm-ablv
        ELSEIF(ablv.GT.sifm) THEN
         sir(j,i) = sifm
        ELSE
         sir(j,i) = ablv
        ENDIF
        abl(j,i) = ablv - sifm
        ! Sanity check
        IF(abl(j,i).lt.0) abl(j,i)=0

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



  SUBROUTINE calculate_pdd_monthly(nx, ny, tma, tmj, pdd, rfr)
    ! Positive degree day model, added by Heiko Goelzer, Jan 2017
    ! PDD model from Huybrechts and De Wolde 1999

    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(1.0D0) ! Kind of double precision numbers.
    REAL, PARAMETER :: triple_point_of_water = 273.16_dp

    REAL, PARAMETER ::  pi = 2.0_dp * ACOS(0.0_dp)

    ! Input variables: 
    INTEGER, INTENT(IN)  :: nx, ny ! grid size

    REAL(dp), INTENT(IN)  :: tma(ny,nx)
    REAL(dp), INTENT(IN)  :: tmj(ny,nx)

    ! Output variables: 
    REAL(dp), INTENT(OUT)  :: pdd(ny,nx)
    REAL(dp), INTENT(OUT)  :: rfr(ny,nx)

    ! Local variables:
    LOGICAL, SAVE             :: first_call = .TRUE.
    INTEGER                   :: i, j, k
    REAL(dp), PARAMETER       :: sigma = 6.3 
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
    allocate(pdd12(12,ny,nx))
    allocate(rfr12(12,ny,nx))

    
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
        pdd(j,i)=0.0
        rfr(j,i)=0.0
        ! Seosonal amplitude
        ampl=abs(tmj(j,i)-tma(j,i))
        DO k=1,12
          ! Current temperature
          ntemp12=tma(j,i)+ampl*cos(fac2*(k-1))
          ntemp12=ntemp12/sigma
          help2=nintx*ntemp12/amax1(abs(ntemp12),valmax)
          ! Monthly PDDs
          pdd12(k,j,i)=help1*amax1(tabepdd(nint(help2)),ntemp12)
          ! Annual PDDS
          pdd(j,i)=pdd(j,i)+pdd12(k,j,i)  
          tempnorm=ntemp12-rainlimit/sigma
          help3=nintx*tempnorm/amax1(abs(tempnorm),valmax)  
          ! Monthly rain fraction
          rfr12(k,j,i)=taberf(nint(help3))+0.5
          ! Annual rain fraction
          rfr(j,i)=rfr(j,i)+rfr12(k,j,i)/12.  
        END DO
      END DO
    END DO


    deallocate(pdd12)
    deallocate(rfr12)


  END SUBROUTINE calculate_pdd_monthly



END MODULE massbalance_module
