program test_massbalance
  ! Unit tests for massbalance_module.f90.
  ! Compiles without NetCDF — no forcing files needed.
  ! Build:  make test_massbalance
  ! Run:    ./test_massbalance.x

  use massbalance_module
  implicit none

  ! dp is imported from massbalance_module
  real(dp), parameter :: ATOL = 1.0e-10_dp   ! absolute tolerance for identities
  real(dp), parameter :: RTOL = 1.0e-4_dp    ! relative tolerance for physical checks

  integer :: nfail
  nfail = 0

  write(*,*) "=== test_massbalance: PDD physics unit tests ==="
  write(*,*)

  call test_cold_dry()
  call test_cold_wet()
  call test_warm_wet()
  call test_conservation()
  call test_bounds()
  call test_rainlimit_boundary()
  call test_retained()
  call test_melt_split()
  call test_2d_3d_equivalence()
  call test_total_runoff()
  call test_seasonalcycle_equivalence()
  call test_pmax_zero()


  write(*,*)
  if (nfail == 0) then
     write(*,*) "ALL TESTS PASSED"
     stop 0
  else
     write(*,'(a,i0,a)') " FAILED: ", nfail, " test(s)"
     stop 1
  end if

contains

  ! ── helpers ──────────────────────────────────────────────────────────────────

  subroutine check(cond, name, actual, expected)
    logical,          intent(in) :: cond
    character(len=*), intent(in) :: name
    real(dp),         intent(in) :: actual, expected
    if (cond) then
       write(*,'(a,a)') "  PASS  ", name
    else
       write(*,'(a,a,a,g16.6,a,g16.6)') "  FAIL  ", name, &
            "  actual=", actual, "  expected≈", expected
       nfail = nfail + 1
    end if
  end subroutine check

  subroutine checkl(cond, name, msg)
    logical,          intent(in) :: cond
    character(len=*), intent(in) :: name, msg
    if (cond) then
       write(*,'(a,a)') "  PASS  ", name
    else
       write(*,'(a,a,a,a)') "  FAIL  ", name, "  ", msg
       nfail = nfail + 1
    end if
  end subroutine checkl

  subroutine call_model(t_c, pr_m, smb, snow, rain, sir, abl, pdd, rfr)
    ! Convenience wrapper: uniform t_c (°C) and pr_m (m/month) over 1×1 grid.
    real(dp), intent(in)  :: t_c, pr_m
    real(dp), intent(out) :: smb, snow, rain, sir, abl, pdd, rfr
    integer,  parameter :: NX=1, NY=1
    real(dp) :: t2m(NX,NY,12), tp(NX,NY,12)
    real(dp) :: smba(NX,NY,12), snowa(NX,NY,12), raina(NX,NY,12)
    real(dp) :: sira(NX,NY,12), abla(NX,NY,12), pdda(NX,NY,12), rfra(NX,NY,12)
    ! Physics defaults (match driver defaults)
    real(dp), parameter :: ddfactorsnow = 0.00297_dp
    real(dp), parameter :: ddfactorice  = 0.00791_dp
    real(dp), parameter :: sigma        = 4.5_dp
    real(dp), parameter :: rainlimit    = 1.0_dp

    t2m = t_c + 273.15_dp   ! K
    tp  = pr_m              ! m ice equivalent per month

    call pdd_model_greenland_total_monthly_inout( &
         NX, NY, ddfactorsnow, ddfactorice, sigma, rainlimit, &
         0.3_dp, 1, tp, t2m, smba, snowa, raina, sira, abla, pdda, rfra)

    ! Reduce to scalars (annual mean of the 1×1 grid)
    smb  = sum(smba)  / 12.0_dp
    snow = sum(snowa) / 12.0_dp
    rain = sum(raina) / 12.0_dp
    sir  = sum(sira)  / 12.0_dp
    abl  = sum(abla)  / 12.0_dp
    pdd  = sum(pdda)  / 12.0_dp
    rfr  = sum(rfra)  / 12.0_dp
  end subroutine call_model

  ! ── Test 1: Cold and dry ─────────────────────────────────────────────────────
  subroutine test_cold_dry()
    real(dp) :: smb, snow, rain, sir, abl, pdd, rfr
    write(*,*) "Test 1: cold dry  (T=-30°C, pr=0)"
    call call_model(-30.0_dp, 0.0_dp, smb, snow, rain, sir, abl, pdd, rfr)
    call check(pdd  >= 0.0_dp,          "pdd >= 0",        pdd,  0.0_dp)
    call check(abs(rain) < ATOL,         "rain = 0",        rain, 0.0_dp)
    call check(abs(snow) < ATOL,         "snow = 0",        snow, 0.0_dp)
    ! At T=-30°C pdd is ~1e-7, giving abl~3e-10 — not exactly 0 since WHERE only fires at pdd<=0
    call check(abs(abl)  < 1.0e-8_dp,   "abl  ≈ 0",        abl,  0.0_dp)
    call check(abs(smb)  < 1.0e-8_dp,   "smb  ≈ 0",        smb,  0.0_dp)
    write(*,*)
  end subroutine test_cold_dry

  ! ── Test 2: Cold and wet (accumulation) ─────────────────────────────────────
  subroutine test_cold_wet()
    real(dp) :: smb, snow, rain, sir, abl, pdd, rfr
    real(dp), parameter :: PR = 0.01_dp  ! m ice/month
    write(*,*) "Test 2: cold wet  (T=-30°C, pr=0.01 m/month)"
    call call_model(-30.0_dp, PR, smb, snow, rain, sir, abl, pdd, rfr)
    call check(pdd  >= 0.0_dp,               "pdd >= 0",       pdd,  0.0_dp)
    call check(pdd  < 0.001_dp,              "pdd ≈ 0 (cold)", pdd,  0.0_dp)
    call check(rfr  < 0.01_dp,               "rfr ≈ 0 (all snow)", rfr, 0.0_dp)
    call check(abs(rain) < ATOL,              "rain = 0",       rain, 0.0_dp)
    call check(abs(snow - PR) < ATOL,         "snow = pr",      snow, PR)
    call check(abs(abl)  < ATOL,              "abl  = 0",       abl,  0.0_dp)
    call check(abs(smb - PR) < 1.0e-6_dp,    "smb ≈ pr",       smb,  PR)
    write(*,*)
  end subroutine test_cold_wet

  ! ── Test 3: Warm and wet (melt-dominated) ───────────────────────────────────
  subroutine test_warm_wet()
    real(dp) :: smb, snow, rain, sir, abl, pdd, rfr
    real(dp), parameter :: PR = 0.01_dp
    write(*,*) "Test 3: warm wet  (T=+10°C, pr=0.01 m/month)"
    call call_model(10.0_dp, PR, smb, snow, rain, sir, abl, pdd, rfr)
    call check(pdd  > 0.0_dp,           "pdd > 0 (warm)",   pdd, 1.0_dp)
    ! rfr = taberf((T-rainlimit)/sigma)+0.5; at T=+10, sigma=4.5: rfr ≈ 0.977
    call check(rfr  > 0.9_dp,           "rfr > 0.9 (mostly rain)", rfr, 1.0_dp)
    call check(rain > 0.9_dp * PR,      "rain > 0.9*pr",    rain, PR)
    call check(snow < 0.1_dp * PR,      "snow < 0.1*pr",    snow, 0.0_dp)
    call check(abl  > 0.0_dp,           "abl > 0",          abl,  0.0_dp)
    call check(smb  < 0.0_dp,           "smb < 0 (net melt)", smb, 0.0_dp)
    write(*,*)
  end subroutine test_warm_wet

  ! ── Test 4: Conservation identities ─────────────────────────────────────────
  ! These must hold exactly (they are enforced by the code structure).
  subroutine test_conservation()
    integer,  parameter :: NX=3, NY=4
    real(dp), parameter :: ddfactorsnow = 0.00297_dp
    real(dp), parameter :: ddfactorice  = 0.00791_dp
    real(dp), parameter :: sigma = 4.5_dp
    real(dp), parameter :: rainlimit = 1.0_dp
    real(dp) :: t2m(NX,NY,12), tp(NX,NY,12)
    real(dp) :: smb(NX,NY,12), snow(NX,NY,12), rain(NX,NY,12)
    real(dp) :: sir(NX,NY,12), abl(NX,NY,12), pdd(NX,NY,12), rfr(NX,NY,12)
    real(dp) :: max_err
    integer  :: i, j, k
    write(*,*) "Test 4: conservation identities on 3×4 grid with varied T,pr"

    ! Mixed grid: cold to warm across j, varying pr across i
    do j = 1, NY
       do i = 1, NX
          do k = 1, 12
             t2m(i,j,k) = (-20.0_dp + 10.0_dp*(j-1)) + 273.15_dp  ! -20 to +10°C
             tp(i,j,k)  = 0.001_dp * i                              ! 0.001 to 0.003 m/month
          end do
       end do
    end do
    call pdd_model_greenland_total_monthly_inout( &
         NX, NY, ddfactorsnow, ddfactorice, sigma, rainlimit, &
         0.3_dp, 1, tp, t2m, smb, snow, rain, sir, abl, pdd, rfr)

    ! snow + rain = tp (exact by construction: rain=rfr*tp; snow=tp-rain)
    max_err = maxval(abs(snow + rain - tp))
    call check(max_err < ATOL, "snow + rain = tp (all cells)", max_err, 0.0_dp)

    ! smb = tp - abl - rain (exact by construction)
    max_err = maxval(abs(smb - (tp - abl - rain)))
    call check(max_err < ATOL, "smb = tp - abl - rain", max_err, 0.0_dp)

    ! smb = snow - abl (follows from above + snow=tp-rain)
    max_err = maxval(abs(smb - (snow - abl)))
    call check(max_err < ATOL, "smb = snow - abl", max_err, 0.0_dp)

    ! abl >= 0 everywhere
    call checkl(all(abl >= -ATOL), "abl >= 0 everywhere", "negative abl found")

    ! sir >= 0 everywhere
    call checkl(all(sir >= -ATOL), "sir >= 0 everywhere", "negative sir found")

    ! pdd >= 0 everywhere
    call checkl(all(pdd >= -ATOL), "pdd >= 0 everywhere", "negative pdd found")

    ! rfr in [0, 1]
    call checkl(all(rfr >= -ATOL .and. rfr <= 1.0_dp + ATOL), &
         "rfr in [0,1]", "rfr out of range")

    write(*,*)
  end subroutine test_conservation

  ! ── Test 5: Explicit bound checks ───────────────────────────────────────────
  subroutine test_bounds()
    integer,  parameter :: NX=2, NY=2
    real(dp), parameter :: ddfactorsnow = 0.00297_dp
    real(dp), parameter :: ddfactorice  = 0.00791_dp
    real(dp), parameter :: PMAX = 0.3_dp
    real(dp) :: snow(NX,NY,12), pdd(NX,NY,12), tp(NX,NY,12)
    real(dp) :: sir(NX,NY,12), abl(NX,NY,12)
    real(dp) :: pdds, ablv, sifm, err
    integer  :: i, j, k
    write(*,*) "Test 5: refreezing cap  (sir <= min(pmax * snow, ablv))"
    ! Mixed regimes: full melt, marginal melt, tiny melt, no melt
    snow = 0.004_dp
    tp   = 1.0_dp
    pdd  = 0.0_dp
    pdd(1,1,:) = 5.0_dp      ! melt-dominated: ablv > sifm
    pdd(2,1,:) = 0.5_dp      ! marginal melt
    pdd(1,2,:) = 0.0001_dp   ! tiny melt: ablv < sifm
                             ! (2,2) stays cold: ablv = 0
    call melt_cascade_3d(NX, NY, ddfactorsnow, ddfactorice, PMAX, &
         snow, pdd, tp, sir, abl)
    err = 0.0_dp
    do k = 1, 12
       do j = 1, NY
          do i = 1, NX
             pdds = snow(i,j,k) / ddfactorsnow
             sifm = PMAX * snow(i,j,k)
             if (pdds <= pdd(i,j,k)) then
                ablv = (pdd(i,j,k) - pdds) * ddfactorice + snow(i,j,k)
             else
                ablv = pdd(i,j,k) * ddfactorsnow
             end if
             err = max(err, abs(sir(i,j,k) - min(ablv, sifm)), &
                            abs(abl(i,j,k) - max(ablv - sifm, 0.0_dp)))
          end do
       end do
    end do
    call check(err < ATOL, "sir = min(ablv,sifm); abl = max(ablv-sifm,0)", err, 0.0_dp)
    write(*,*)
  end subroutine test_bounds

  ! ── Test 6: Rain fraction at rainlimit boundary ──────────────────────────────
  subroutine test_rainlimit_boundary()
    real(dp) :: smb, snow, rain, sir, abl, pdd, rfr
    real(dp), parameter :: RAINLIMIT = 1.0_dp   ! default value
    real(dp), parameter :: PR = 0.01_dp
    write(*,*) "Test 6: rain fraction at T = rainlimit (rfr ≈ 0.5)"
    call call_model(RAINLIMIT, PR, smb, snow, rain, sir, abl, pdd, rfr)
    ! The erf distribution centred at rainlimit → rfr ≈ 0.5
    call check(abs(rfr - 0.5_dp) < 0.1_dp, "rfr ≈ 0.5 at T=rainlimit", rfr, 0.5_dp)
    write(*,*)
  end subroutine test_rainlimit_boundary

  ! ── Test 7: retained-amount identity  sir = ablv - abl ──────────────────────
  subroutine test_retained()
    integer,  parameter :: NX=4, NY=3
    real(dp), parameter :: ddfactorsnow = 0.00297_dp
    real(dp), parameter :: ddfactorice  = 0.00791_dp
    real(dp), parameter :: sigma = 4.5_dp
    real(dp), parameter :: rainlimit = 1.0_dp
    real(dp), parameter :: PMAX = 0.3_dp
    real(dp) :: t2m(NX,NY,12), tp(NX,NY,12)
    real(dp) :: smb(NX,NY,12), snow(NX,NY,12), rain(NX,NY,12)
    real(dp) :: sir(NX,NY,12), abl(NX,NY,12), pdd(NX,NY,12), rfr(NX,NY,12)
    real(dp) :: pdds, ablv, sifm, e1, e2
    integer  :: i, j, k
    write(*,*) "Test 7: retained amount  sir = ablv - abl = min(ablv, pmax*snow)"

    ! Varied grid: cold to warm across j, dry to wet across i
    do j = 1, NY
       do i = 1, NX
          do k = 1, 12
             t2m(i,j,k) = (-15.0_dp + 12.0_dp*(j-1)) + 273.15_dp
             tp(i,j,k)  = 0.001_dp * (i + j)
          end do
       end do
    end do
    call pdd_model_greenland_total_monthly_inout( &
         NX, NY, ddfactorsnow, ddfactorice, sigma, rainlimit, &
         0.3_dp, 1, tp, t2m, smb, snow, rain, sir, abl, pdd, rfr)

    ! Recompute ablv and sifm independently from the model's snow/pdd
    e1 = 0.0_dp
    e2 = 0.0_dp
    do k = 1, 12
       do j = 1, NY
          do i = 1, NX
             pdds = snow(i,j,k) / ddfactorsnow
             sifm = PMAX * snow(i,j,k)
             if (pdds <= pdd(i,j,k)) then
                ablv = (pdd(i,j,k) - pdds) * ddfactorice + snow(i,j,k)
             else
                ablv = pdd(i,j,k) * ddfactorsnow
             end if
             e1 = max(e1, abs(sir(i,j,k) - (ablv - abl(i,j,k))))
             e2 = max(e2, abs(sir(i,j,k) - min(ablv, sifm)))
          end do
       end do
    end do
    call check(e1 < ATOL, "sir = ablv - abl (all cells)", e1, 0.0_dp)
    call check(e2 < ATOL, "sir = min(ablv, pmax*snow) (all cells)", e2, 0.0_dp)
    write(*,*)
  end subroutine test_retained

  ! ── Test 8: melt split  sir + abl = ablv  (pdd > 0) ─────────────────────────
  subroutine test_melt_split()
    integer,  parameter :: NX=3, NY=2
    real(dp), parameter :: ddfactorsnow = 0.00297_dp
    real(dp), parameter :: ddfactorice  = 0.00791_dp
    real(dp), parameter :: PMAX = 0.3_dp
    real(dp) :: snow(NX,NY,12), pdd(NX,NY,12), tp(NX,NY,12)
    real(dp) :: sir(NX,NY,12), abl(NX,NY,12)
    real(dp) :: pdds, ablv, err
    integer  :: i, j, k
    write(*,*) "Test 8: melt split  sir + abl = ablv  (pdd > 0)"

    ! pdd spans the pdds threshold so both ablv branches are exercised
    do k = 1, 12
       do j = 1, NY
          do i = 1, NX
             snow(i,j,k) = 0.002_dp * (i + 0.5_dp*j)
             pdd(i,j,k)  = 0.5_dp * (i*j) + 0.1_dp*k
             tp(i,j,k)   = 10.0_dp * snow(i,j,k)
          end do
       end do
    end do
    call melt_cascade_3d(NX, NY, ddfactorsnow, ddfactorice, PMAX, &
         snow, pdd, tp, sir, abl)
    err = 0.0_dp
    do k = 1, 12
       do j = 1, NY
          do i = 1, NX
             if (pdd(i,j,k) > 0.0_dp) then
                pdds = snow(i,j,k) / ddfactorsnow
                if (pdds <= pdd(i,j,k)) then
                   ablv = (pdd(i,j,k) - pdds) * ddfactorice + snow(i,j,k)
                else
                   ablv = pdd(i,j,k) * ddfactorsnow
                end if
                err = max(err, abs(sir(i,j,k) + abl(i,j,k) - ablv))
             end if
          end do
       end do
    end do
    call check(err < ATOL, "sir + abl = ablv where pdd > 0", err, 0.0_dp)
    write(*,*)
  end subroutine test_melt_split

  ! ── Test 9: 2D/3D cascade equivalence ───────────────────────────────────────
  subroutine test_2d_3d_equivalence()
    integer,  parameter :: NX=3, NY=2
    real(dp), parameter :: ddfactorsnow = 0.00297_dp
    real(dp), parameter :: ddfactorice  = 0.00791_dp
    real(dp), parameter :: PMAX = 0.3_dp
    real(dp) :: snow2d(NX,NY), pdd2d(NX,NY), tpa(NX,NY)
    real(dp) :: sir2d(NX,NY), abl2d(NX,NY)
    real(dp) :: snow3d(NX,NY,12), pdd3d(NX,NY,12), tp3d(NX,NY,12)
    real(dp) :: sir3d(NX,NY,12), abl3d(NX,NY,12)
    real(dp) :: err
    integer  :: k
    write(*,*) "Test 9: 2D/3D cascade equivalence (same cellwise equations)"

    do k = 1, NY
       snow2d(:,k) = 0.001_dp * (1 + k)
       pdd2d(:,k)  = 0.5_dp * (1 + k)   ! spans pdds threshold
       tpa(:,k)    = 12.0_dp * snow2d(:,k)
    end do
    do k = 1, 12
       snow3d(:,:,k) = snow2d
       pdd3d(:,:,k)  = pdd2d
       tp3d(:,:,k)   = tpa
    end do
    call melt_cascade_2d(NX, NY, ddfactorsnow, ddfactorice, PMAX, &
         snow2d, pdd2d, tpa, sir2d, abl2d)
    call melt_cascade_3d(NX, NY, ddfactorsnow, ddfactorice, PMAX, &
         snow3d, pdd3d, tp3d, sir3d, abl3d)
    err = 0.0_dp
    do k = 1, 12
       err = max(err, maxval(abs(sir3d(:,:,k) - sir2d)), &
                       maxval(abs(abl3d(:,:,k) - abl2d)))
    end do
    call check(err == 0.0_dp, "2D and 3D cascades identical", err, 0.0_dp)
    write(*,*)
  end subroutine test_2d_3d_equivalence

  ! ── Test 10: total-runoff output convention  mrro = abl + rain ──────────────
  ! The drivers write mrro as meltwater runoff + rainfall (rain drains
  ! directly). The SMB closure acabf = pr - mrro then holds exactly.
  subroutine test_total_runoff()
    integer,  parameter :: NX=3, NY=4
    real(dp), parameter :: ddfactorsnow = 0.00297_dp
    real(dp), parameter :: ddfactorice  = 0.00791_dp
    real(dp), parameter :: sigma = 4.5_dp
    real(dp), parameter :: rainlimit = 1.0_dp
    real(dp) :: t2m(NX,NY,12), tp(NX,NY,12)
    real(dp) :: smb(NX,NY,12), snow(NX,NY,12), rain(NX,NY,12)
    real(dp) :: sir(NX,NY,12), abl(NX,NY,12), pdd(NX,NY,12), rfr(NX,NY,12)
    real(dp) :: max_err
    integer  :: i, j, k
    write(*,*) "Test 10: total-runoff convention  mrro = abl + rain"

    ! Same mixed grid as Test 4: cold to warm across j, varying pr across i
    do j = 1, NY
       do i = 1, NX
          do k = 1, 12
             t2m(i,j,k) = (-20.0_dp + 10.0_dp*(j-1)) + 273.15_dp  ! -20 to +10°C
             tp(i,j,k)  = 0.001_dp * i                            ! 0.001 to 0.003 m/month
          end do
       end do
    end do
    call pdd_model_greenland_total_monthly_inout( &
         NX, NY, ddfactorsnow, ddfactorice, sigma, rainlimit, &
         0.3_dp, 1, tp, t2m, smb, snow, rain, sir, abl, pdd, rfr)

    ! smb = tp - (abl + rain): the driver-level identity acabf = pr - mrro
    max_err = maxval(abs(smb - (tp - (abl + rain))))
    call check(max_err < ATOL, "smb = tp - (abl + rain)  [acabf = pr - mrro]", &
         max_err, 0.0_dp)

    ! On the warm wet row both melt and rain contribute to total runoff
    if (maxval(pdd(:,NY,:)) > 0.0_dp .and. minval(rain(:,NY,:)) > 0.0_dp) then
       call checkl(all(abl(:,NY,:) > 0.0_dp), &
            "melt runoff > 0 on warm row (mrro > rain there)", "no melt runoff")
       call checkl(all(rain(:,NY,:) > 0.0_dp), &
            "rain > 0 on warm row (mrro > abl there)", "no rain")
    else
       call checkl(.false., "warm wet row present", &
            "test grid has no warm wet cells — extend grid")
    end if
    write(*,*)
  end subroutine test_total_runoff

  ! ── Test 11: seasonalcycle=0 (actual monthly temps) ≡ taj for a cosine ──────
  ! taj reconstructs the annual cycle from the annual mean and the July
  ! temperature. For a sinusoid peaking in July the reconstruction is exact,
  ! so both seasonal-cycle branches must agree to round-off.
  subroutine test_seasonalcycle_equivalence()
    integer,  parameter :: NX=2, NY=3
    real(dp), parameter :: ddfactorsnow = 0.00297_dp
    real(dp), parameter :: ddfactorice  = 0.00791_dp
    real(dp), parameter :: sigma = 4.5_dp
    real(dp), parameter :: rainlimit = 1.0_dp
    real(dp), parameter :: PMAX = 0.3_dp
    real(dp) :: t2m(NX,NY,12), tp(NX,NY,12)
    real(dp) :: smb1(NX,NY,12), snow1(NX,NY,12), rain1(NX,NY,12)
    real(dp) :: sir1(NX,NY,12), abl1(NX,NY,12), pdd1(NX,NY,12), rfr1(NX,NY,12)
    real(dp) :: smb0(NX,NY,12), snow0(NX,NY,12), rain0(NX,NY,12)
    real(dp) :: sir0(NX,NY,12), abl0(NX,NY,12), pdd0(NX,NY,12), rfr0(NX,NY,12)
    real(dp) :: err
    integer  :: i, j, k
    write(*,*) "Test 11: seasonalcycle=0 (actual temps) = taj for July-peak cosine"

    ! NH cosine climatology peaking in July: T = Tmean + ampl*cos(2pi(k-7)/12)
    do j = 1, NY
       do i = 1, NX
          do k = 1, 12
             t2m(i,j,k) = (-5.0_dp + 2.0_dp*(j-1) + 8.0_dp* &
                  cos(2.0_dp*3.141592653589793_dp*dble(k-7)/12.0_dp)) + 273.15_dp
             tp(i,j,k)  = 0.002_dp * (i + j)
          end do
       end do
    end do
    call pdd_model_greenland_total_monthly_inout( &
         NX, NY, ddfactorsnow, ddfactorice, sigma, rainlimit, &
         PMAX, 1, tp, t2m, smb1, snow1, rain1, sir1, abl1, pdd1, rfr1)
    call pdd_model_greenland_total_monthly_inout( &
         NX, NY, ddfactorsnow, ddfactorice, sigma, rainlimit, &
         PMAX, 0, tp, t2m, smb0, snow0, rain0, sir0, abl0, pdd0, rfr0)

    err = 0.0_dp
    do k = 1, 12
       err = max(err, maxval(abs(pdd1(:,:,k) - pdd0(:,:,k))), &
                       maxval(abs(rfr1(:,:,k) - rfr0(:,:,k))), &
                       maxval(abs(smb1(:,:,k) - smb0(:,:,k))))
    end do
    call check(err < 1.0e-8_dp, "seasonalcycle 0 vs 1 identical for July-peak cosine", &
         err, 0.0_dp)
    write(*,*)
  end subroutine test_seasonalcycle_equivalence

  ! ── Test 12: pmax=0 disables retention ──────────────────────────────────────
  ! With pmax=0 no meltwater is retained (sir=0 exactly) and melt runoff is
  ! at least as large as with pmax=0.3.
  subroutine test_pmax_zero()
    integer,  parameter :: NX=2, NY=3
    real(dp), parameter :: ddfactorsnow = 0.00297_dp
    real(dp), parameter :: ddfactorice  = 0.00791_dp
    real(dp), parameter :: sigma = 4.5_dp
    real(dp), parameter :: rainlimit = 1.0_dp
    real(dp) :: t2m(NX,NY,12), tp(NX,NY,12)
    real(dp) :: smb1(NX,NY,12), snow1(NX,NY,12), rain1(NX,NY,12)
    real(dp) :: sir1(NX,NY,12), abl1(NX,NY,12), pdd1(NX,NY,12), rfr1(NX,NY,12)
    real(dp) :: smb0(NX,NY,12), snow0(NX,NY,12), rain0(NX,NY,12)
    real(dp) :: sir0(NX,NY,12), abl0(NX,NY,12), pdd0(NX,NY,12), rfr0(NX,NY,12)
    integer  :: i, j, k
    write(*,*) "Test 12: pmax=0 disables retention (sir=0, runoff >= pmax=0.3 case)"

    do j = 1, NY
       do i = 1, NX
          do k = 1, 12
             t2m(i,j,k) = (-10.0_dp + 8.0_dp*(j-1)) + 273.15_dp  ! -10 to +6 C
             tp(i,j,k)  = 0.002_dp * i
          end do
       end do
    end do
    call pdd_model_greenland_total_monthly_inout( &
         NX, NY, ddfactorsnow, ddfactorice, sigma, rainlimit, &
         0.3_dp, 1, tp, t2m, smb1, snow1, rain1, sir1, abl1, pdd1, rfr1)
    call pdd_model_greenland_total_monthly_inout( &
         NX, NY, ddfactorsnow, ddfactorice, sigma, rainlimit, &
         0.0_dp, 1, tp, t2m, smb0, snow0, rain0, sir0, abl0, pdd0, rfr0)

    call check(maxval(abs(sir0)) == 0.0_dp, "pmax=0 gives sir=0 exactly", &
         maxval(abs(sir0)), 0.0_dp)
    call checkl(all(abl0 >= abl1 - 1.0e-12_dp), &
         "pmax=0 melt runoff >= pmax=0.3 melt runoff", "runoff decreased with pmax=0")
    write(*,*)
  end subroutine test_pmax_zero

end program test_massbalance
