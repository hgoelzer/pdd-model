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
         tp, t2m, smba, snowa, raina, sira, abla, pdda, rfra)

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
         tp, t2m, smb, snow, rain, sir, abl, pdd, rfr)

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
    real(dp) :: smb, snow, rain, sir, abl, pdd, rfr
    real(dp), parameter :: PR = 0.01_dp
    real(dp), parameter :: PMAX = 0.3_dp
    write(*,*) "Test 5: refreezing cap  (sir <= pmax * snow)"
    ! Mildly cold: some snow, some melt, test sir cap
    call call_model(-2.0_dp, PR, smb, snow, rain, sir, abl, pdd, rfr)
    call check(sir <= PMAX * snow + ATOL, "sir <= pmax * snow", sir, PMAX*snow)
    call check(sir >= 0.0_dp - ATOL,     "sir >= 0",           sir, 0.0_dp)
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

end program test_massbalance
