!DEC$ FREEFORM
!=======================================================================
!  Abaqus/Standard User Subroutine: DLOAD + UEXTERNALDB
!  Model: External Helical Gear Pair - Wheel (Crown) Moving Contact Load
!         Crown Gear & Differential Carrier Assembly (Gear_Shaft_geometry.inp)
!  Author: Victor Maia (vhfm08@gmail.com)
!
!  References:
!    - ISO 21771: Cylindrical gears - Concepts and geometry
!    - ISO 6336-1: Calculation of load capacity of spur and helical gears
!
!  Conventions:
!    - Wheel axis: Global Z (centre at X = 0, Y = 0)
!    - Theta: Right-hand rule about +Z, from +X towards +Y
!    - Face width: Z in [0.0, 25.0] mm
!    - Units: mm, N, MPa, s (Torque: N.mm)
!    - Flank surfaces: SURF_TOOTH (drive flank), SURF_TOOTH_COAST (coast flank)
!    - Load cycle: one flank per analysis step (step 1 drive, step 2 coast /
!      regeneration at the same rotation direction); a single-step analysis
!      with SURF_TOOTH only reproduces the unidirectional drive case
!=======================================================================

module gear_contact_mod
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none

    ! ------------------------------------------------------------------
    ! 1. USER INPUT DATA ONLY (No derived values here)
    ! ------------------------------------------------------------------
    ! Gear pair nominal data
    real(real64), parameter :: z_wheel      = 60.0_real64    ! z2 [-]
    real(real64), parameter :: z_pinion     = 17.0_real64    ! z1 [-]
    real(real64), parameter :: mn           = 2.631213_real64 ! Normal module [mm] (mt = 180/62 = 2.9032 mm)
    real(real64), parameter :: alpha_n_deg  = 18.2565_real64 ! Normal pressure angle [deg] (alpha_t = 20 deg)
    real(real64), parameter :: beta_deg     = 25.0_real64    ! Helix angle at reference d [deg]
    real(real64), parameter :: hand_helix   = 1.0_real64     ! +1.0 right-hand, -1.0 left-hand

    ! Tip diameters & operating centre distance (set 0.0 for standard d + 2*mn)
    real(real64), parameter :: da_wheel_in  = 180.0_real64   ! Wheel tip diameter da2 [mm] (0.0 = auto)
    real(real64), parameter :: da_pinion_in = 0.0_real64     ! Pinion tip diameter da1 [mm] (0.0 = auto)
    real(real64), parameter :: a_centre_in  = 0.0_real64     ! Working centre distance a [mm] (0.0 = auto)

    ! Wheel axial bounds & pinion position
    real(real64), parameter :: z_face_1     = 0.0_real64     ! Face 1 [mm]
    real(real64), parameter :: z_face_2     = 25.0_real64    ! Face 2 [mm]
    real(real64), parameter :: theta_p_deg  = 0.0_real64     ! Pinion centre azimuth [deg]

    ! Operating condition
    integer(int32), parameter :: rot_dir    = 1              ! +1: CCW, -1: CW about +Z

    ! Load cycle per analysis step (kstep): +1 drive flank, -1 coast flank (regeneration)
    integer(int32), parameter :: n_cycle    = 2              ! Number of tabulated steps [-]
    integer(int32), parameter :: drive_mode_step(n_cycle) = [1, -1]
    real(real64), parameter :: torque_step_nm(n_cycle)    = [270.0_real64, 270.0_real64] ! Wheel torque [N.m]

    ! Contact profile band & time kinematics
    real(real64), parameter :: hb           = 1.0_real64     ! Contact band half-width [mm]
    real(real64), parameter :: w_end        = 0.5_real64     ! Entry/exit relief factor [-]
    real(real64), parameter :: z_taper      = 0.25_real64    ! Relief taper fraction [-]
    real(real64), parameter :: t_ramp       = 0.20_real64    ! Stationary load ramp duration [s]
    real(real64), parameter :: t_roll       = 1.00_real64    ! Rolling motion duration [s]
    real(real64), parameter :: n_pitches    = 60.0_real64    ! Full 360-deg revolution (60 pitches) [-]
    character(len=*), parameter :: surf_tag  = 'SURF_TOOTH'       ! Drive flank surface
    character(len=*), parameter :: coast_tag = 'SURF_TOOTH_COAST' ! Coast flank surface

    ! Mathematical constants
    real(real64), parameter :: pi     = 4.0_real64 * atan(1.0_real64)
    real(real64), parameter :: two_pi = 2.0_real64 * pi
    real(real64), parameter :: d2r    = pi / 180.0_real64

    ! ------------------------------------------------------------------
    ! 2. RUNTIME-DERIVED VARIABLES (Computed at t=0 by init_contact_tables)
    ! ------------------------------------------------------------------
    real(real64) :: mt, alpha_t, tan_betab, theta_p
    real(real64) :: r_wheel, r_pinion, rb_wheel, rb_pinion, ra_wheel, ra_pinion
    real(real64) :: a_dist, alpha_wt, u_a, u_e, g_alpha, p_bt, b_face, z_mid
    real(real64) :: z_lo, z_hi
    real(real64) :: rb2, inv_p_bt, inv_ga, inv_rb_hb, r_min2, r_max2
    real(real64) :: chx(2), p_unit, omega

    ! Flank index: 1 = drive flank, 2 = coast flank
    integer(int32), parameter :: n_tab = 128
    real(real64)              :: w_tab(0:n_tab, 2)
    logical                   :: is_init = .false.

contains

    !-------------------------------------------------------------------
    ! Pure function: signed parabolic profile coordinate q
    !-------------------------------------------------------------------
    pure function profile_q(u, dist) result(q)
        real(real64), intent(in) :: u, dist
        real(real64)             :: q
        q = dist * (u + 0.5_real64 * dist) * inv_rb_hb
    end function profile_q

    !-------------------------------------------------------------------
    ! Calculates all gear geometry, lines of action, and normalisation
    !-------------------------------------------------------------------
    subroutine init_contact_tables()
        integer(int32) :: i, j, k, m, ny, kmax
        real(real64)   :: beta, alpha_n, da_w, da_p, aw_calc, cos_awt
        real(real64)   :: dz, phi, z_pos, u_pos, w_sum, q_hi, q_lo, frac_cut

        if (is_init) return

        ! Angular conversions
        beta    = beta_deg * d2r
        alpha_n = alpha_n_deg * d2r
        theta_p = theta_p_deg * d2r

        ! Transverse module and pressure angle
        mt         = mn / cos(beta)
        alpha_t    = atan(tan(alpha_n) / cos(beta))
        tan_betab  = tan(beta) * cos(alpha_t)

        ! Reference & base radii
        r_wheel    = 0.5_real64 * mt * z_wheel
        r_pinion   = 0.5_real64 * mt * z_pinion
        rb_wheel   = r_wheel * cos(alpha_t)
        rb_pinion  = r_pinion * cos(alpha_t)

        ! Tip radii (use user inputs if provided, else standard addendum)
        da_w = da_wheel_in
        if (da_w <= 0.0_real64) da_w = 2.0_real64 * (r_wheel + mn)
        ra_wheel = 0.5_real64 * da_w

        da_p = da_pinion_in
        if (da_p <= 0.0_real64) da_p = 2.0_real64 * (r_pinion + mn)
        ra_pinion = 0.5_real64 * da_p

        ! Operating centre distance and working transverse pressure angle
        aw_calc = a_centre_in
        if (aw_calc <= 0.0_real64) aw_calc = r_wheel + r_pinion
        a_dist = aw_calc

        cos_awt = (rb_wheel + rb_pinion) / a_dist
        if (cos_awt >= 1.0_real64 .or. cos_awt <= 0.0_real64) then
            call stdb_abqerr(-3, "gear_contact: invalid centre distance (cos alpha_wt out of bounds)", 0, 0.0d0, " ")
            call xit
        end if
        alpha_wt = acos(cos_awt)

        ! Active path of contact on wheel line of action (u_A to u_E)
        u_e = sqrt(ra_wheel**2 - rb_wheel**2)
        u_a = a_dist * sin(alpha_wt) - sqrt(ra_pinion**2 - rb_pinion**2)
        if (u_a <= 0.0_real64) then
            call stdb_abqerr(-3, "gear_contact: pinion tip penetrates wheel base cylinder", 0, 0.0d0, " ")
            call xit
        end if

        g_alpha = u_e - u_a
        if (g_alpha <= 0.0_real64) then
            call stdb_abqerr(-3, "gear_contact: zero or negative path of contact", 0, 0.0d0, " ")
            call xit
        end if

        ! Base pitch and face geometry (cached limits for fast axial check)
        p_bt   = two_pi * rb_wheel / z_wheel
        b_face = abs(z_face_2 - z_face_1)
        z_mid  = 0.5_real64 * (z_face_1 + z_face_2)
        z_lo   = min(z_face_1, z_face_2)
        z_hi   = max(z_face_1, z_face_2)

        ! Pre-computed factors for fast integration point evaluation
        rb2        = rb_wheel**2
        inv_p_bt   = 1.0_real64 / p_bt
        inv_ga     = 1.0_real64 / g_alpha
        inv_rb_hb  = 1.0_real64 / (rb_wheel * hb)

        ! Mathematically exact radial envelope bounds (|q| <= 1.0) with 5% safety margin
        ! Captures 100% of the active band at tooth root while filtering non-active dedendum
        r_min2     = rb2 + max(u_a**2 - 2.05_real64 * rb_wheel * hb, 0.0_real64)
        r_max2     = ra_wheel**2 + 2.05_real64 * rb_wheel * hb

        ! Helical phase slope per flank (mirrored on the coast flank) and
        ! peak pressure per unit transmitted torque [MPa / (N.m)]
        chx(1)     = real(rot_dir, real64) * hand_helix * tan_betab
        chx(2)     = -chx(1)
        p_unit     = 1.0e3_real64 / (rb_wheel * (4.0_real64 / 3.0_real64) * hb)

        ! Rotational velocity: exactly n_pitches rolled over duration t_roll
        if (t_roll <= 0.0_real64) then
            call stdb_abqerr(-3, "gear_contact: t_roll must be positive", 0, 0.0d0, " ")
            call xit
        end if
        omega      = n_pitches * (two_pi / z_wheel) / t_roll

        ! Contact line length normalisation table with tip/root boundary truncation
        ! Guarantees constant transmitted torque (residual ripple < 0.3%)
        ny   = 200
        dz   = b_face / real(ny, real64)
        kmax = int((g_alpha + b_face * tan_betab) * inv_p_bt) + 2

        do m = 1, 2
            do i = 0, n_tab - 1
                phi   = p_bt * real(i, real64) / real(n_tab, real64)
                w_sum = 0.0_real64

                do j = 1, ny
                    z_pos = -0.5_real64 * b_face + (real(j, real64) - 0.5_real64) * dz
                    do k = -kmax, kmax
                        u_pos = u_a + phi - chx(m) * z_pos + real(k, real64) * p_bt

                        ! Boundary truncation fraction accounting for parabola cut at tip/root
                        q_hi = min(max((u_e**2 - u_pos**2) * 0.5_real64 * inv_rb_hb, -1.0_real64), 1.0_real64)
                        q_lo = min(max((r_min2 - rb2 - u_pos**2) * 0.5_real64 * inv_rb_hb, -1.0_real64), 1.0_real64)
                        frac_cut = 0.75_real64 * ((q_hi - q_hi**3 / 3.0_real64) - (q_lo - q_lo**3 / 3.0_real64))

                        w_sum = w_sum + flank_weight(u_pos) * frac_cut * dz
                    end do
                end do

                w_tab(i, m) = 1.0_real64 / max(w_sum, 1.0e-12_real64)
            end do

            w_tab(n_tab, m) = w_tab(0, m)
        end do
        is_init = .true.
    end subroutine init_contact_tables

    pure function flank_weight(u) result(w)
        real(real64), intent(in) :: u
        real(real64)             :: w, xi, s

        xi = (u - u_a) * inv_ga
        if (xi < 0.0_real64 .or. xi > 1.0_real64) then
            w = 0.0_real64
        else
            s = min(xi, 1.0_real64 - xi) / z_taper
            if (s >= 1.0_real64) then
                w = 1.0_real64
            else
                w = w_end + (1.0_real64 - w_end) * s * s * (3.0_real64 - 2.0_real64 * s)
            end if
        end if
    end function flank_weight

    pure function time_ramp(t) result(amp)
        real(real64), intent(in) :: t
        real(real64)             :: amp, tau

        if (t <= 0.0_real64) then
            amp = 0.0_real64
        else if (t < t_ramp) then
            tau = t / t_ramp
            amp = tau * tau * (3.0_real64 - 2.0_real64 * tau)
        else
            amp = 1.0_real64
        end if
    end function time_ramp

end module gear_contact_mod


!=======================================================================
!  UEXTERNALDB: One-time setup at step initialization
!=======================================================================
subroutine uexternaldb(lop, lrst, time, dtime, kstep, kinc)
    use gear_contact_mod, only: init_contact_tables
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none

    integer(int32), intent(in) :: lop, lrst, kstep, kinc
    real(real64),   intent(in) :: time(2), dtime

    if (lop == 0 .or. lop == 4) call init_contact_tables()
end subroutine uexternaldb


!=======================================================================
!  DLOAD: Surface distributed pressure evaluation
!=======================================================================
subroutine dload(f, kstep, kinc, time, noel, npt, layer, kspt, coords, jltyp, sname)
    use gear_contact_mod
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none

    real(real64),     intent(out) :: f
    integer(int32),   intent(in)  :: kstep, kinc, noel, npt, layer, kspt, jltyp
    real(real64),     intent(in)  :: time(2), coords(3)
    character(len=*), intent(in)  :: sname

    real(real64)   :: x, y, z, r2, u_p, alpha_p, t_step, phi_rot
    real(real64)   :: theta, delta_th, dist_flank, u_c, q, g_flank
    real(real64)   :: phase, s_tab, frac, w_inv
    integer(int32) :: idx, mode, istep, m

    f = 0.0_real64

    ! 1. Load type check (jltyp = 0 for surface distributed load PNU) & flank
    !    identification from the surface name (the coast tag contains the drive tag)
    if (jltyp /= 0) return
    if (index(sname, coast_tag) > 0) then
        mode = -1
    else if (index(sname, surf_tag) > 0) then
        mode = 1
    else
        return
    end if

    ! Load case of the current step: only the flank matching its drive mode is loaded
    istep = min(max(kstep, 1), n_cycle)
    if (drive_mode_step(istep) /= mode) return
    m = merge(1, 2, mode == 1)

    ! 2. Safety check: the tables are built by UEXTERNALDB before the first call; this guard
    !    only covers a missing call, and the critical section takes effect when compiled with OpenMP
    if (.not. is_init) then
        !$omp critical
        if (.not. is_init) call init_contact_tables()
        !$omp end critical
    end if

    ! 3. Coordinates & fast axial face-width guard
    x = coords(1)
    y = coords(2)
    z = coords(3)
    if (z < z_lo .or. z > z_hi) return

    ! 4. Exact radial bounding-box filter
    r2 = x * x + y * y
    if (r2 < r_min2 .or. r2 > r_max2) return

    ! 5. Involute roll coordinates
    u_p     = sqrt(r2 - rb2)
    alpha_p = atan(u_p / rb_wheel)

    ! 6. Kinematics: stationary ramp during t_ramp, then rolls for duration t_roll
    t_step  = time(1)
    phi_rot = omega * min(max(t_step - t_ramp, 0.0_real64), t_roll)

    ! 7. Line-of-action angular deviation wrapped to [-pi, pi] (fast branchless anint)
    theta    = atan2(y, x) + real(rot_dir, real64) * phi_rot
    delta_th = real(rot_dir * mode, real64) * (theta_p - theta) + alpha_wt - alpha_p
    delta_th = delta_th - two_pi * anint(delta_th / two_pi)

    ! 8. Contact point along LoA
    dist_flank = rb_wheel * delta_th
    u_c        = u_p + dist_flank
    if (u_c < u_a .or. u_c > u_e) return

    ! 9. Parabolic pressure band across profile (modular pure function)
    q = profile_q(u_p, dist_flank)
    if (abs(q) >= 1.0_real64) return
    g_flank = 1.0_real64 - q * q

    ! 10. Normalised contact length interpolation
    phase = modulo(u_c - u_a + chx(m) * (z - z_mid), p_bt)
    s_tab = phase * inv_p_bt * real(n_tab, real64)
    idx   = min(int(s_tab, int32), n_tab - 1)
    frac  = s_tab - real(idx, real64)
    w_inv = w_tab(idx, m) + frac * (w_tab(idx + 1, m) - w_tab(idx, m))

    ! 11. Pressure evaluation with smoothstep ramp
    f = p_unit * torque_step_nm(istep) * time_ramp(t_step) * flank_weight(u_c) * g_flank * w_inv

end subroutine dload
