!DEC$ FREEFORM
!=======================================================================
!  Abaqus/Standard User Subroutine: DFLUX
!  Model: Circumferential Laser Welding - Moving Keyhole Heat Source:
!         double conical Gaussian (head cone over body cone) with a
!         hemispherical bottom
!         Crown Gear & Differential Carrier Assembly with Ni 6082 filler
!         (bimaterial weld model: THERMAL_RLS/bimaterial_RLS_Thermal.inp)
!  Author: Victor Maia (vhfm08@gmail.com)
!
!  Target shape: "nail head" keyhole weld, a wide and shallow head at the
!  surface (Marangoni flow) over a narrow, nearly parallel body. Both parts
!  are truncated cones of Gaussian radius r0(z), uniform along the depth
!  (TDC model of Wu et al.), so the power of the head stays inside its
!  depth d_up instead of the Gaussian tail of a Goldak head:
!    q = q0_i * exp(-3 r2 / r0(z)2) * g(z)
!  - Head cone (fraction f_up of the absorbed power): 0 <= z <= d_up,
!      r0 = r_top -> r_mid
!  - Body cone (fraction 1 - f_up): d_up < z <= depth_cone,
!      r0 = r_mid -> r_bottom
!  - Bottom (keyhole tip, part of the body power): z > depth_cone,
!      r0 = r_bottom, g = exp(-3 (z - depth_cone)2 / r_cap2). With
!      r_cap = r_bottom the iso-surfaces close as hemispheres.
!  The radius is continuous at z = d_up; the peak flux steps from q0_up to
!  q0_low there, which sharpens the neck between head and body.
!  Each part is normalised over its own volume: the integral of the whole
!  source is Q = eta_absorb * qtot_nominal for any f_up. The penetration of
!  the body, at 5 % of its peak flux, is depth_total = depth_cone + r_cap.
!  The weld pool is also widened by the liquid conductivity of the materials
!  (3 x the solid value above the liquidus).
!  References:
!    - Wu, C.S., Wang, H.G., Zhang, Y.M. (2006): A new heat source model for
!      keyhole plasma arc welding in FEM analysis, Welding Journal, 85(12), 284-291.
!    - Farrokhi, F., Endelt, B., Kristiansen, M. (2019): A numerical model
!      for full and partial penetration hybrid laser welding of thick-section
!      steels, Optics and Laser Technology, 109, 629-642 (combined sources).
!  Conventions & Units:
!    - Global coordinates: mm, s, tonne, mJ, C
!    - Volumetric flux: mJ/(s.mm3) = mW/mm3
!    - Trajectory: Circular weld in the X-Y plane about the Z axis, on the joint
!      radius R = 55 mm (centre of the 0.3 mm filler gap), with the beam axis
!      offset by 0.2 mm towards the steel crown (r = 55.2 mm), as in production.
!    - Ramping: 10 deg ramp-up, 360 deg steady weld, 10 deg ramp-down
!    - Process: 4 kW at 1500 mm/min (25 mm/s); 160 J/mm nominal, 82.6 J/mm
!      absorbed in the 5 mm weld land with eta = 0.516: in a full-penetration
!      keyhole part of the beam leaves through the open root into the gas
!      escape chamber; eta is the share that stays in the land.
!    - First calibration (2026-10-15): heat_source_calibrator.py (conduction
!      surrogate, k = 0.0349 W/mm.K, alpha = 8.0 mm2/s) fitted to the
!      fusion zone of Nogueira (2021), Fig. 6.13(a), CC BY 4.0: 2.42 mm wide
!      at the surface, 5.21 mm deep (full penetration of the weld land into
!      the gas escape chamber). The source power below the 5 mm weld land
!      falls in the chamber and is lost, in the fit and on the mesh: it is
!      the power transmitted through the root of the open keyhole. The
!      shape stays close to the previous one (depth_total 7.1 mm).
!      Widths fitted / measured [mm] at depth 0.15 / 1.12 / 2.24 / 3.51 /
!      4.79: 2.22/2.38, 1.67/1.68, 1.16/0.83, 1.05/0.93, 0.62/1.29; the
!      widening of the root at the chamber is not reproduced. The Abaqus run
!      with that source was too wide in the body (1.32 mm at 2.24 mm, +59 %):
!      the surrogate has no latent heat, no liquid conductivity, no chamber.
!    - Second calibration (2026-10-16, these values): the ratio of the FE
!      solidus width to the surrogate width at each depth, from the first
!      run, kept fixed (space mapping), and the source fitted so that
!      ratio x surrogate matches the macrograph down to 4.2 mm; depth_total
!      fixed at the measured penetration (5.2 mm), eta_absorb free. Widths
!      predicted / measured [mm] at 0.15 / 0.63 / 1.12 / 1.60 / 2.24 / 2.88 /
!      3.51 / 4.15: 2.26/2.38, 2.05/1.92, 1.68/1.68, 1.23/1.29, 0.95/0.83,
!      0.88/0.94, 0.85/0.93, 1.06/1.01. The ratios come from a run with the
!      first source: the next thermal run checks them.
!    - mesh_factor = 0.960: 2054 W of the 2064 W absorbed fall in the weld
!      land; integrated on the DC3D8 mesh (2x2x2 Gauss points), beam at
!      r = 55.2 mm on a node line as at every increment end, the source gives
!      2139 W. Between two node lines it would give 1955 W (-8.6 %): the body
!      radius (0.42-0.45 mm) is smaller than the circumferential element
!      (0.58 mm), keep the increments aligned.
!=======================================================================

module laser_dflux_mod
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none

    ! ------------------------------------------------------------------
    ! 1. PROCESS PARAMETERS
    ! ------------------------------------------------------------------
    real(real64), parameter :: qtot_nominal = 4.00d6        ! Nominal beam power [mW = mJ/s] (4000 W)
    real(real64), parameter :: eta_absorb   = 0.516d0       ! Share of the beam absorbed in the weld land (full penetration) [-]
    real(real64), parameter :: weld_radius  = 55.00d0       ! Circular joint radius [mm] (centre of the filler gap)
    real(real64), parameter :: beam_offset  = 0.20d0        ! Radial offset of the beam axis towards the steel crown [mm]
    real(real64), parameter :: travel_speed = 1500.0d0      ! Travel speed [mm/min] (25 mm/s)
    real(real64), parameter :: z_surface    = 25.00d0       ! Top irradiated surface Z coordinate [mm]
    real(real64), parameter :: cutoff_tol   = 1.0d-8        ! Gaussian tail cut-off, fraction of the peak [-]
    real(real64), parameter :: mesh_factor  = 0.960d0       ! Correction of the power integrated on the mesh [-]

    ! ------------------------------------------------------------------
    ! 2. SOURCE GEOMETRY AND POWER SPLIT
    ! ------------------------------------------------------------------
    real(real64), parameter :: depth_total  = 5.200d0        ! Penetration of the source, depth_cone + r_cap [mm]
    real(real64), parameter :: d_up         = 0.750d0        ! Depth of the head cone [mm]
    real(real64), parameter :: r_top        = 1.016d0        ! Head Gaussian radius at the surface [mm]
    real(real64), parameter :: r_mid        = 0.454d0        ! Gaussian radius at the head / body junction [mm]
    real(real64), parameter :: r_bottom     = 0.420d0        ! Body Gaussian radius at the bottom of the cone [mm] (the body narrows downwards)
    real(real64), parameter :: r_cap        = 0.420d0        ! Depth radius of the bottom (= r_bottom: hemisphere) [mm]
    real(real64), parameter :: f_up         = 0.418d0        ! Fraction of the absorbed power in the head cone [-]

    ! ------------------------------------------------------------------
    ! 3. ANGULAR POWER SCHEDULE PARAMETERS
    ! ------------------------------------------------------------------
    real(real64), parameter :: deg_ramp_up   = 10.00d0      ! Power ramp-up segment [deg]
    real(real64), parameter :: deg_weld      = 360.00d0     ! Full nominal power weld segment [deg]
    real(real64), parameter :: deg_ramp_down = 10.00d0      ! Power ramp-down overlap segment [deg]

    ! Derived constants
    real(real64), parameter :: pi_const   = 3.14159265358979323846d0
    real(real64), parameter :: q_abs      = mesh_factor * eta_absorb * qtot_nominal
    real(real64), parameter :: depth_cone = depth_total - r_cap
    ! Volume integrals of exp(-3 r2 / r0(z)2) * g(z), divided by pi / 3:
    !   truncated cone of height h, radii a -> b:  h (a2 + a b + b2) / 3
    !   bottom: r_bottom2 * integral of exp(-3 u2 / r_cap2) du over u > 0 = r_bottom2 r_cap sqrt(pi / 3) / 2
    real(real64), parameter :: v_up       = d_up * (r_top**2 + r_top * r_mid + r_mid**2) / 3.0d0
    real(real64), parameter :: v_low      = (depth_cone - d_up) * (r_mid**2 + r_mid * r_bottom + r_bottom**2) / 3.0d0
    real(real64), parameter :: v_cap      = r_bottom**2 * r_cap * sqrt(pi_const / 3.0d0) / 2.0d0
    real(real64), parameter :: q0_up      = 3.0d0 * f_up * q_abs / (pi_const * v_up)
    real(real64), parameter :: q0_low     = 3.0d0 * (1.0d0 - f_up) * q_abs / (pi_const * (v_low + v_cap))
    ! Distance, in units of the Gaussian radius, where the flux falls to cutoff_tol
    real(real64), parameter :: cut_factor = sqrt(-log(cutoff_tol) / 3.0d0)

    real(real64), parameter :: v_sec      = travel_speed / 60.0d0
    ! The angular speed is taken on the joint radius, not on the offset beam radius, so the beam
    ! axis stays on a node line at every increment end (dt = 0.6 deg / omega); the travel speed
    ! of the beam itself is then v (weld_radius + beam_offset) / weld_radius, +0.4 %
    real(real64), parameter :: omega      = v_sec / weld_radius
    real(real64), parameter :: r_beam     = weld_radius + beam_offset
    real(real64), parameter :: rad_to_deg = 180.0d0 / pi_const

contains

    !> Power factor of the angular schedule (ramp-up, steady weld, ramp-down, off)
    pure function power_factor(theta_deg) result(amp)
        real(real64), intent(in) :: theta_deg
        real(real64) :: amp

        if (theta_deg < deg_ramp_up) then
            amp = theta_deg / deg_ramp_up
        else if (theta_deg <= (deg_ramp_up + deg_weld)) then
            amp = 1.0d0
        else if (theta_deg <= (deg_ramp_up + deg_weld + deg_ramp_down)) then
            amp = 1.0d0 - (theta_deg - (deg_ramp_up + deg_weld)) / deg_ramp_down
        else
            amp = 0.0d0
        end if
    end function power_factor

    !> Volumetric flux of the double cone + hemispherical bottom at a point
    pure subroutine compute_flux(time_val, coords, q_vol)
        real(real64), intent(in)  :: time_val
        real(real64), intent(in)  :: coords(3)
        real(real64), intent(out) :: q_vol

        real(real64) :: theta_rad, amp, z_local, r0_z, r_dist_sq, g_z, q0_z

        q_vol = 0.0d0

        ! Beam position on the joint and power schedule
        theta_rad = omega * time_val
        amp = power_factor(theta_rad * rad_to_deg)
        if (amp <= 0.0d0) return

        ! Horizontal distance to the beam axis and depth below the irradiated face
        r_dist_sq = (coords(1) - r_beam * cos(theta_rad))**2 + &
                    (coords(2) - r_beam * sin(theta_rad))**2
        z_local = max(z_surface - coords(3), 0.0d0)

        ! Head cone, body cone, hemispherical bottom
        if (z_local <= d_up) then
            r0_z = r_top + (r_mid - r_top) * (z_local / d_up)
            q0_z = q0_up
            g_z  = 1.0d0
        else if (z_local <= depth_cone) then
            r0_z = r_mid + (r_bottom - r_mid) * ((z_local - d_up) / (depth_cone - d_up))
            q0_z = q0_low
            g_z  = 1.0d0
        else
            r0_z = r_bottom
            q0_z = q0_low
            g_z  = 0.0d0
            if (z_local - depth_cone <= cut_factor * r_cap) &
                g_z = exp(-3.0d0 * (z_local - depth_cone)**2 / r_cap**2)
        end if
        if (g_z > 0.0d0 .and. r_dist_sq <= (cut_factor * r0_z)**2) then
            q_vol = amp * q0_z * g_z * exp(-3.0d0 * r_dist_sq / r0_z**2)
        end if

    end subroutine compute_flux

end module laser_dflux_mod

!=======================================================================
!  Abaqus/Standard Entry Point Subroutine
!=======================================================================
subroutine dflux(flux, sol, kstep, kinc, time, noel, npt, coords, &
                 jltyp, temp, press, sname)
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use laser_dflux_mod, only: compute_flux
    implicit none

    ! Abaqus interface arguments
    real(real64), intent(out)     :: flux(2)
    real(real64), intent(in)      :: sol
    integer(int32), intent(in)    :: kstep
    integer(int32), intent(in)    :: kinc
    real(real64), intent(in)      :: time(2)
    integer(int32), intent(in)    :: noel
    integer(int32), intent(in)    :: npt
    real(real64), intent(in)      :: coords(3)
    integer(int32), intent(in)    :: jltyp
    real(real64), intent(in)      :: temp
    real(real64), intent(in)      :: press
    character(len=80), intent(in) :: sname

    real(real64) :: q_val

    ! time(2) is the total time; the welding steps start at total time 0 and follow
    ! each other (filler sector steps), so it is the time since the start of the weld
    call compute_flux(time(2), coords, q_val)

    flux(1) = q_val
    flux(2) = 0.0d0    ! df/dtheta (zero for purely prescribed thermal load)

end subroutine dflux
