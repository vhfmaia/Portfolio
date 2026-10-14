!DEC$ FREEFORM
!=======================================================================
!  Abaqus/Standard User Subroutine: UMATHT
!  Model: Hydrogen transport for cold (delayed) cracking of welds, solved
!         with the heat transfer analogy: the nodal temperature DOF is the
!         lattice hydrogen concentration C_L.
!  Author: Victor Maia (vhfm08@gmail.com)
!
!  Balance (Sofronis & McMeeking 1989, Oriani 1970 local equilibrium):
!      d(C_L + C_T)/dt + div J = 0
!      J   = -D_L grad C_L + D_L C_L V_H / (R T) grad sigma_h
!      C_T = N_T theta_T,  theta_T / (1 - theta_T) = K_T theta_L / (1 - theta_L)
!      K_T = exp(W_B / (R T)),  theta_L = C_L / N_L
!  Heat transfer analogy (density = 1):
!      U    = C_L + C_T            (internal energy per unit mass)
!      FLUX = J                    (heat flux vector, -k grad T in Abaqus)
!      DUDT = 1 + dC_T/dC_L        (= D_L / D_eff)
!  The increment of U is computed from the exact trapped amounts at the start
!  and at the end of the increment, so the hydrogen balance holds for any trap
!  occupancy, and traps created by plastic strain or released by heating act
!  as sinks and sources.
!
!  Lattice diffusivity: Arrhenius per phase, mixed linearly with the austenite
!  fraction:  D_L = (1 - f_g) D0_a exp(-Q_a/RT) + f_g D0_g exp(-Q_g/RT).
!  The solubility jump between austenite and ferrite / martensite is not
!  modelled (single N_L); see the README.
!  Trap density: constant N_T, plus optionally the dislocation traps of
!  Kumnick & Johnson (1980) as a function of the equivalent plastic strain:
!      log10 N_T [sites/m3] = 23.26 - 2.33 exp(-5.5 eps_p)
!
!  Units: mm, s, MPa, mol, K for the Arrhenius terms (C for field variable 1).
!    C_L, C_T, N_L, N_T in mol/m3 (any one unit for all four);
!    D in mm2/s; V_H in mm3/mol; W_B and Q in J/mol.
!
!  PROPS (*USER MATERIAL, TYPE=THERMAL, CONSTANTS=9):
!    1  D0_a   pre-exponential diffusivity, ferrite / martensite [mm2/s]
!    2  Q_a    activation energy, ferrite / martensite [J/mol]
!    3  D0_g   pre-exponential diffusivity, austenite [mm2/s]
!    4  Q_g    activation energy, austenite [J/mol]
!    5  N_L    lattice site density [mol/m3]
!    6  V_H    partial molar volume of hydrogen [mm3/mol]
!    7  W_B    trap binding energy [J/mol]
!    8  N_T    constant trap density [mol/m3]
!    9  KJ     1 = add the Kumnick-Johnson traps of field variable 6, 0 = off
!  Field variables (all six must be defined):
!    1  temperature [C]
!    2  austenite fraction [-]
!    3-5  gradient of the hydrostatic stress, x y z [MPa/mm]
!    6  equivalent plastic strain [-]
!  STATEV (*DEPVAR 4): 1 C_L, 2 C_T, 3 D_L [mm2/s], 4 N_T [mol/m3]
!  Density: *DENSITY 1.0. No *CONDUCTIVITY or *SPECIFIC HEAT.
!
!  References:
!    - Oriani, R.A. (1970): The diffusion and trapping of hydrogen in steel,
!      Acta Metallurgica, 18, 147-157.
!    - Sofronis, P., McMeeking, R.M. (1989): Numerical analysis of hydrogen
!      transport near a blunting crack tip, J. Mech. Phys. Solids, 37, 317-350.
!    - Kumnick, A.J., Johnson, H.H. (1980): Deep trapping states for hydrogen
!      in deformed iron, Acta Metallurgica, 28, 33-39.
!    - Diaz, A., Alegre, J.M., Cuesta, I.I. (2016): Coupled hydrogen diffusion
!      simulation using a heat transfer analogy, Int. J. Mech. Sci., 115-116,
!      360-369.
!=======================================================================

module hydrogen_umatht_mod
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none

    real(real64), parameter :: r_gas_j    = 8.314462618d0     ! Gas constant [J/(mol.K)]
    real(real64), parameter :: r_gas_mj   = 8314.462618d0     ! Gas constant [mJ/(mol.K)] = [N.mm/(mol.K)]
    real(real64), parameter :: t_abs_zero = 273.15d0          ! 0 C in K
    real(real64), parameter :: n_avogadro = 6.02214076d23     ! [1/mol]

contains

    !> Trap density [mol/m3]: constant part plus Kumnick-Johnson dislocation traps
    pure function trap_density(nt_const, kj_flag, peeq) result(n_t)
        real(real64), intent(in) :: nt_const, kj_flag, peeq
        real(real64) :: n_t

        n_t = nt_const
        if (kj_flag > 0.5d0) then
            n_t = n_t + 10.0d0**(23.26d0 - 2.33d0 * exp(-5.5d0 * max(peeq, 0.0d0))) / n_avogadro
        end if
    end function trap_density

    !> Trapped concentration in Oriani equilibrium with the lattice
    pure function trapped(c_l, n_l, n_t, k_t) result(c_t)
        real(real64), intent(in) :: c_l, n_l, n_t, k_t
        real(real64) :: c_t, theta_l

        ! Negative C_L (small oscillations of the solution) are not trapped
        theta_l = max(c_l, 0.0d0) / n_l
        c_t = n_t * k_t * theta_l / (1.0d0 - theta_l + k_t * theta_l)
    end function trapped

    !> Derivative dC_T/dC_L of the Oriani equilibrium
    pure function dtrapped(c_l, n_l, n_t, k_t) result(dct)
        real(real64), intent(in) :: c_l, n_l, n_t, k_t
        real(real64) :: dct, theta_l

        theta_l = max(c_l, 0.0d0) / n_l
        dct = n_t / n_l * k_t / (1.0d0 - theta_l + k_t * theta_l)**2
    end function dtrapped

end module hydrogen_umatht_mod

!=======================================================================
!  Abaqus Entry Point Subroutine
!=======================================================================
subroutine umatht(u, dudt, dudg, flux, dfdt, dfdg, statev, temp, dtemp, &
                  dtemdx, time, dtime, predef, dpred, cmname, ntgrd, nstatv, &
                  props, nprops, coords, pnewdt, noel, npt, layer, kspt, &
                  kstep, kinc)
    use, intrinsic :: iso_fortran_env, only: real64
    use hydrogen_umatht_mod
    implicit none

    ! Abaqus interface arguments
    character(len=80), intent(in) :: cmname
    integer, intent(in)           :: ntgrd, nstatv, nprops, noel, npt, layer, kspt, kstep, kinc
    real(real64), intent(inout)   :: u
    real(real64), intent(out)     :: dudt, dudg(ntgrd), flux(ntgrd), dfdt(ntgrd), dfdg(ntgrd, ntgrd)
    real(real64), intent(inout)   :: statev(nstatv), pnewdt
    real(real64), intent(in)      :: temp, dtemp, dtemdx(ntgrd), time(2), dtime
    real(real64), intent(in)      :: predef(*), dpred(*), props(nprops), coords(3)

    real(real64) :: t_old, t_new, f_gam, peeq_old, peeq_new
    real(real64) :: d_alpha, d_gamma, d_l, n_l, v_h, w_b
    real(real64) :: k_old, k_new, nt_old, nt_new, c_old, c_new, drift
    real(real64) :: grad_sh(3)
    integer :: i

    ! Field variables at the start (old) and at the end (new) of the increment
    t_old    = predef(1) + t_abs_zero
    t_new    = predef(1) + dpred(1) + t_abs_zero
    f_gam    = min(max(predef(2) + dpred(2), 0.0d0), 1.0d0)
    grad_sh  = predef(3:5) + dpred(3:5)
    peeq_old = predef(6)
    peeq_new = predef(6) + dpred(6)

    ! Material constants
    n_l = props(5)
    v_h = props(6)
    w_b = props(7)

    ! Lattice diffusivity, end of increment (backward Euler)
    d_alpha = props(1) * exp(-props(2) / (r_gas_j * t_new))
    d_gamma = props(3) * exp(-props(4) / (r_gas_j * t_new))
    d_l     = (1.0d0 - f_gam) * d_alpha + f_gam * d_gamma

    ! Trap equilibrium constant and trap density, start and end of increment
    k_old  = exp(w_b / (r_gas_j * t_old))
    k_new  = exp(w_b / (r_gas_j * t_new))
    nt_old = trap_density(props(8), props(9), peeq_old)
    nt_new = trap_density(props(8), props(9), peeq_new)

    ! Total hydrogen (lattice + trapped) as the internal energy
    c_old = temp
    c_new = temp + dtemp
    u     = u + (c_new + trapped(c_new, n_l, nt_new, k_new)) &
              - (c_old + trapped(c_old, n_l, nt_old, k_old))
    dudt  = 1.0d0 + dtrapped(c_new, n_l, nt_new, k_new)
    dudg  = 0.0d0

    ! Flux: Fick + drift up the hydrostatic stress gradient
    drift = v_h / (r_gas_mj * t_new)
    dfdg  = 0.0d0
    do i = 1, ntgrd
        flux(i)    = -d_l * dtemdx(i) + d_l * c_new * drift * grad_sh(i)
        dfdt(i)    = d_l * drift * grad_sh(i)
        dfdg(i, i) = -d_l
    end do

    ! Output
    if (nstatv >= 4) then
        statev(1) = c_new
        statev(2) = trapped(c_new, n_l, nt_new, k_new)
        statev(3) = d_l
        statev(4) = nt_new
    end if

end subroutine umatht
