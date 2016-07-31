!> This module defines classes for managing Nosé-Hoover Chain thermostats for
!! contant-temperature molecular dynamics.
!! @remark Reference:
!!         Martyna, Klein, and Tuckerman, J. Chem. Phys. 97 (4), 1992.
!! 
!! @author Charlles R. A. Abreu (abreu@eq.ufrj.br)
!! @date Sept 12, 2013
!!
!! @todo Implement reversible integrator described by
!!       Abrams, Rosso, and Tuckerman, J. Chem. Phys. 125, 074115 (2006).
module mThermostat

use mGlobal

implicit none

!===================================================================================================
!> An abstract class for Nosé-Hoover chain integrators.
type, abstract, private :: nhc

  integer,  private :: M, loop = 1, nsy = 1
  real(rb), private :: kT, gKT, dt, dt_2
  real(rb), allocatable :: Q(:), eta(:), p_eta(:), wsy(:), TwoQ(:)

  contains

    !> Performs initialization of the Nosé-Hoover chain.
    !! @param[in] nchain (integer) number of thermostats in the chain.
    !! @param[in] kT (real) target temperature (in energy units).
    !! @param[in] tdamp (real) thermostat period parameter (inverse of
    !!            frequency parameter).
    !! @param[in] dof (integer) number of degrees of freedom in the system.
    procedure :: setup => nhc_setup

    !> Performs allocation of the Suzuki-Yoshida weights:
    procedure :: allocate_suzuky_yoshida => nhc_allocate_suzuki_yoshida

    procedure :: update_direct => nhc_update_direct
    procedure :: update_reverse => nhc_update_reverse

    !> Computes the thermostat-related "energy".
    !!
    !! @details A Nosé-Hoover chain with M thermostats contributes to the
    !!          system's conserved "energy" with:
    !! \f[ H_{therm} = \sum_{i=1}^M \frac{p_{\eta_i}^2}{2Q_i} + gKT\eta_1 +
    !!                 \sum_{i=2}^M kT\eta_i \f]
    procedure :: energy => nhc_energy

    procedure :: positions => nhc_positions
    procedure :: momenta => nhc_momenta

end type

!===================================================================================================
!> A class for Nosé-Hoover chain integration using the iterative,
!! irreversible Velocity Verlet algorithm.
!! @remark Reference:
!!         Martyna, Klein, and Tuckerman, J. Chem. Phys. 97 (4), 1992.
type, extends(nhc) :: nhc_ivv
  real(rb) :: tol = 1.0e-16_rb
  real(rb) :: TwoKE
  real(rb), allocatable :: u(:), um(:), uf(:)
  contains
    ! Deferred bindings:
    procedure :: extra_setup => nhc_ivv_extra_setup
    procedure :: initial_integrate => nhc_ivv_initial_integrate
    procedure :: final_integrate => nhc_ivv_final_integrate
end type nhc_ivv

private :: nhc_setup
private :: nhc_ivv_initial_integrate, nhc_ivv_final_integrate

!===================================================================================================
!> A class for Nosé-Hoover chain integration using a reversible scheme:
type, extends(nhc) :: nhc_rev
  contains
    procedure :: nhc_rev_integrate_atoms, nhc_rev_integrate_rigid_bodies
    generic :: integrate => nhc_rev_integrate_atoms, nhc_rev_integrate_rigid_bodies
    procedure :: integrate_momenta => nhc_rev_integrate_momenta
end type nhc_rev

!===================================================================================================
!> A class for Nosé-Hoover chain integration using a reversible scheme:
type, extends(nhc) :: nhc_mkt
  contains
    procedure :: nhc_mkt_integrate_atoms, nhc_mkt_integrate_rigid_bodies
    generic :: integrate => nhc_mkt_integrate_atoms, nhc_mkt_integrate_rigid_bodies
    procedure :: integrate_momenta => nhc_mkt_integrate_momenta
end type nhc_mkt

contains
  !=================================================================================================
  !                            G E N E R A L     T H E R M O S T A T
  !=================================================================================================
  subroutine nhc_setup( chain, nchain, kT, tdamp, dof, timestep, loop, nsy )
    class(nhc), intent(inout)        :: chain
    integer,    intent(in)           :: nchain
    real(rb),   intent(in)           :: kT, tdamp, timestep
    integer,    intent(in)           :: dof
    integer,    intent(in), optional :: loop, nsy
    chain%M = nchain
    chain%kT = kT
    chain%gKT = dof*kT
    chain%dt = timestep
    chain%dt_2 = 0.5_rb*timestep
    if (present(loop)) chain%loop = loop
    allocate( chain%Q(nchain), chain%TwoQ(chain%M), chain%eta(nchain), chain%p_eta(nchain) )
    chain%Q(1) = chain%gKT*tdamp**2
    chain%Q(2:nchain) = kT*tdamp**2
    chain%TwoQ = two*chain%Q
    chain%eta = zero
    chain%p_eta = zero
    if (present(nsy)) chain%nsy = nsy
    call chain % allocate_suzuky_yoshida()
  end subroutine nhc_setup
  !-------------------------------------------------------------------------------------------------
  subroutine nhc_allocate_suzuki_yoshida( chain )
    class(nhc), intent(inout) :: chain
    allocate( chain%wsy(chain%nsy) )
    select case (chain%nsy)
      case (1)
        chain%wsy = 1.0_rb
      case (3)
        chain%wsy([1,3]) =  1.3512071919596578_rb
        chain%wsy(  2  ) = -1.7024143839193155_rb
      case (7)
        chain%wsy([1,7]) =  0.784513610477560_rb
        chain%wsy([2,6]) =  0.235573213359357_rb
        chain%wsy([3,5]) = -1.177679984178872_rb
        chain%wsy(  4  ) =  1.315186320683910_rb
      case default
        call error( "Suzuki-Yoshida method with", int2str(chain%nsy), "terms is not supported." )
    end select
  end subroutine nhc_allocate_suzuki_yoshida
  !-------------------------------------------------------------------------------------------------
  subroutine nhc_extra_setup( chain )
    class(nhc), intent(inout) :: chain
  end subroutine nhc_extra_setup
  !-------------------------------------------------------------------------------------------------
  subroutine nhc_update_direct( chain, KE, dt )
    class(nhc), intent(inout) :: chain
    real(rb),   intent(in)    :: KE, dt
    integer :: i
    real(rb) :: G
    G = two*KE - chain%gKT
    do i = 1, chain%M-1
      call update_thermostat( chain%eta(i), chain%p_eta(i), chain%TwoQ(i), &
                              G, chain%p_eta(i+1)/chain%Q(i+1), dt )
      G = chain%p_eta(i)**2/chain%Q(i) - chain%kT
    end do
    call update_thermostat( chain%eta(chain%M), chain%p_eta(chain%M), &
                            chain%TwoQ(chain%M), G, 0.0_rb, dt )
  end subroutine nhc_update_direct
  !-------------------------------------------------------------------------------------------------
  subroutine nhc_update_reverse( chain, KE, dt )
    class(nhc), intent(inout) :: chain
    real(rb),   intent(in)    :: KE, dt
    integer :: i
    real(rb) :: V
    V = 0.0_rb
    do i = chain%M, 2, -1
      call update_thermostat( chain%eta(i), chain%p_eta(i), chain%TwoQ(i),       &
                              chain%p_eta(i-1)**2/chain%Q(i-1) - chain%kT, V, dt )
      V = chain%p_eta(i)/chain%Q(i)
    end do
    call update_thermostat( chain%eta(1), chain%p_eta(1), chain%TwoQ(1), two*KE-chain%gKT, V, dt )
  end subroutine nhc_update_reverse
  !-------------------------------------------------------------------------------------------------
  function nhc_energy( chain ) result( H )
    class(nhc), intent(in) :: chain
    real(rb)               :: H
    H = half*sum(chain%p_eta**2/chain%Q) + chain%gKT * chain%eta(1) + &
        chain%kT * sum( chain%eta(2:chain%M) )
  end function nhc_energy
  !-------------------------------------------------------------------------------------------------
  function nhc_positions( chain ) result( eta )
    class(nhc), intent(in) :: chain
    real(rb)               :: eta(chain%M)
    eta = chain % eta
  end function nhc_positions
  !-------------------------------------------------------------------------------------------------
  function nhc_momenta( chain ) result( p_eta )
    class(nhc), intent(in) :: chain
    real(rb)               :: p_eta(chain%M)
    p_eta = chain % p_eta
  end function nhc_momenta
  !=================================================================================================
  !                     I T E R A T I V E    V E L O C I T Y    V E R L E T
  !=================================================================================================
  subroutine nhc_ivv_extra_setup( chain )
    class(nhc_ivv), intent(inout) :: chain
    allocate( chain%u(chain%M), chain%um(chain%M), chain%uf(chain%M) )
    chain%u = zero
  end subroutine nhc_ivv_extra_setup
  !-------------------------------------------------------------------------------------------------
  subroutine nhc_ivv_initial_integrate( chain, Px, Py, Pz, InvMass, Fx, Fy, Fz )
    class(nhc_ivv), intent(inout) :: chain
    real(rb),       intent(inout) :: Px(:), Py(:), Pz(:)
    real(rb),       intent(in)    :: InvMass(:), Fx(:), Fy(:), Fz(:)
    integer  :: j
    real(rb) :: factor, du, force
    factor = one - chain%dt_2 * chain%u(1)
    Px = Px * factor
    Py = Py * factor
    Pz = Pz * factor
    chain%TwoKE = sum(InvMass*(Px**2 + Py**2 + Pz**2))
    force = chain%TwoKE - chain%gKT
    do j = 1, chain%M-1
      du = chain%dt_2 * (force / chain%Q(j) - chain%u(j) * chain%u(j+1))
      chain%um(j) = chain%u(j) + du
      force = chain%Q(j) * chain%u(j)**2 - chain%kT
      chain%u(j) = chain%um(j) + du
    end do
    du = chain%dt_2 * force / chain%Q(chain%M)
    chain%um(chain%M) = chain%u(chain%M) + du
    chain%u(chain%M) = chain%um(chain%M) + du
    chain%eta = chain%eta + chain%dt * chain%um
    Px = Px + chain%dt_2*Fx
    Py = Py + chain%dt_2*Fy
    Pz = Pz + chain%dt_2*Fz
  end subroutine nhc_ivv_initial_integrate
  !-------------------------------------------------------------------------------------------------
  subroutine nhc_ivv_final_integrate( chain, Px, Py, Pz, InvMass, Fx, Fy, Fz )
    class(nhc_ivv), intent(inout) :: chain
    real(rb),       intent(inout) :: Px(:), Py(:), Pz(:)
    real(rb),       intent(in)    :: InvMass(:), Fx(:), Fy(:), Fz(:)
    integer  :: j
    real(rb) :: factor, du, force
    logical  :: not_converged
    not_converged = .true.
    do while (not_converged)
      chain%uf = chain%u
      force = chain%TwoKE / (one + chain%dt_2 * chain%uf(1))**2 - chain%gKT
      do j = 1, chain%M-1
        du = chain%dt_2 * (force / chain%Q(j) - chain%uf(j) * chain%uf(j+1))
        chain%u(j) = chain%um(j) + du
        force = chain%Q(j) * chain%uf(j)**2 - chain%kT
      end do
      chain%u(chain%M) = chain%um(chain%M) + chain%dt_2 / chain%Q(chain%m) * force
      not_converged = any(abs(chain%u - chain%uf) > chain%tol)
    end do
    factor = one/(one + chain%dt_2 * chain%u(1))
    chain % p_eta = chain%Q * chain%u
    Px = Px*factor + chain%dt_2*Fx
    Py = Py*factor + chain%dt_2*Fy
    Pz = Pz*factor + chain%dt_2*Fz
  end subroutine nhc_ivv_final_integrate
  !=================================================================================================
  !                        R E V E R S I B L E     I N T E G R A T O R
  !=================================================================================================
  subroutine nhc_rev_integrate_atoms( chain, Px, Py, Pz, InvMass, Fx, Fy, Fz )
    class(nhc_rev), intent(inout) :: chain
    real(rb),       intent(inout) :: Px(:), Py(:), Pz(:)
    real(rb),       intent(in)    :: InvMass(:), Fx(:), Fy(:), Fz(:)
    integer :: k, m
    real(rb) :: V, alpha, beta, dt_loop, delta_t
    dt_loop = half*chain%dt_2/chain%loop
    do k = 1, chain%loop
      do m = 1, chain%nsy
        delta_t = dt_loop*chain%wsy(m)
        call chain % update_reverse( sum(InvMass*(Px**2 + Py**2 + Pz**2)), delta_t )
        V = chain%p_eta(1)/chain%Q(1)
        alpha = phi(V*two*delta_t)*two*delta_t
        beta = one - V*alpha
        Px = beta*Px + alpha*Fx
        Py = beta*Py + alpha*Fy
        Pz = beta*Pz + alpha*Fz
        call chain % update_direct( sum(InvMass*(Px**2 + Py**2 + Pz**2)), delta_t )
      end do ! m
    end do ! k
  end subroutine nhc_rev_integrate_atoms
  !-------------------------------------------------------------------------------------------------
  subroutine nhc_rev_integrate_rigid_bodies( chain, body )
    use mRigidBody
    class(nhc_rev),   intent(inout) :: chain
    type(tRigidBody), intent(inout) :: body(:)
    integer :: k, i, m, N, l
    real(rb) :: V, alpha, beta, dt_loop, delta_t
    N = size(body) 
    dt_loop = half*chain%dt_2/chain%loop  
    do k = 1, chain%loop
      do m = 1, chain%nsy
        delta_t = dt_loop*chain%wsy(m)
        call chain % update_reverse( sum(body%kinetic_energy()), delta_t )
        do l = 1, chain%M, 1
	 chain%eta(l) = chain%eta(l) + (chain%p_eta(l)/chain%Q(l))*delta_t*2.0
	end do
        call chain % update_direct( sum(body%kinetic_energy()), delta_t )
      end do ! m
    end do ! k
  end subroutine nhc_rev_integrate_rigid_bodies
  !-------------------------------------------------------------------------------------------------
  subroutine nhc_rev_integrate_momenta( chain, body, delta_t )
    use mRigidBody
    class(nhc_rev),   intent(inout) :: chain
    type(tRigidBody), intent(inout) :: body(:)
    real(rb),          intent(in)    :: delta_t
    integer :: i
    real(rb) :: V, alpha, beta
    V = chain%p_eta(1)/chain%Q(1)
    alpha = phi(V*delta_t)*delta_t
    beta = one - V*alpha
    forall (i=1:size(body))
      body(i)%Pcm = beta*body(i)%Pcm + alpha*body(i)%F
      body(i)%Pi = beta*body(i)%Pi + alpha*body(i)%Two_C_Tau
    end forall
  end subroutine nhc_rev_integrate_momenta
  !=================================================================================================
  !                               M K T      I N T E G R A T O R
  !=================================================================================================
  subroutine nhc_mkt_integrate_atoms( chain, Px, Py, Pz, InvMass, Fx, Fy, Fz )
    class(nhc_mkt), intent(inout) :: chain
    real(rb),       intent(inout) :: Px(:), Py(:), Pz(:)
    real(rb),       intent(in)    :: InvMass(:), Fx(:), Fy(:), Fz(:)
    integer :: k, m, l
    real(rb) :: V, alpha, beta, dt_loop, delta_t
    dt_loop = half*chain%dt_2/chain%loop
    do k = 1, chain%loop
      do m = 1, chain%nsy
        delta_t = dt_loop*chain%wsy(m)
        call chain % update_reverse( sum(InvMass*(Px**2 + Py**2 + Pz**2)), delta_t )
	do l = 1, chain%M, 1
	 chain%eta(l) = chain%eta(l) + (chain%p_eta(l)/chain%Q(l))*delta_t*2.0
	end do
	V = chain%p_eta(1)/chain%Q(1)
        alpha = exp(-V*two*delta_t)            
        Px = alpha*Px 
        Py = alpha*Py 
        Pz = alpha*Pz 
        call chain % update_direct( sum(InvMass*(Px**2 + Py**2 + Pz**2)), delta_t )
      end do ! m
    end do ! k


  end subroutine nhc_mkt_integrate_atoms
  !-------------------------------------------------------------------------------------------------
  subroutine nhc_mkt_integrate_rigid_bodies( chain, body )
    use mRigidBody
    class(nhc_mkt),   intent(inout) :: chain
    type(tRigidBody), intent(inout) :: body(:)
    integer :: k, i, m, N, l
    real(rb) :: V, alpha, beta, dt_loop, delta_t
    N = size(body)
    dt_loop = half*chain%dt_2/chain%loop
    do k = 1, chain%loop
      do m = 1, chain%nsy
        delta_t = dt_loop*chain%wsy(m)
        call chain % update_reverse( sum(body%kinetic_energy()), delta_t )
	do l = 1, chain%M, 1
	 chain%eta(l) = chain%eta(l) + (chain%p_eta(l)/chain%Q(l))*delta_t*2.0
	end do
        V = chain%p_eta(1)/chain%Q(1)
        alpha = exp(-V*two*delta_t)        
        forall (i=1:N)
          body(i)%Pcm = alpha*body(i)%Pcm 
          body(i)%Pi = alpha*body(i)%Pi 
        end forall
        call chain % update_direct( sum(body%kinetic_energy()), delta_t )
      end do ! m
    end do ! k


  end subroutine nhc_mkt_integrate_rigid_bodies
  !-------------------------------------------------------------------------------------------------
  subroutine nhc_mkt_integrate_momenta( chain, body, delta_t )
    use mRigidBody
    class(nhc_mkt),   intent(inout) :: chain
    type(tRigidBody), intent(inout) :: body(:)
    real(rb),          intent(in)    :: delta_t
    integer :: i

    forall (i=1:size(body))
    body(i)%Pcm = body(i)%Pcm + delta_t*body(i)%F
    body(i)%Pi = body(i)%Pi + delta_t*body(i)%Two_C_Tau
    end forall
  end subroutine nhc_mkt_integrate_momenta

  !=================================================================================================
  !                             A U X I L I A R Y     F U N C T I O N S
  !=================================================================================================
  subroutine update_thermostat( eta, p_eta, TwoQ, G, V, delta_t )
    real(rb), intent(inout) :: eta, p_eta
    real(rb), intent(in)    :: TwoQ, G, V, delta_t
    p_eta = p_eta + (G - V*p_eta)*phi(V*delta_t)*delta_t
  end subroutine update_thermostat
  !-------------------------------------------------------------------------------------------------
  elemental real(rb) function phi( x )
    real(rb), intent(in) :: x
    if (abs(x) > 1e-4_rb ) then
      phi = (one - exp(-x))/x
    else
      phi = one + half*x*(third*x*(one - fourth*x) - one)
    end if
  end function phi
  !-------------------------------------------------------------------------------------------------
end module mThermostat
