module mRigidBody

use mGlobal
use mMath

implicit none

#ifndef nsy
#  define nsy 1
#endif

#if nsy == 1
  real(rb), parameter, private :: wsy(1) = [ 1.0_rb ]
#elif nsy == 3
  real(rb), parameter, private :: wsy(3) = [ 1.3512071919596578_rb, &
                                            -1.7024143839193155_rb, &
                                             1.3512071919596578_rb  ]
#elif nsy == 7
  real(rb), parameter, private :: wsy(7) = [ 0.784513610477560_rb,  &
                                             0.235573213359357_rb,  &
                                            -1.177679984178872_rb,  &
                                             1.315186320683910_rb,  &
                                            -1.177679984178872_rb,  &
                                             0.235573213359357_rb,  &
                                             0.784513610477560_rb   ]
#endif

type tRigidBody
  integer :: NP
  real(rb) :: Mass, Rcm(3), Pcm(3), F(3), Two_C_Tau(0:3)
  real(rb), allocatable :: M(:), delta(:,:), d(:,:)
  real(rb) :: MoI(3)
  real(rb) :: Q(0:3), Pi(0:3)
  real(rb) :: InvMass, InvMoI(3)
  contains
    procedure :: setup => tRigidBody_setup
    procedure :: assign_momenta => tRigidBody_assign_momenta

    procedure :: translational_energy => tRigidBody_translational_energy
    procedure :: rotational_energy => tRigidBody_rotational_energy
    procedure :: kinetic_energy => tRigidBody_kinetic_energy

    procedure :: angular_velocity => tRigidBody_angular_velocity
    procedure :: angular_momentum => tRigidBody_angular_momentum
    procedure :: space_fixed_angular_momentum => tRigidBody_space_fixed_angular_momentum
    procedure :: maximum_distance => tRigidBody_maximum_distance

    procedure :: force_and_torque => tRigidBody_force_and_torque

    procedure :: boost => tRigidBody_boost
    procedure :: rotate => tRigidBody_rotate
    procedure :: displace => tRigidBody_displace

    procedure :: get_positions => tRigidBody_get_positions
    procedure :: get_momenta => tRigidBody_get_momenta
    procedure :: set_momenta => tRigidBody_set_momenta

    procedure :: quaternion_unity_sqdev => tRigidBody_quaternion_unity_sqdev

end type tRigidBody

private :: uniaxial_rotation, quaternion, matrix_B, matrix_C, add_inertia !, rotation_matrix

contains
  !-------------------------------------------------------------------------------------------------
  subroutine tRigidBody_setup( me, Mass, Rx, Ry, Rz )
    class(tRigidBody), intent(inout) :: me
    real(rb),          intent(in)    :: Mass(:), Rx(:), Ry(:), Rz(:)
    integer :: i
    real(rb) :: inertia(3,3), A(3,3), R(3,size(Mass))
    me%NP = size(Mass)
    allocate( me%M(me%NP), me%delta(3,me%NP), me%d(3,me%NP) )
    ! Save particle masses and positions:
    me%M = Mass
    R(1,:) = Rx
    R(2,:) = Ry
    R(3,:) = Rz
    ! Compute total mass and center-of-mass position:
    me%Mass = sum(me%M)
    me%InvMass = one/me%Mass
    forall (i=1:3) me%Rcm(i) = sum(me%M*R(i,:))*me%InvMass
    ! Compute inertia tensor:
    forall (i=1:3) me%delta(i,:) = R(i,:) - me%Rcm(i)
    inertia = zero
    do i = 1, me%NP
      call add_inertia( inertia, me%M(i), me%delta(:,i) )
    end do
    ! Diagonalize the inertia tensor:
    me%MoI = eigenvalues( inertia )
    me%InvMoI = one/me%MoI
    A = transpose(eigenvectors( inertia, me%MoI ))
    ! Compute quaternion:
    me%Q = quaternion( A )
    ! Calculate position in the body-fixed frame:
    forall (i=1:me%NP) me%d(:,i) = matmul( A, me%delta(:,i) )
    ! Zero center-of-mass and quaternion-conjugated momenta:
    me%Pcm = zero
    me%Pi = zero
  end subroutine tRigidBody_setup
  !-------------------------------------------------------------------------------------------------
  ! This subroutine assigns the linear momentum and quaternion-conjugated momentum of a rigid body:
  subroutine tRigidBody_assign_momenta( me, kBT, random, rotation )
    use mRandom
    class(tRigidBody), intent(inout) :: me
    real(rb),          intent(in)    :: kBT
    class(i32rng),     intent(inout) :: random
    logical,           intent(in)    :: rotation
    integer :: k
    real(rb) :: factor, vcm(3), omega(3)
    if (rotation) then
      factor = sqrt(kBT*me%InvMass)
    else
      factor = sqrt(two*kBT*me%InvMass)
    end if
    do k = 1, 3
      vcm(k) = factor*random%normal()
    end do
    me%Pcm = me%Mass*vcm
    if (rotation) then
      do k = 1, 3
        omega(k) = sqrt(kBT*me%InvMoI(k))*random%normal()
      end do
      me%Pi = two*matmul( matrix_B(me%Q), me%MoI*omega )
    else
      me%Pi = zero
    end if
  end subroutine tRigidBody_assign_momenta
  !-------------------------------------------------------------------------------------------------
  ! This function returns the translational kinetic energy of a rigid body:
  elemental function tRigidBody_translational_energy( me ) result( KEt )
    class(tRigidBody), intent(in) :: me
    real(rb)                      :: KEt
    KEt = half*me%InvMass*sum(me%Pcm**2)
  end function tRigidBody_translational_energy
  !-------------------------------------------------------------------------------------------------
  ! This function returns the rotational kinetic energy of a rigid body:
  elemental function tRigidBody_rotational_energy( me, dir ) result( KEr )
    class(tRigidBody), intent(in) :: me
    real(rb)                      :: KEr
    integer, intent(in), optional :: dir
    real(rb) :: omega(3)
    omega = me%angular_momentum()
    if (present(dir)) then
      KEr = half*me%InvMoI(dir)*omega(dir)**2
    else
      KEr = half*sum(me%InvMoI*omega**2)
    end if
  end function tRigidBody_rotational_energy
  !-------------------------------------------------------------------------------------------------
  ! This function returns the (translational + rotational) kinetic energy of a rigid body:
  elemental function tRigidBody_kinetic_energy( me ) result( KE )
    class(tRigidBody), intent(in) :: me
    real(rb)                      :: KE
    KE = half*( me%InvMass*sum(me%Pcm**2) + sum(me%InvMoI * me%angular_momentum()**2) )
  end function tRigidBody_kinetic_energy
  !-------------------------------------------------------------------------------------------------
  ! This function returns the body-fixed angular velocity of a rigid body:
  pure function tRigidBody_angular_velocity( me ) result( omega )
    class(tRigidBody), intent(in) :: me
    real(rb)                      :: omega(3)
    omega = half * me%InvMoI * matmul(transpose(matrix_B(me%Q)),me%Pi)
  end function tRigidBody_angular_velocity
  !-------------------------------------------------------------------------------------------------
  ! This function returns the body-fixed angular momentum of a rigid body:
  pure function tRigidBody_angular_momentum( me ) result( Iw )
    class(tRigidBody), intent(in) :: me
    real(rb)                      :: Iw(3)
    Iw = half * matmul(transpose(matrix_B(me%Q)),me%Pi)
  end function tRigidBody_angular_momentum
  !-------------------------------------------------------------------------------------------------
  ! This function returns the space-fixed angular momentum of a rigid body:
  pure function tRigidBody_space_fixed_angular_momentum( me ) result( L )
    class(tRigidBody), intent(in) :: me
    real(rb)                      :: L(3)
    L = half * matmul(transpose(matrix_C(me%Q)),me%Pi)
  end function tRigidBody_space_fixed_angular_momentum
  !-------------------------------------------------------------------------------------------------
  ! This function returns the maximum influence length of a rigid body:
  elemental function tRigidBody_maximum_distance( me ) result( L )
    class(tRigidBody), intent(in) :: me
    real(rb)                      :: L
    integer :: i, j
    real(rb) :: Lsq
    Lsq = zero
    do i = 1, me%NP-1
      do j = i+1, me%NP
        Lsq = max(Lsq,sum((me%d(:,i) - me%d(:,j))**2))
      end do
    end do
    L = sqrt(Lsq)
  end function tRigidBody_maximum_distance
  !-------------------------------------------------------------------------------------------------
  ! This subroutine zeroes the total linear momentum and the total angular momentum of a set of
  ! rigid bodies and rescales each body's linear and angular momentum so that the total kinetic
  ! energy equals a specified value. Note: current version considers that bodies do not rotate.
  subroutine setup_microcanonical_ensemble( body, KE )
    class(tRigidBody), intent(inout) :: body(:)
    real(rb),          intent(in)    :: KE
    integer :: NB, i
    real(rb) :: Mass, P(3), Rcm(3), L(3), Inertia(3,3), Omega(3), delta(3), factor
    ! Compute total mass, center-of-mass location, and total linear momemtum:
    NB = size(body)
    Mass = zero
    Rcm = zero
    P = zero
    do i = 1, NB
      Mass = Mass + body(i)%Mass
      Rcm = Rcm + body(i)%Mass * body(i)%Rcm
      P = P + body(i)%Pcm
    end do
    Rcm = Rcm/Mass
    ! Zero the total linear momentum:
    do i = 1, NB
      body(i)%Pcm = body(i)%Pcm - (body(i)%Mass/Mass)*P
    end do
    ! Compute total angular momentum and inertia tensor with respect to the center of mass:
    L = zero
    inertia = zero
    do i = 1, NB
      delta = body(i)%Rcm - Rcm
      L = L + cross_product( delta, body(i)%Pcm )
      call add_inertia( inertia, body(i)%Mass, delta )
    end do
    ! Compute angular velocity of the system:
    omega = matmul( inv3x3(Inertia), L )
    ! Zero angular momentum:
    do i = 1, NB
      body(i)%Pcm = body(i)%Pcm - body(i)%Mass*cross_product( omega, body(i)%Rcm - Rcm )
    end do
    ! Rescale velocities to adjust kinetic energy:
    factor = sqrt(KE/sum(body%kinetic_energy()))
    forall (i=1:NB)
      body(i)%Pcm = factor*body(i)%Pcm
      body(i)%Pi = factor*body(i)%Pi
    end forall
  end subroutine setup_microcanonical_ensemble
  !-------------------------------------------------------------------------------------------------
  subroutine tRigidBody_force_and_torque( me, Fx, Fy, Fz )
    class(tRigidBody), intent(inout) :: me
    real(rb),          intent(in)    :: Fx(:), Fy(:), Fz(:)
    integer :: i
    real(rb) :: Fi(3), torque(3)
    me%F = zero
    torque = zero
    do i = 1, me%NP
      Fi = [Fx(i), Fy(i), Fz(i)]
      me%F = me%F + Fi
      torque = torque + cross_product( me%delta(:,i), Fi )
    end do
    me%Two_C_Tau = two*matmul( matrix_C( me%Q ), torque )
  end subroutine tRigidBody_force_and_torque
  !-------------------------------------------------------------------------------------------------
  elemental subroutine tRigidBody_boost( me, dt )
    class(tRigidBody), intent(inout) :: me
    real(rb),          intent(in)    :: dt
    me%Pcm = me%Pcm + dt*me%F
    me%Pi = me%Pi + dt*me%Two_C_Tau
  end subroutine tRigidBody_boost
  !-------------------------------------------------------------------------------------------------

#if nsy == 1

  elemental subroutine tRigidBody_rotate( me, dt ) !_total )
    class(tRigidBody), intent(inout) :: me
    real(rb),          intent(in)    :: dt !_total
    real(rb) :: half_dt
    half_dt = half*dt
    call uniaxial_rotation( 3, half_dt, me%MoI(3), me%Q, me%Pi )
    call uniaxial_rotation( 2, half_dt, me%MoI(2), me%Q, me%Pi )
    call uniaxial_rotation( 1, dt, me%MoI(1), me%Q, me%Pi )
    call uniaxial_rotation( 2, half_dt, me%MoI(2), me%Q, me%Pi )
    call uniaxial_rotation( 3, half_dt, me%MoI(3), me%Q, me%Pi )
  end subroutine tRigidBody_rotate

#else

  elemental subroutine tRigidBody_rotate( me, dt_total )
    class(tRigidBody), intent(inout) :: me
    real(rb),          intent(in)    :: dt_total
    real(rb) :: half_dt, dt
    integer :: k
    do k = 1, nsy
      dt = dt_total*wsy(k)
      half_dt = half*dt
      call uniaxial_rotation( 3, half_dt, me%MoI(3), me%Q, me%Pi )
      call uniaxial_rotation( 2, half_dt, me%MoI(2), me%Q, me%Pi )
      call uniaxial_rotation( 1, dt, me%MoI(1), me%Q, me%Pi )
      call uniaxial_rotation( 2, half_dt, me%MoI(2), me%Q, me%Pi )
      call uniaxial_rotation( 3, half_dt, me%MoI(3), me%Q, me%Pi )
    end do
    !me%Q = me%Q/sqrt(sum(me%Q**2))
  end subroutine tRigidBody_rotate

#endif

  !-------------------------------------------------------------------------------------------------
  elemental subroutine tRigidBody_displace( me, dt )
    class(tRigidBody), intent(inout) :: me
    real(rb),          intent(in)    :: dt
    me%Rcm = me%Rcm + dt*me%InvMass*me%Pcm
  end subroutine tRigidBody_displace
  !-------------------------------------------------------------------------------------------------
  subroutine tRigidBody_get_positions( me, Rx, Ry, Rz )
    class(tRigidBody), intent(inout) :: me
    real(rb),          intent(out)   :: Rx(me%NP), Ry(me%NP), Rz(me%NP)
    integer :: i
    real(rb) :: At(3,3)
    At = inverse_rotation_matrix( me%Q )
    forall (i=1:me%NP) me%delta(:,i) = matmul(At,me%d(:,i))
    Rx = me%Rcm(1) + me%delta(1,:)
    Ry = me%Rcm(2) + me%delta(2,:)
    Rz = me%Rcm(3) + me%delta(3,:)
  end subroutine tRigidBody_get_positions
  !-------------------------------------------------------------------------------------------------
  subroutine tRigidBody_get_momenta( me, Px, Py, Pz )
    class(tRigidBody), intent(inout) :: me
    real(rb),          intent(out)   :: Px(me%NP), Py(me%NP), Pz(me%NP)
    integer :: i
    real(rb) :: At(3,3), omega(3), P(3)
    At = inverse_rotation_matrix( me%Q )
    omega = me % angular_velocity()
    do i = 1, me%NP
      P = me%M(i)*(me%Pcm/me%Mass + matmul(At,cross_product(omega,me%d(:,i))))
      Px(i) = P(1)
      Py(i) = P(2)
      Pz(i) = P(3)
    end do
  end subroutine tRigidBody_get_momenta
  !-------------------------------------------------------------------------------------------------
  subroutine tRigidBody_set_momenta( me, Px, Py, Pz )
    class(tRigidBody), intent(inout) :: me
    real(rb),          intent(in)    :: Px(me%NP), Py(me%NP), Pz(me%NP)
    integer :: i
    real(rb) :: P(3), L(3)
    me%Pcm = zero
    L = zero
    do i = 1, me%NP
      P = [Px(i),Py(i),Pz(i)]
      me%Pcm = me%Pcm + P
      L = L + cross_product(me%delta(:,i),P)
    end do
    me%Pi = two*matmul(matrix_C( me%Q ),L)
  end subroutine tRigidBody_set_momenta
  !-------------------------------------------------------------------------------------------------
  elemental function tRigidBody_quaternion_unity_sqdev( me ) result( sqdev )
    class(tRigidBody), intent(in) :: me
    real(rb)                      :: sqdev
    sqdev = (sqrt(sum(me%Q**2)) - one)**2
  end function tRigidBody_quaternion_unity_sqdev
  !-------------------------------------------------------------------------------------------------
  pure subroutine uniaxial_rotation( k, dt, MoI, Q, Pi )
    integer,  intent(in)    :: k
    real(rb), intent(in)    :: dt, MoI
    real(rb), intent(inout) :: Q(4), Pi(4)
    real(rb) :: BkQ(4), omega_dt_over_2, vsin, vcos
    BkQ = Permutation(Q,k)
    omega_dt_over_2 = dt*sum(Pi*BkQ)/(4.0_rb*MoI)
    vsin = sin(omega_dt_over_2)
    vcos = cos(omega_dt_over_2)
    Q = vcos*Q + vsin*BkQ
    Pi = vcos*Pi + vsin*Permutation(Pi,k)
    contains
      pure function Permutation( Q, k ) result( BkQ )
        real(rb), intent(in) :: Q(0:3)
        integer,  intent(in) :: k
        real(rb)             :: BkQ(4)
        select case (k)
          case (1); BkQ = [-Q(1),  Q(0),  Q(3), -Q(2)]
          case (2); BkQ = [-Q(2), -Q(3),  Q(0),  Q(1)]
          case (3); BkQ = [-Q(3),  Q(2), -Q(1),  Q(0)]
        end select
     end function Permutation
  end subroutine uniaxial_rotation
  !-------------------------------------------------------------------------------------------------
  function quaternion( A ) result( Q )
    real(rb), intent(in) :: A(3,3)
    real(rb)             :: Q(4)
    integer, parameter :: B(4,4) = reshape([1,1,1,1 ,1,1,-1,-1, 1,-1,1,-1, 1,-1,-1,1],[4,4])
    integer :: imax
    real(rb) :: Q2(4)
    Q2 = 0.25_rb*matmul(real(B,rb),[one,A(1,1),A(2,2),A(3,3)])
    imax = maxloc(Q2,1)
    Q(imax) = sqrt(Q2(imax))
    select case(imax)
      case (1)
        Q(2) = 0.25_rb*(A(2,3) - A(3,2))/Q(1)
        Q(3) = 0.25_rb*(A(3,1) - A(1,3))/Q(1)
        Q(4) = 0.25_rb*(A(1,2) - A(2,1))/Q(1)
      case (2)
        Q(1) = 0.25_rb*(A(2,3) - A(3,2))/Q(2)
        Q(3) = 0.25_rb*(A(1,2) + A(2,1))/Q(2)
        Q(4) = 0.25_rb*(A(1,3) + A(3,1))/Q(2)
      case (3)
        Q(1) = 0.25_rb*(A(3,1) - A(1,3))/Q(3)
        Q(2) = 0.25_rb*(A(1,2) + A(2,1))/Q(3)
        Q(4) = 0.25_rb*(A(2,3) + A(3,2))/Q(3)
      case (4)
        Q(1) = 0.25_rb*(A(1,2) - A(2,1))/Q(4)
        Q(2) = 0.25_rb*(A(1,3) + A(3,1))/Q(4)
        Q(3) = 0.25_rb*(A(2,3) + A(3,2))/Q(4)
    end select
    Q = Q/sqrt(sum(Q**2))
  end function quaternion
  !-------------------------------------------------------------------------------------------------
!  function rotation_matrix( Q ) result( A )
!    real(rb), intent(in) :: Q(4)
!    real(rb)             :: A(3,3)
!    real(rb) :: Q2(4), B12, B13, B14, B23, B24, B34
!    Q2 = Q*Q
!    B12 = two*Q(1)*Q(2)
!    B13 = two*Q(1)*Q(3)
!    B14 = two*Q(1)*Q(4)
!    B23 = two*Q(2)*Q(3)
!    B24 = two*Q(2)*Q(4)
!    B34 = two*Q(3)*Q(4)
!    A(1,:) = [ Q2(1)+Q2(2)-Q2(3)-Q2(4), B23+B14,                 B24-B13                 ]
!    A(2,:) = [ B23-B14,                 Q2(1)-Q2(2)+Q2(3)-Q2(4), B34+B12                 ]
!    A(3,:) = [ B24+B13,                 B34-B12,                 Q2(1)-Q2(2)-Q2(3)+Q2(4) ]
!  end function rotation_matrix
  !-------------------------------------------------------------------------------------------------
  function inverse_rotation_matrix( Q ) result( A )
    real(rb), intent(in) :: Q(4)
    real(rb)             :: A(3,3)
    real(rb) :: Q2(4), B12, B13, B14, B23, B24, B34
    Q2 = Q*Q
    B12 = two*Q(1)*Q(2)
    B13 = two*Q(1)*Q(3)
    B14 = two*Q(1)*Q(4)
    B23 = two*Q(2)*Q(3)
    B24 = two*Q(2)*Q(4)
    B34 = two*Q(3)*Q(4)
    A(:,1) = [ Q2(1)+Q2(2)-Q2(3)-Q2(4), B23+B14,                 B24-B13                 ]
    A(:,2) = [ B23-B14,                 Q2(1)-Q2(2)+Q2(3)-Q2(4), B34+B12                 ]
    A(:,3) = [ B24+B13,                 B34-B12,                 Q2(1)-Q2(2)-Q2(3)+Q2(4) ]
  end function inverse_rotation_matrix
  !-------------------------------------------------------------------------------------------------
  pure function matrix_B( Q ) result( B )
    real(rb), intent(in) :: Q(0:3)
    real(rb)             :: B(4,3)
    B = reshape( [-Q(1),  Q(0),  Q(3), -Q(2),  &
                  -Q(2), -Q(3),  Q(0),  Q(1),  &
                  -Q(3),  Q(2), -Q(1),  Q(0)], [4,3] )
  end function matrix_B
  !-------------------------------------------------------------------------------------------------
  pure function matrix_C( Q ) result( C )
    real(rb), intent(in) :: Q(0:3)
    real(rb)             :: C(4,3)
    C = reshape( [-Q(1),  Q(0), -Q(3),  Q(2),  &
                  -Q(2),  Q(3),  Q(0), -Q(1),  &
                  -Q(3), -Q(2),  Q(1),  Q(0)], [4,3] )
  end function matrix_C
  !-------------------------------------------------------------------------------------------------
  subroutine add_inertia( inertia, mass, delta )
    real(rb), intent(inout) :: inertia(3,3)
    real(rb), intent(in)    :: mass, delta(3)
    inertia(1,1) = inertia(1,1) + mass*(delta(2)**2 + delta(3)**2)
    inertia(2,2) = inertia(2,2) + mass*(delta(1)**2 + delta(3)**2)
    inertia(3,3) = inertia(3,3) + mass*(delta(1)**2 + delta(2)**2)
    inertia(1,2) = inertia(1,2) - mass*delta(1)*delta(2)
    inertia(1,3) = inertia(1,3) - mass*delta(1)*delta(3)
    inertia(2,3) = inertia(2,3) - mass*delta(2)*delta(3)
    inertia(2,1) = inertia(1,2)
    inertia(3,1) = inertia(1,3)
    inertia(3,2) = inertia(2,3)
  end subroutine add_inertia
  !-------------------------------------------------------------------------------------------------
end module mRigidBody
