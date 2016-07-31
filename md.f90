program lj_md_nve

use mGlobal
use mRandom
use mConfig
use mThermostat
use mRigidBody
use EmDee
use omp_lib

implicit none

real(rb), parameter :: mvv2e = 2390.057364_rb         ! Da*A²/fs² to kcal/mol
real(rb), parameter :: kB = 8.31451E-7_rb             ! Boltzmann constant in Da*A²/(fs²*K)
real(rb), parameter :: Pconv = 1.6388244954E+8_rb     ! Da/(A*fs²) to atm
real(rb), parameter :: kCoul = 0.13893545755135628_rb ! Coulomb constant in Da*A³/(fs*e)²

! Simulation specifications:
character(sl) :: Base
integer :: i, j, N, seed, Nchain, Nloop, Nsy
real(rb) :: T, Rc, dt, skin, tdamp, tconf, thermo, tequil, tprod
character(3) :: ensemble
logical :: nvt = .false.

! System properties:
real(rb) :: Lx, Ly, Lz, Volume, PE, KE, Virial
integer :: NB, dof
integer, allocatable :: first(:), last(:)
type(tRigidBody), allocatable :: body(:)

! Other variables:
integer  :: step, Nconf, Nprop, Nequil, Nprod
logical  :: Compute
real(rb) :: InvLx, InvLy, InvLz
real(rb) :: Rc2, half_dt, KE_sp
real(rb) :: Acc_Press = zero
real(rb) :: total_time
type(shr3) :: random
type(nhc_rev) :: NHC

integer :: argcount, threads
character(256) :: line
integer, allocatable, target :: indexes(:)
type(tEmDee), target :: md
type(c_ptr) :: system, forces, coords
type(EmDee_Model), pointer :: model(:)

! Executable code:
#ifdef coul
  call writeln( "md/lj/coul ("//__DATE__//")" )
#else
  call writeln( "md/lj ("//__DATE__//")" )
#endif

argcount = command_argument_count()
if (argcount == 1) then
  threads = 1
  call get_command_argument( 1, line )
else if (argcount == 2) then
  call get_command_argument( 1, line )
  read(line,*) threads
  call get_command_argument( 2, line )
else
  write(0,'("Usage: md_lj/md_lj_coul [number-of-threads] input-file")')
  stop
end if

total_time = omp_get_wtime()
call Read_Specifications( line )

call init_log( trim(Base)//".log" )
call Config % Read( trim(Base)//".lmp" )
call Setup_Simulation

md = EmDee_system( threads, Rc, skin, Config%natoms, c_loc(Config%Type), c_loc(Config%mass) )
system = c_loc(md)
forces = c_loc(Config%F)
coords = c_loc(Config%R)

allocate( model(Config%ntypes) )
do i = 1, Config%ntypes
  if (abs(Config%epsilon(i)) < epsilon(1.0_rb)) then
    model(i) = EmDee_pair_none()
  else
!    model(i) = EmDee_pair_lj( Config%epsilon(i)/mvv2e, Config%sigma(i) )
    model(i) = EmDee_pair_lj_sf( Config%epsilon(i)/mvv2e, Config%sigma(i), Rc )
  end if
  call EmDee_set_pair_type( system, i, i, c_loc(model(i)) )
end do
do i = 1, maxval(Config%Mol)
  indexes = pack( [(j,j=1,N)], Config%Mol == i )
  call EmDee_add_rigid_body( system, size(indexes), c_loc(indexes), c_loc(Config%R), Config%Lx )
end do
#ifdef coul
  call EmDee_set_charges( system, c_loc(Config%Charge) )
#endif

call Identify_Rigid_Bodies( NB )
call Assign_Momenta
call writeln( "Ensemble: ", ensemble )
call Config % Save_XYZ( trim(Base)//".xyz" )

call NHC % Setup( Nchain, kB*T, tdamp, dof, dt, loop = Nloop, nsy = Nsy )

Compute = .true.
call Compute_Forces
KE = sum(body%kinetic_energy())
call writeln( "Step Temp KinEng KinEng_t KinEng_r PotEng TotEng Press H_nhc" )
step = 0
call writeln( properties() )
do step = 1, NEquil
  Compute = mod(step,Nprop) == 0
  call Verlet_Step
  if (Compute) call writeln( properties() )
end do

Acc_Press = zero
Compute = .true.
call Compute_Forces
call writeln( )
call writeln( "Memory usage" )
call writeln( "Step Temp KinEng KinEng_t KinEng_r PotEng TotEng Press H_nhc" )
step = 0
call writeln( properties() )
do step = NEquil+1, NEquil+NProd
  Compute = mod(step,Nprop) == 0
  call Verlet_Step
  if (mod(step,Nconf)==0) call Config % Save_XYZ( trim(Base)//".xyz", append = .true. )
  if (Compute) call writeln( properties() )
end do

total_time = omp_get_wtime() - total_time
call writeln( "Loop time of", real2str(total_time), "s." )
call writeln()
call writeln( "Pair time  =", real2str(md%time), "s." )
call writeln( "Other time =", real2str(total_time-md%time), "s." )
call writeln()
call writeln( "Neighbor list builds =", int2str(md%builds) )

call Get_Momenta
call Config % Write( trim(Base)//"_out.lmp", velocities = .true. )
call stop_log

contains
  !-------------------------------------------------------------------------------------------------
  character(sl) function properties()
    real(rb) :: Temp
    Temp = (KE/KE_sp)*T
    Acc_Press = Acc_Press + Pconv*(NB*kB*Temp + Virial)/Volume

    properties = trim(adjustl(int2str(step))) // " " // &
                 join(real2str([ Temp, &
                                 mvv2e*KE, &
                                 mvv2e*sum(body%translational_energy()), &
                                 mvv2e*sum(body%rotational_energy()), &
                                 mvv2e*PE, &
                                 mvv2e*(PE + KE), &
                                 Pconv*((NB-1)*kB*Temp + Virial)/Volume, &
                                 mvv2e*(PE + KE + NHC % energy()) ]))
  end function properties
  !-------------------------------------------------------------------------------------------------
  subroutine Read_Specifications( file )
    character(*), intent(in) :: file
    integer :: inp
    real(rb) :: alpha
    character(sl) :: damp
    open( newunit=inp, file = file, status = "old" )
    read(inp,*); read(inp,*) Base
    read(inp,*); read(inp,*) T
    read(inp,*); read(inp,*) Rc
    read(inp,*); read(inp,*) alpha, damp
    read(inp,*); read(inp,*) seed
    read(inp,*); read(inp,*) dt
    read(inp,*); read(inp,*) Nchain, Nloop, Nsy
    read(inp,*); read(inp,*) tdamp
    read(inp,*); read(inp,*) skin
    read(inp,*); read(inp,*) tconf
    read(inp,*); read(inp,*) thermo
    read(inp,*); read(inp,*) tequil, tprod
    read(inp,*); read(inp,*) ensemble
    close(inp)
    call writeln()
    call writeln( "Base for file names: ", Base )
    call writeln( "Temperature: ", real2str(T), "K" )
    call writeln( "Cutoff distance: ", real2str(Rc), "Å" )
    call writeln( "Coulombic damping model: ", damp )
    if (damp /= "none") &
      call writeln( "Coulombic damping constant: ", real2str(alpha), "Å^(-1)" )
    call writeln( "Seed for random numbers: ", int2str(seed) )
    call writeln( "Time step size: ", real2str(dt), "fs" )
    call writeln( "Thermostat chain size: ", int2str(Nchain) )
    call writeln( "Thermostat integration loops: ", int2str(Nloop) )
    call writeln( "Thermostat Suzuki-Yoshida weights: ", int2str(Nsy) )
    call writeln( "Thermostat time constant: ", real2str(tdamp), "fs" )
    call writeln( "Skin size for neighbor lists: ", real2str(skin), "Å" )
    call writeln( "Time interval for saving configurations: ", real2str(tconf), "fs" )
    call writeln( "Time interval for printing properties: ", real2str(thermo), "fs" )
    call writeln( "Total equilibration time: ", real2str(tequil), "fs" )
    call writeln( "Total production time: ", real2str(tprod), "fs" )
    call writeln( "Simulated ensemble: ", ensemble )
    call writeln()
  end subroutine Read_Specifications
  !-------------------------------------------------------------------------------------------------
  subroutine Setup_Simulation
    if ((ensemble /= "nvt").and.(ensemble /= "nve")) call error( "wrong ensemble definition" )
    Config%Charge = sqrt(kCoul)*Config%Charge
    N = Config % natoms
    Lx = Config % Lx
    Ly = Config % Ly
    Lz = Config % Lz
    if (Rc+skin >= half*min(Lx,Ly,Lz)) call error( "minimum image convention failed!" )
    InvLx = 1.0_rb/Lx
    InvLy = 1.0_rb/Ly
    InvLz = 1.0_rb/Lz
    Rc2 = Rc**2
    half_dt = half*dt
    NB = maxval(Config%Mol)
    nvt = ensemble == "nvt"
    dof = 6*NB - 3
    KE_sp = half*dof*kB*T
    Volume = Lx*Ly*Lz
    Nconf = nint(tconf/dt)
    Nprop = nint(thermo/dt)
    Nequil = nint(tequil/dt)
    Nprod = nint(tprod/dt)
  end subroutine Setup_Simulation
  !-------------------------------------------------------------------------------------------------
  subroutine Identify_Rigid_Bodies( NB )
    integer, intent(in) :: NB
    integer :: i, j, k
    real(rb) :: Rijx, Rijy, Rijz
    allocate( first(NB), last(NB), body(NB) )
    first = 0
    do i = 1, N
      j = Config%Mol(i)
      if (first(j) == 0) first(j) = i
      last(j) = i
    end do
    do k = 1, NB
      i = first(k)
      do j = i+1, last(k)
        Rijx = Config%Rx(i) - Config%Rx(j)
        Rijy = Config%Ry(i) - Config%Ry(j)
        Rijz = Config%Rz(i) - Config%Rz(j)
        Rijx = Rijx - Lx*nint(InvLx*Rijx)
        Rijy = Rijy - Ly*nint(InvLx*Rijy)
        Rijz = Rijz - Lz*nint(InvLx*Rijz)
        Config%Rx(j) = Config%Rx(i) - Rijx
        Config%Ry(j) = Config%Ry(i) - Rijy
        Config%Rz(j) = Config%Rz(i) - Rijz
      end do
      j = last(k)
      call body(k) % setup( Config%Mass(Config%Type(i:j)), Config%Rx(i:j), Config%Ry(i:j), Config%Rz(i:j) )
    end do
    if (maxval(body%maximum_distance()) + two*Rc > min(Lx,Ly,Lz)) then
      call warning( "minimum image convention holds for atoms, but not for molecules" )
      call warning( "pressure computation might be incorrect" )
    end if
  end subroutine Identify_Rigid_Bodies
  !-------------------------------------------------------------------------------------------------
  subroutine Assign_Momenta
    integer  :: i, p1, pN
    real(rb) :: kBT, P(3), Mass
    if (Config % velocity_input) then
      do i = 1, NB
        p1 = first(i)
        pN = last(i)
        call body(i)%set_momenta( Config%Px(p1:pN), Config%Py(p1:pN), Config%Pz(p1:pN) )
      end do
    else
      call random%setup( seed )
      kBT = kB*T
      do i = 1, NB
        call body(i) % assign_momenta( kBT, random, rotation = .false. )
      end do
    end if
    P = [sum(body%Pcm(1)),sum(body%Pcm(2)),sum(body%Pcm(3))]
    Mass = sum(body%Mass)
    do i = 1, NB
      body(i)%Pcm = body(i)%Pcm - (body(i)%Mass/Mass)*P
    end do
    KE = sum(body%kinetic_energy())
  end subroutine Assign_Momenta
  !-------------------------------------------------------------------------------------------------
  subroutine Get_Momenta
    integer :: i, p1, pN
    do i = 1, NB
      p1 = first(i)
      pN = last(i)
      call body(i)%get_momenta( Config%Px(p1:pN), Config%Py(p1:pN), Config%Pz(p1:pN) )
    end do
  end subroutine Get_Momenta
  !-------------------------------------------------------------------------------------------------
  subroutine Compute_Forces
    integer  :: i, k, p1, pN
    real(rb) :: W
    call EmDee_compute( system, forces, coords, Config%Lx )
    if (Compute) then
      PE = md%Energy
      W = zero
      k = 0
      do i = 1, NB
        do j = 1, body(i)%NP
          k = k + 1
          W = W + sum(body(i)%delta(:,j)*Config%F(:,k))
        end do
      end do
      Virial = md%Virial - third*W
    end if
    do i = 1, NB
      p1 = first(i)
      pN = last(i)
      call body(i) % force_and_torque( Config%Fx(p1:pN), Config%Fy(p1:pN), Config%Fz(p1:pN) )
    end do
  end subroutine Compute_Forces
  !-------------------------------------------------------------------------------------------------
  subroutine Verlet_Step
    integer :: i, p1, pN
    if (nvt) then
      call NHC % integrate( body )
      call NHC % integrate_momenta( body,half_dt )
    else   
     call body % boost( half_dt )
    end if
    call body % rotate( dt )
    call body % displace( dt )
    do i = 1, NB
      p1 = first(i)
      pN = last(i)
      call body(i)%get_positions( Config%Rx(p1:pN), Config%Ry(p1:pN), Config%Rz(p1:pN) )
    end do
    call Compute_Forces
    if (nvt) then
      call NHC % integrate_momenta( body,half_dt )
      call NHC % integrate( body )      
    else   
     call body % boost( half_dt )
    end if
    if (compute) KE = sum(body%kinetic_energy())
    call Get_Momenta
  end subroutine Verlet_Step
  !-------------------------------------------------------------------------------------------------
end program lj_md_nve
