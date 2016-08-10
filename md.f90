program lj_md_nve

use mGlobal
use mConfig
use mThermostat
use EmDee

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

! System properties:
real(rb) :: Lx, Ly, Lz, Volume
integer :: dof

! Other variables:
integer  :: step, Nconf, Nprop, Nequil, Nprod
real(rb) :: InvLx, InvLy, InvLz
real(rb) :: Rc2, half_dt, KE_sp
real(rb) :: Acc_Press = zero
type(nhc_rev) :: NHC

integer :: argcount, threads
character(256) :: line
integer, allocatable, target :: indexes(:)
type(tEmDee), target :: md
type(c_ptr) :: system
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

call Read_Specifications( line )

call init_log( trim(Base)//".log" )
call Config % Read( trim(Base)//".lmp" )
call Setup_Simulation

md = EmDee_system( threads, Rc, skin, Config%natoms, c_loc(Config%Type), c_loc(Config%mass), seed )
system = c_loc(md)
md%forces = c_loc(Config%F)
md%coords = c_loc(Config%R)
md%momenta = c_loc(Config%P)

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
  call EmDee_add_rigid_body( system, c_loc(indexes), size(indexes) )
end do
#ifdef coul
  call EmDee_set_charges( system, c_loc(Config%Charge) )
#endif
call EmDee_upload( system, c_loc(Config%Lx), c_loc(Config%R), c_null_ptr )
call EmDee_random_momenta( system, kB*T, 1 )

call Config % Save_XYZ( trim(Base)//".xyz" )

call NHC % Setup( Nchain, kB*T, tdamp, dof, dt, loop = Nloop, nsy = Nsy )

call EmDee_compute( system )
call writeln( "Step Temp KinEng KinEng_t KinEng_r PotEng TotEng Press H_nhc" )
step = 0
call writeln( properties() )
do step = 1, NEquil
  call Verlet_Step
  if (mod(step,Nprop) == 0) call writeln( properties() )
end do

Acc_Press = zero
call EmDee_compute( system )
call writeln( )
call writeln( "Memory usage" )
call writeln( "Step Temp KinEng KinEng_t KinEng_r PotEng TotEng Press H_nhc" )
step = 0
call writeln( properties() )
do step = NEquil+1, NEquil+NProd
  call Verlet_Step
  if (mod(step,Nconf)==0) call Config % Save_XYZ( trim(Base)//".xyz", append = .true. )
  if (mod(step,Nprop) == 0) call writeln( properties() )
end do

call writeln( "Loop time of", real2str(md%totalTime), "s." )
call writeln()
call writeln( "Pair time  =", real2str(md%pairTime), "s." )
call writeln( "Other time =", real2str(md%totalTime-md%pairTime), "s." )
call writeln()
call writeln( "Neighbor list builds =", int2str(md%builds) )

call Config % Write( trim(Base)//"_out.lmp", velocities = .true. )
call stop_log

contains
  !-------------------------------------------------------------------------------------------------
  character(sl) function properties()
    real(rb) :: Temp
    Temp = (md%Kinetic/KE_sp)*T
    Acc_Press = Acc_Press + Pconv*(md%nbodies*kB*Temp + md%Virial)/Volume

    properties = trim(adjustl(int2str(step))) // " " // &
                 join(real2str([ Temp, &
                                 mvv2e*md%Kinetic, &
                                 mvv2e*(md%Kinetic - md%Rotational), &
                                 mvv2e*md%Rotational, &
                                 mvv2e*md%Potential, &
                                 mvv2e*(md%Potential + md%Kinetic), &
                                 Pconv*((md%nbodies-1)*kB*Temp + md%Virial)/Volume, &
                                 mvv2e*(md%Potential + md%Kinetic + NHC % energy()) ]))
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
    dof = 6*maxval(Config%Mol) - 3
    KE_sp = half*dof*kB*T
    Volume = Lx*Ly*Lz
    Nconf = nint(tconf/dt)
    Nprop = nint(thermo/dt)
    Nequil = nint(tequil/dt)
    Nprod = nint(tprod/dt)
  end subroutine Setup_Simulation
  !-------------------------------------------------------------------------------------------------
  subroutine Verlet_Step
    call EmDee_boost( system, half_dt )
    call EmDee_move( system, dt )
    call EmDee_compute( system )
    call EmDee_boost( system, half_dt )
  end subroutine Verlet_Step
  !-------------------------------------------------------------------------------------------------
end program lj_md_nve
