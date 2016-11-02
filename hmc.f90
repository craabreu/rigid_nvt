program lj_hmc

use mGlobal
use mConfig
use mRandom
use iso_c_binding

implicit none

real(rb), parameter :: mvv2e = 2390.057364_rb         ! Da*A²/fs² to kcal/mol
real(rb), parameter :: kB = 8.31451E-7_rb             ! Boltzmann constant in Da*A²/(fs²*K)
real(rb), parameter :: Pconv = 1.6388244954E+8_rb     ! Da/(A*fs²) to atm
real(rb), parameter :: kCoul = 0.13893545755135628_rb ! Coulomb constant in Da*A³/(fs*e)²

! Simulation specifications:
character(sl) :: Base
integer       :: i, j, N, nbodies, seed, MDsteps, Nconf, thermo, Nequil, Nprod
real(rb)      :: T, Rc, dt, skin

! System properties:
integer  :: dof
real(rb) :: Volume

! Other variables:
integer  :: step, nReject
real(rb) :: half_dt, KE_sp, kT
character(256) :: filename, configFile

#include "emdee.f03"

integer :: threads
integer, allocatable, target :: indexes(:)
type(tEmDee) :: system
type(c_ptr), allocatable :: model(:)
type(kiss) :: random

! Executable code:
#ifdef coul
  call writeln( "md/lj/coul ("//__DATE__//")" )
#else
  call writeln( "md/lj ("//__DATE__//")" )
#endif

call Get_Command_Line_Args( threads, filename )
call Read_Specifications( filename )

call init_log( trim(Base)//".log" )
call Config % Read( configFile )
call Setup_Simulation

system = EmDee_system( threads, 1, Rc, skin, Config%natoms, c_loc(Config%Type), c_loc(Config%mass) )

allocate( model(Config%ntypes) )
do i = 1, Config%ntypes
  if (abs(Config%epsilon(i)) < 1.e-8_rb) then
    model(i) = EmDee_model_none()
  else
    model(i) = EmDee_pair_lj( Config%epsilon(i)/mvv2e, Config%sigma(i) )
  end if
  call EmDee_set_pair_type( system, i, i, model(i) )
end do
nbodies = maxval(Config%Mol)
do i = 1, nbodies
  indexes = pack( [(j,j=1,N)], Config%Mol == i )
  call EmDee_add_rigid_body( system, size(indexes), c_loc(indexes) )
end do
#ifdef coul
  call EmDee_set_charges( system, c_loc(Config%Charge) )
#endif
call EmDee_upload( system, c_loc(Config%Lx), c_loc(Config%R), c_loc(Config%P), c_null_ptr )
!call EmDee_random_momenta( system, kT, 1, seed )
call Config % Save_XYZ( trim(Base)//".xyz" )

call writeln( "Step Temp KinEng KinEng_t KinEng_r PotEng TotEng Press Acceptance" )
step = 0
call writeln( properties(1) )
stop
nReject = 0
do step = 1, NEquil
  call Monte_Carlo_Step
  if (mod(step,thermo) == 0) call writeln( properties(step) )
end do
call Report( NEquil )
stop
call writeln( )
call writeln( "Memory usage" )
call writeln( "Step Temp KinEng KinEng_t KinEng_r PotEng TotEng Press Acceptance" )
step = NEquil
call writeln( properties( 1 ) )
nReject = 0
do step = NEquil+1, NEquil+NProd
  call Monte_Carlo_Step
  if (mod(step,Nconf)==0) then
    call EmDee_download( system, c_null_ptr, c_loc(Config%R), c_null_ptr, c_null_ptr )
    call Config % Save_XYZ( trim(Base)//".xyz", append = .true. )
  end if
  if (mod(step,thermo) == 0) call writeln( properties(step-Nequil) )
end do
call Report( NProd )

call EmDee_download( system, c_null_ptr, c_loc(Config%R), c_loc(Config%P), c_null_ptr )
call Config % Write( trim(Base)//"_out.lmp", velocities = .true. )
call stop_log

contains
  !-------------------------------------------------------------------------------------------------
  subroutine Report( Ntotal )
    integer, intent(in) :: Ntotal
    call writeln( "Loop time of", real2str(system%totalTime), "s." )
    call writeln()
    call writeln( "Pair time  =", real2str(system%pairTime), "s." )
    call writeln( "Other time =", real2str(system%totalTime-system%pairTime), "s." )
    call writeln()
    call writeln( "Neighbor list builds =", int2str(system%builds) )
    call writeln( )
    call writeln( "Acceptance ratio = ", real2str(one-real(nReject,rb)/Ntotal) )
  end subroutine Report
  !-------------------------------------------------------------------------------------------------
  character(sl) function properties( Ntotal )
    integer, intent(in) :: Ntotal
    real(rb) :: Temp
    Temp = (system%Kinetic/KE_sp)*T
    properties = trim(adjustl(int2str(step))) // " " // &
                 join(real2str([ Temp, &
                                 mvv2e*system%Kinetic, &
                                 mvv2e*(system%Kinetic - system%Rotational), &
                                 mvv2e*system%Rotational, &
                                 mvv2e*system%Potential, &
                                 mvv2e*(system%Potential + system%Kinetic), &
                                 Pconv*(nbodies*kB*Temp + system%Virial)/Volume, &
                                 one-real(nReject,rb)/Ntotal]))
  end function properties
  !-------------------------------------------------------------------------------------------------
  subroutine Get_Command_Line_Args( threads, filename )
    integer,        intent(out) :: threads
    character(256), intent(out) :: filename
    integer :: argcount
    character(256) :: line
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
    filename = line
  end subroutine Get_Command_Line_Args
  !-------------------------------------------------------------------------------------------------
  subroutine Read_Specifications( file )
    character(*), intent(in) :: file
    integer :: inp
    open( newunit=inp, file = file, status = "old" )
    read(inp,*); read(inp,*) Base
    read(inp,*); read(inp,*) configFile
    read(inp,*); read(inp,*) T
    read(inp,*); read(inp,*) Rc
    read(inp,*); read(inp,*) seed
    read(inp,*); read(inp,*) dt
    read(inp,*); read(inp,*) MDsteps
    read(inp,*); read(inp,*) skin
    read(inp,*); read(inp,*) Nconf
    read(inp,*); read(inp,*) thermo
    read(inp,*); read(inp,*) Nequil, Nprod
    close(inp)
    call writeln()
    call writeln( "Base for file names: ", Base )
    call writeln( "Name of configuration file:", configFile )
    call writeln( "Temperature: ", real2str(T), "K" )
    call writeln( "Cutoff distance: ", real2str(Rc), "Å" )
    call writeln( "Seed for random numbers: ", int2str(seed) )
    call writeln( "Time step size: ", real2str(dt), "fs" )
    call writeln( "Number of MD steps per MC step: ", int2str(MDsteps) )
    call writeln( "Skin size for neighbor lists: ", real2str(skin), "Å" )
    call writeln( "Interval for saving configurations: ", int2str(Nconf) )
    call writeln( "Interval for printing properties: ", int2str(thermo) )
    call writeln( "Number of equilibration steps: ", int2str(Nequil) )
    call writeln( "Number of production steps: ", int2str(Nprod) )
    call writeln()
  end subroutine Read_Specifications
  !-------------------------------------------------------------------------------------------------
  subroutine Setup_Simulation
    real(rb) :: Lx, Ly, Lz
    Config%Charge = sqrt(kCoul)*Config%Charge
    N = Config % natoms
    Lx = Config % Lx
    Ly = Config % Ly
    Lz = Config % Lz
    if (Rc+skin >= half*min(Lx,Ly,Lz)) call error( "minimum image convention failed!" )
    half_dt = half*dt
    dof = 6*maxval(Config%Mol)
    KE_sp = half*dof*kB*T
    Volume = Lx*Ly*Lz
    kT = kB*T
    call random % setup( seed )
  end subroutine Setup_Simulation
  !-------------------------------------------------------------------------------------------------
  subroutine Verlet_Step
    call EmDee_boost( system, 1.0_rb, 0.0_rb, half_dt, 1, 1 )
    call EmDee_move( system, 1.0_rb, 0.0_rb, dt )
    call EmDee_boost( system, 1.0_rb, 0.0_rb, half_dt, 1, 1 )
  end subroutine Verlet_Step
  !-------------------------------------------------------------------------------------------------
  subroutine Monte_Carlo_Step
    integer :: step
    real(rb) :: DeltaE, Potential, Virial
    real(rb), target :: Rsave(3,N), Fsave(3,N)
    DeltaE = system%Potential + system%Kinetic
    Potential = system%Potential
    Virial = system%Virial
    call EmDee_download( system, c_null_ptr, c_loc(Rsave), c_null_ptr, c_loc(Fsave) )
    do step = 1, MDsteps
      call Verlet_Step
    end do
    DeltaE = system%Potential + system%Kinetic - DeltaE
    if (DeltaE >= zero) then
      if (random % uniform() > exp(-DeltaE/kT)) then
        call EmDee_upload( system, c_null_ptr, c_loc(Rsave), c_null_ptr, c_loc(Fsave) )
        system%Potential = Potential
        system%Virial = Virial
        nReject = nReject + 1
      end if
    end if
    call EmDee_random_momenta( system, kT, 0, seed )
  end subroutine Monte_Carlo_Step
  !-------------------------------------------------------------------------------------------------
end program lj_hmc
