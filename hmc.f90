program lj_hmc

use mGlobal
use mConfig
use mRandom
use EmDee

implicit none

real(rb), parameter :: mvv2e = 2390.057364_rb         ! Da*A²/fs² to kcal/mol
real(rb), parameter :: kB = 8.31451E-7_rb             ! Boltzmann constant in Da*A²/(fs²*K)
real(rb), parameter :: Pconv = 1.6388244954E+8_rb     ! Da/(A*fs²) to atm
real(rb), parameter :: kCoul = 0.13893545755135628_rb ! Coulomb constant in Da*A³/(fs*e)²

! Simulation specifications:
character(sl) :: Base
integer       :: i, j, N, seed, MDsteps, Nconf, thermo, Nequil, Nprod
real(rb)      :: T, Rc, dt, skin

! System properties:
integer  :: dof
real(rb) :: Volume

! Other variables:
integer  :: step
real(rb) :: half_dt, KE_sp, kT
real(rb) :: Acc_Press = zero
character(256) :: filename

integer :: threads
integer, allocatable, target :: indexes(:)
type(tEmDee), target :: md
type(EmDee_Model), pointer :: model(:)
type(c_ptr) :: system
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
call Config % Read( trim(Base)//".lmp" )
call Setup_Simulation

md = EmDee_system( threads, Rc, skin, Config%natoms, c_loc(Config%Type), c_loc(Config%mass), seed )
system = c_loc(md)

allocate( model(Config%ntypes) )
do i = 1, Config%ntypes
  if (abs(Config%epsilon(i)) < epsilon(1.0_rb)) then
    model(i) = EmDee_pair_none()
  else
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
call EmDee_random_momenta( system, kT, 1 )

call Config % Save_XYZ( trim(Base)//".xyz" )

call EmDee_compute( system )
call writeln( "Step Temp KinEng KinEng_t KinEng_r PotEng TotEng Press" )
step = 0
call writeln( properties() )
do step = 1, NEquil
  call Monte_Carlo_Step
!  call Verlet_Step
  if (mod(step,thermo) == 0) call writeln( properties() )
end do
call Report_Times()

Acc_Press = zero
call EmDee_compute( system )
call writeln( )
call writeln( "Memory usage" )
call writeln( "Step Temp KinEng KinEng_t KinEng_r PotEng TotEng Press" )
step = NEquil
call writeln( properties() )
do step = NEquil+1, NEquil+NProd
  call Monte_Carlo_Step
!  call Verlet_Step
  if (mod(step,Nconf)==0) then
    call EmDee_download( system, c_null_ptr, c_loc(Config%R), c_null_ptr, c_null_ptr )
    call Config % Save_XYZ( trim(Base)//".xyz", append = .true. )
  end if
  if (mod(step,thermo) == 0) call writeln( properties() )
end do
call Report_Times()

call Config % Write( trim(Base)//"_out.lmp", velocities = .true. )
call stop_log

contains
  !-------------------------------------------------------------------------------------------------
  subroutine Report_Times
    call writeln( "Loop time of", real2str(md%totalTime), "s." )
    call writeln()
    call writeln( "Pair time  =", real2str(md%pairTime), "s." )
    call writeln( "Other time =", real2str(md%totalTime-md%pairTime), "s." )
    call writeln()
    call writeln( "Neighbor list builds =", int2str(md%builds) )
  end subroutine Report_Times
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
                                 Pconv*((md%nbodies-1)*kB*Temp + md%Virial)/Volume]))
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
    dof = 6*maxval(Config%Mol) - 3
    KE_sp = half*dof*kB*T
    Volume = Lx*Ly*Lz
    kT = kB*T
    call random % setup( seed )
  end subroutine Setup_Simulation
  !-------------------------------------------------------------------------------------------------
  subroutine Verlet_Step
    call EmDee_boost( system, 1.0_rb, 0.0_rb, half_dt )
    call EmDee_move( system, 1.0_rb, 0.0_rb, dt )
    call EmDee_compute( system )
    call EmDee_boost( system, 1.0_rb, 0.0_rb, half_dt )
  end subroutine Verlet_Step
  !-------------------------------------------------------------------------------------------------
  subroutine Monte_Carlo_Step
    integer :: step
    real(rb) :: DeltaE
    real(rb), target :: Rsave(3,N)

    call EmDee_download( system, c_null_ptr, c_loc(Rsave), c_null_ptr, c_null_ptr )
    call EmDee_random_momenta( system, kB*T, 0 )
    DeltaE = md%Potential + md%Kinetic
    do step = 1, MDsteps
      call Verlet_Step
    end do
    DeltaE = md%Potential + md%Kinetic - DeltaE
    if (DeltaE >= zero) then
      if (random % uniform() > exp(-DeltaE/kT)) then
        call EmDee_upload( system, c_null_ptr, c_loc(Rsave), c_null_ptr )
        call EmDee_compute( system )
      end if
    end if
  end subroutine Monte_Carlo_Step
  !-------------------------------------------------------------------------------------------------
end program lj_hmc
