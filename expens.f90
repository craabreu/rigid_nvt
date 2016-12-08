program expanded_ensemble

use mGlobal
use mConfig
use mRandom
use mThermostat
use mSamplingStats
use iso_c_binding

implicit none

#define isZero(X) (abs(X) < epsilon(1.0_8))

real(rb), parameter :: mvv2e = 2390.057364_rb         ! Da*A²/fs² to kcal/mol
real(rb), parameter :: kB = 8.31451E-7_rb             ! Boltzmann constant in Da*A²/(fs²*K)
real(rb), parameter :: Pconv = 1.6388244954E+8_rb     ! Da/(A*fs²) to atm
real(rb), parameter :: kCoul = 0.13893545755135628_rb ! Coulomb constant in Da*A³/(fs*e)²

! Simulation specifications:
character(sl) :: Base
integer       :: i, j, k, N, NB, seed, MDsteps, Nconf, thermo, Nequil, Nprod, rotationMode
real(rb)      :: T, Rc, dt, skin, beta

! System properties:
integer  :: dof
real(rb) :: Volume

! Thermostat variables:
integer :: M, ndamp, nloops
type(nhc_pscaling) :: thermostat

! Expanded ensemble variables:
logical :: iInSolute
integer :: Nlayers, layer
integer,  allocatable :: soluteType(:)
real(rb), allocatable :: lambda(:), eta(:)
real(rb), pointer :: energy(:)

! Other variables:
integer  :: step
real(rb) :: dt_2, dt_4, KE_sp, kT
real(rb), pointer :: charges(:)
character(256) :: filename, configFile

#include "emdee.f03"

integer :: threads
integer, allocatable, target :: indexes(:)
type(tEmDee) :: md
type(c_ptr) :: model
type(c_ptr), allocatable :: multimodel(:)
type(kiss) :: random
type(SamplingStats) :: stats

! Executable code:
call Get_Command_Line_Args( threads, filename )
call Read_Specifications( filename )

call Config % Read( configFile )
call Setup_Simulation

if (.not.allocated( Config%epsilon )) call error( "Pair coefficients were not found" )

md = EmDee_system( threads, Nlayers, Rc, skin, Config%natoms, c_loc(Config%Type), c_loc(Config%mass) )
md%Options%rotationMode = rotationMode
call c_f_pointer( md%layerEnergy, energy, [Nlayers] )

do i = 1, Config%ntypes
  if (.not.any(soluteType == i)) then
    if (isZero(Config%epsilon(i))) then
      model = EmDee_pair_coul_sf()
    else
      model = EmDee_pair_lj_sf_coul_sf( Config%epsilon(i)/mvv2e, Config%sigma(i) )
    end if
    call EmDee_set_pair_model( md, i, i, model )
  end if
end do

! VALID ONLY FOR A NON-POLAR SOLUTE:
allocate( multimodel(Nlayers) )
do i = 1, Config%ntypes-1
  iInSolute = any(soluteType == i)
  do j = i+1, Config%ntypes
    if (.not.( (iInSolute .eqv. any(soluteType == j)) .or. &
               isZero(Config%epsilon(i))   .or. &
               isZero(Config%epsilon(j))        )) then
      do k = 1, Nlayers
        multimodel(k) = EmDee_pair_softcore_cut( sqrt(Config%epsilon(i)*Config%epsilon(j))/mvv2e, &
                                         half*(Config%sigma(i) + Config%sigma(j)), lambda(k) )
      end do
      call EmDee_set_pair_multimodel( md, i, j, multimodel )
    end if
  end do
end do

call EmDee_set_charges( md, c_loc(charges) )
call EmDee_switch_model_layer( md, layer )

do i = 1, maxval(Config%Mol)
  indexes = pack( [(j,j=1,N)], Config%Mol == i )
  if (size(indexes) > 1) call EmDee_add_rigid_body( md, size(indexes), c_loc(indexes) )
end do

call EmDee_upload( md, c_loc(Config%Lx), c_loc(Config%R), c_null_ptr, c_null_ptr )
call EmDee_random_momenta( md, kT, 1, seed )
call Config % Save_XYZ( trim(Base)//".xyz" )

call writeln( "Step Temp KinEng KinEng_t KinEng_r PotEng TotEng Virial Press H_nhc node", &
              join([("E"//int2str(i),i=1,Nlayers)]) )
step = 0
call writeln( properties() )
do step = 1, NEquil
  call Pscaling_Step
  call Periodic_Tasks( step )
end do
call Report

call writeln( )
call writeln( "Memory usage" )
call writeln( "Step Temp KinEng KinEng_t KinEng_r PotEng TotEng Virial Press H_nhc node", &
              join([("E"//int2str(i),i=1,Nlayers)]) )
step = 0
call writeln( trim(properties())//" "//join([int2str(layer),real2str(mvv2e*energy)]) )

call stats % initialize( Nlayers )
do step = 1, NProd
  call Pscaling_Step
  call Periodic_Tasks( step )
end do
call Report

call Config % Write( trim(Base)//"_out.lmp", velocities = .true. )
call stats % save( trim(Base)//"_stats.csv", lambda, eta )
call stop_log

contains
  !-------------------------------------------------------------------------------------------------
  subroutine Periodic_Tasks( step )
    integer, intent(in) :: step

    logical  :: print_props
    integer  :: newLayer
    real(rb) :: u(Nlayers)

    print_props = mod(step,thermo) == 0
    if (print_props) then
      call writeln( trim(properties())//" "//join([int2str(layer),real2str(mvv2e*energy)]) )
    end if

    if (mod(step,MDsteps) == 0) then
      u = eta - beta*energy
      newLayer = random % label( exp(u - maxval(u)) )
      if (newLayer /= layer) then
        layer = newLayer
        call EmDee_switch_model_layer( md, layer )
      end if
      if (stats % initialized) call stats % sample( layer )
    end if

  end subroutine Periodic_Tasks
  !-------------------------------------------------------------------------------------------------
  subroutine Report( )
    call writeln( "Loop time of", real2str(md%totalTime), "s." )
    call writeln()
    call writeln( "Pair time  =", real2str(md%pairTime), "s." )
    call writeln( "Other time =", real2str(md%totalTime-md%pairTime), "s." )
    call writeln()
    call writeln( "Neighbor list builds =", int2str(md%builds) )
  end subroutine Report
  !-------------------------------------------------------------------------------------------------
  character(sl) function properties()
    real(rb) :: Temp, Etotal
    Temp = (md%Kinetic/KE_sp)*T
    Etotal = md%Potential + md%Kinetic
    properties = trim(adjustl(int2str(step))) // " " // &
                 join(real2str([ Temp, &
                                 mvv2e*md%Kinetic, &
                                 mvv2e*(md%Kinetic - md%Rotational), &
                                 mvv2e*md%Rotational, &
                                 mvv2e*md%Potential, &
                                 mvv2e*Etotal, &
                                 mvv2e*md%Virial, &
                                 Pconv*((NB-1)*kB*Temp + md%Virial)/Volume, &
                                 mvv2e*(Etotal + thermostat%energy()) ]))
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
    integer :: inp, n
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
    read(inp,*); read(inp,*) ndamp, M, nloops
    read(inp,*); read(inp,*) rotationMode
    read(inp,*); read(inp,*) n
    allocate( soluteType(n) )
    read(inp,*); read(inp,*) soluteType
    read(inp,*); read(inp,*) Nlayers
    allocate( lambda(Nlayers), eta(Nlayers) )
    read(inp,*); read(inp,*) lambda
    read(inp,*); read(inp,*) eta
    read(inp,*); read(inp,*) layer
    close(inp)
    call init_log( trim(Base)//".log" )
    call writeln()
    call writeln( "Base for file names:", Base )
    call writeln( "Name of configuration file:", configFile )
    call writeln( "Temperature:", real2str(T), "K" )
    call writeln( "Cutoff distance:", real2str(Rc), "Å" )
    call writeln( "Seed for random numbers:", int2str(seed) )
    call writeln( "Time step size:", real2str(dt), "fs" )
    call writeln( "Number of MD steps per MC step:", int2str(MDsteps) )
    call writeln( "Skin size for neighbor lists:", real2str(skin), "Å" )
    call writeln( "Interval for saving configurations:", int2str(Nconf) )
    call writeln( "Interval for printing properties:", int2str(thermo) )
    call writeln( "Number of equilibration steps:", int2str(Nequil) )
    call writeln( "Number of production steps:", int2str(Nprod) )
    call writeln( "Thermostat parameters:", int2str(ndamp), int2str(M), int2str(nloops) )
    if (rotationMode == 0) then
      call writeln( "Rotation mode: exact solution" )
    else
      call writeln( "Rotation mode: Miller with", int2str(rotationMode), "respa steps" )
    end if
    call writeln()
  end subroutine Read_Specifications
  !-------------------------------------------------------------------------------------------------
  subroutine Setup_Simulation
    real(rb) :: Lx, Ly, Lz
    N = Config % natoms
    allocate( charges(N) )
    charges = sqrt(kCoul)*Config%Charge
    Lx = Config % Lx
    Ly = Config % Ly
    Lz = Config % Lz
    if (Rc+skin >= half*min(Lx,Ly,Lz)) call error( "minimum image convention failed!" )
    dt_2 = half*dt
    dt_4 = 0.25_rb*dt
    NB = maxval(Config%Mol)
    dof = 6*NB - 3
    kT = kB*T
    beta = one/kT
    KE_sp = half*dof*kT
    Volume = Lx*Ly*Lz
    call random % setup( seed )
    call thermostat % setup( M, kT, ndamp*dt, 6*NB-3, nloops )
  end subroutine Setup_Simulation
  !-------------------------------------------------------------------------------------------------
  subroutine Pscaling_Step
    call thermostat % integrate( dt_2, two*md%Kinetic )
    call EmDee_boost( md, zero, thermostat%damping, dt_2 )
    call EmDee_boost( md, one, zero, dt_2 )
    call EmDee_move( md, one, zero, dt )
    call EmDee_boost( md, one, zero, dt_2 )
    call thermostat % integrate( dt_2, two*md%Kinetic )
    call EmDee_boost( md, zero, thermostat%damping, dt_2 )
  end subroutine Pscaling_Step
  !-------------------------------------------------------------------------------------------------
  subroutine Attempt_Lambda_Change
  end subroutine Attempt_Lambda_Change
  !-------------------------------------------------------------------------------------------------
end program expanded_ensemble
