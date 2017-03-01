#include "emdee.f03"

program lj_nvt

use mGlobal
use mConfig
use mRandom
use mThermostat
use iso_c_binding
use EmDee

implicit none

real(rb), parameter :: mvv2e = 2390.057364_rb         ! Da*A²/fs² to kcal/mol
real(rb), parameter :: kB = 8.31451E-7_rb             ! Boltzmann constant in Da*A²/(fs²*K)
real(rb), parameter :: Pconv = 1.6388244954E+8_rb     ! Da/(A*fs²) to atm
real(rb), parameter :: kCoul = 0.13893545755135628_rb ! Coulomb constant in Da*A³/(fs*e)²

! Simulation specifications:
character(sl) :: Base
integer       :: i, N, NB, seed, MDsteps, Nconf, thermo, Nequil, Nprod, rotationMode, Nrespa
real(rb)      :: T, Rc, dt, skin, respaRc

! System properties:
integer  :: dof
real(rb) :: Volume

! Thermostat variables:
integer :: method, M, ndamp, nloops, nts
logical :: single
class(nhc), pointer :: thermostat(:)

! Other variables:
integer  :: step
real(rb) :: dt_2, dt_4, KE_sp, kT
character(256) :: filename, configFile

integer :: threads
type(tEmDee) :: md
type(c_ptr), allocatable :: model(:)
type(kiss) :: random

#define trans_rot md%Options%translate = .true.;  md%Options%rotate = .true.
#define transOnly md%Options%translate = .true.;  md%Options%rotate = .false.
#define rotOnly   md%Options%translate = .false.; md%Options%rotate = .true.

! Executable code:
#ifdef coul
  call writeln( "md/lj/coul ("//__DATE__//")" )
#else
  call writeln( "md/lj ("//__DATE__//")" )
#endif

call Get_Command_Line_Args( threads, filename )
call Read_Specifications( filename )

call Config % Read( configFile )
call Setup_Simulation

md = EmDee_system( threads, 1, Rc, skin, Config%natoms, &
                   c_loc(Config%Type), c_loc(Config%mass), c_loc(Config%Mol) )
md%Options%rotationMode = rotationMode

allocate( model(Config%ntypes) )
do i = 1, Config%ntypes
  if (abs(Config%epsilon(i)) < epsilon(1.0_rb)) then
    model(i) = EmDee_pair_none()
  else
    model(i) = EmDee_shifted_force( EmDee_pair_lj_cut( Config%epsilon(i)/mvv2e, Config%sigma(i) ) )
  end if
  call EmDee_set_pair_model( md, i, i, model(i), kCoul )
end do
#ifdef coul
call EmDee_set_coul_model( md, EmDee_coul_sf() )
call EmDee_upload( md, "charges"//c_null_char, c_loc(Config%Charge) )
#endif

if (Nrespa > 0) call EmDee_set_respa( md, respaRc, Nrespa, 0 )

call EmDee_upload( md, "box"//c_null_char, c_loc(Config%Lx) )
call EmDee_upload( md, "coordinates"//c_null_char, c_loc(Config%R) )
call EmDee_random_momenta( md, kT, .true._1, seed )
call Config % Save_XYZ( trim(Base)//".xyz" )

call writeln( "Step Temp KinEng KinEng_t KinEng_r KinEng_r1 KinEng_r2 KinEng_r3 "// &
              "PotEng DispEng CoulEng TotEng Virial BodyVirial Press H_nhc" )
step = 0
call writeln( properties() )
do step = 1, NEquil
  select case (method)
    case (0); call Verlet_Step
    case (1); call Pscaling_Step
    case (2); call Boosting_Step
    case (3); call Kamberaj_Step
    case (4); call Hybrid_Step
    case (5); call New_Hybrid_Step
  end select
  if (mod(step,thermo) == 0) call writeln( properties() )
end do
call Report

call writeln( )
call writeln( "Memory usage" )
call writeln( "Step Temp KinEng KinEng_t KinEng_r KinEng_r1 KinEng_r2 KinEng_r3 "// &
              "PotEng DispEng CoulEng TotEng Virial BodyVirial Press H_nhc" )
step = NEquil
call writeln( properties() )
do step = NEquil+1, NEquil+NProd
  select case (method)
    case (0); call Verlet_Step
    case (1); call Pscaling_Step
    case (2); call Boosting_Step
    case (3); call Kamberaj_Step
    case (4); call Hybrid_Step
    case (5); call New_Hybrid_Step
  end select
  if (mod(step,Nconf)==0) then
    call EmDee_download( md, "coordinates"//c_null_char, c_loc(Config%R) )
    call Config % Save_XYZ( trim(Base)//".xyz", append = .true. )
  end if
  if (mod(step,thermo) == 0) call writeln( properties() )
end do
call Report

call Config % Write( trim(Base)//"_out.lmp", velocities = .true. )
call stop_log

contains
  !-------------------------------------------------------------------------------------------------
  subroutine Report
    call writeln( "Loop time of", real2str(md%totalTime), "s." )
    call writeln()
    call writeln( "Pair time  =", real2str(md%pairTime), "s." )
    call writeln( "Other time =", real2str(md%totalTime-md%pairTime), "s." )
    call writeln()
    call writeln( "Neighbor list builds =", int2str(md%builds) )
  end subroutine Report
  !-------------------------------------------------------------------------------------------------
  character(sl) function properties()
    real(rb) :: Temp
    real(rb) :: Etotal
    Temp = (md%Energy%Kinetic/KE_sp)*T
    Etotal = md%Energy%Potential + md%Energy%Kinetic
    properties = trim(adjustl(int2str(step))) // " " // &
                 join(real2str([ Temp, &
                                 mvv2e*md%Energy%Kinetic, &
                                 mvv2e*(md%Energy%Kinetic - md%Energy%Rotational), &
                                 mvv2e*md%Energy%Rotational, &
                                 mvv2e*md%Energy%RotPart, &
                                 mvv2e*md%Energy%Potential, &
                                 mvv2e*md%Energy%Dispersion, &
                                 mvv2e*md%Energy%Coulomb, &
                                 mvv2e*Etotal, &
                                 mvv2e*md%Virial, &
                                 mvv2e*md%BodyVirial, &
                                 Pconv*((NB-1)*kB*Temp + (1.0_rb/3.0_rb)*md%Virial)/Volume, &
                                 mvv2e*(Etotal + sum(thermostat%energy())) ]))
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
    read(inp,*); read(inp,*) method
    read(inp,*); read(inp,*) ndamp, M, nloops
    read(inp,*); read(inp,*) nts
    read(inp,*); read(inp,*) rotationMode
    read(inp,*); read(inp,*) Nrespa
    read(inp,*); read(inp,*) respaRc
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
    call writeln( "Thermostat method:", int2str(method) )
    call writeln( "Thermostat parameters:", int2str(ndamp), int2str(M), int2str(nloops) )
    call writeln( "Number of thermostat chains:", int2str(nts) )
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
    Lx = Config % Lx
    Ly = Config % Ly
    Lz = Config % Lz
    if (Rc+skin >= half*min(Lx,Ly,Lz)) call error( "minimum image convention failed!" )
    dt_2 = half*dt
    dt_4 = 0.25_rb*dt
    NB = maxval(Config%Mol)
    dof = 6*NB - 3
    kT = kB*T
    KE_sp = half*dof*kT
    Volume = Lx*Ly*Lz
    call random % setup( seed )
    if ((nts < 1).or.(nts > 2)) call error( "wrong translation/rotation thermostat scheme" )
    single = (nts == 1)
    select case (method)
      case (0,1,4); allocate( nhc_pscaling :: thermostat(nts) )
      case (2,5); allocate( nhc_boosting :: thermostat(nts) )
      case (3); allocate( nhc_kamberaj :: thermostat(nts) )
      case default; call error( "unknown thermostat method" )
    end select
    if (method == 5) then
      call thermostat(1) % setup( M, kT, ndamp*dt, (6/nts)*NB-3, 1 )
      if (nts == 2) call thermostat(2) % setup( M, kT, ndamp*dt, 3*NB, 1 )
    else
      call thermostat(1) % setup( M, kT, ndamp*dt, (6/nts)*NB-3, nloops )
      if (nts == 2) call thermostat(2) % setup( M, kT, ndamp*dt, 3*NB, nloops )
    end if
  end subroutine Setup_Simulation
  !-------------------------------------------------------------------------------------------------
  subroutine Verlet_Step
!    call EmDee_boost( md, one, zero, dt_2 )
!    call EmDee_displace( md, one, zero, dt )
!    call EmDee_boost( md, one, zero, dt_2 )
    call EmDee_advance( md, zero, zero, dt )
  end subroutine Verlet_Step
  !-------------------------------------------------------------------------------------------------
  subroutine Pscaling_Step
    if (single) then
      call thermostat(1) % integrate( dt_2, two*md%Energy%Kinetic )
      call EmDee_boost( md, zero, thermostat(1)%damping, dt_2 )
      call EmDee_boost( md, one, zero, dt_2 )
      call EmDee_displace( md, one, zero, dt )
      call EmDee_boost( md, one, zero, dt_2 )
      call thermostat(1) % integrate( dt_2, two*md%Energy%Kinetic )
      call EmDee_boost( md, zero, thermostat(1)%damping, dt_2 )
    else
      call thermostat(1) % integrate( dt_2, two*(md%Energy%Kinetic - md%Energy%Rotational) )
      call thermostat(2) % integrate( dt_2, two*md%Energy%Rotational )
      transOnly
      call EmDee_boost( md, zero, thermostat(1)%damping, dt_2 )
      rotOnly
      call EmDee_boost( md, zero, thermostat(2)%damping, dt_2 )
      trans_rot
      call EmDee_boost( md, one, zero, dt_2 )
      call EmDee_displace( md, one, zero, dt )
      call EmDee_boost( md, one, zero, dt_2 )
      call thermostat(1) % integrate( dt_2, two*(md%Energy%Kinetic - md%Energy%Rotational) )
      call thermostat(2) % integrate( dt_2, two*md%Energy%Rotational )
      transOnly
      call EmDee_boost( md, zero, thermostat(1)%damping, dt_2 )
      rotOnly
      call EmDee_boost( md, zero, thermostat(2)%damping, dt_2 )
      trans_rot
    end if
  end subroutine Pscaling_Step
  !-------------------------------------------------------------------------------------------------
  subroutine New_Hybrid_Step
    integer :: i
    real(rb) :: alpha1, alpha2, dt_2N, dt_4N
    dt_2N = dt_2/nloops
    dt_4N = dt_4/nloops
    if (single) then
      do i = 1, nloops
        call thermostat(1) % integrate( dt_4N, two*md%Energy%Kinetic )
        alpha1 = thermostat(1)%p(1)*thermostat(1)%InvQ(1)
        thermostat(1)%eta(1) = thermostat(1)%eta(1) + alpha1*dt_2N
        call EmDee_boost( md, one, alpha1, dt_2N )
        call thermostat(1) % integrate( dt_4N, two*md%Energy%Kinetic )
      end do
      call EmDee_displace( md, one, zero, dt )
      do i = 1, nloops
        call thermostat(1) % integrate( dt_4N, two*md%Energy%Kinetic )
        alpha1 = thermostat(1)%p(1)*thermostat(1)%InvQ(1)
        thermostat(1)%eta(1) = thermostat(1)%eta(1) + alpha1*dt_2N
        call EmDee_boost( md, one, alpha1, dt_2N )
        call thermostat(1) % integrate( dt_4N, two*md%Energy%Kinetic )
      end do
    else
      do i = 1, nloops
        call thermostat % integrate( dt_4N, two*[md%Energy%Kinetic - md%Energy%Rotational, md%Energy%Rotational] )
        alpha1 = thermostat(1)%p(1)*thermostat(1)%InvQ(1)
        thermostat(1)%eta(1) = thermostat(1)%eta(1) + alpha1*dt_2N
        transOnly
        call EmDee_boost( md, one, alpha1, dt_2N )
        alpha2 = thermostat(2)%p(1)*thermostat(2)%InvQ(1)
        thermostat(2)%eta(1) = thermostat(2)%eta(1) + alpha2*dt_2N
        rotOnly
        call EmDee_boost( md, one, alpha2, dt_2N )
        call thermostat % integrate( dt_4N, two*[md%Energy%Kinetic - md%Energy%Rotational, md%Energy%Rotational] )
      end do
      trans_rot
      call EmDee_displace( md, one, zero, dt )
      do i = 1, nloops
        call thermostat % integrate( dt_4N, two*[md%Energy%Kinetic - md%Energy%Rotational, md%Energy%Rotational] )
        alpha1 = thermostat(1)%p(1)*thermostat(1)%InvQ(1)
        thermostat(1)%eta(1) = thermostat(1)%eta(1) + alpha1*dt_2N
        transOnly
        call EmDee_boost( md, one, alpha1, dt_2N )
        alpha2 = thermostat(2)%p(1)*thermostat(2)%InvQ(1)
        thermostat(2)%eta(1) = thermostat(2)%eta(1) + alpha2*dt_2N
        rotOnly
        call EmDee_boost( md, one, alpha2, dt_2N )
        call thermostat % integrate( dt_4N, two*[md%Energy%Kinetic - md%Energy%Rotational, md%Energy%Rotational] )
      end do
      trans_rot
    end if
  end subroutine New_Hybrid_Step
  !-------------------------------------------------------------------------------------------------
  subroutine Hybrid_Step
    if (single) then
      call thermostat(1) % integrate( dt_4, two*md%Energy%Kinetic )
      call EmDee_boost( md, zero, thermostat(1)%damping, dt_4 )

      call EmDee_boost( md, one, zero, dt_2 )

      call thermostat(1) % integrate( dt_4, two*md%Energy%Kinetic )
      call EmDee_boost( md, zero, thermostat(1)%damping, dt_4 )

      call EmDee_displace( md, one, zero, dt )

      call thermostat(1) % integrate( dt_4, two*md%Energy%Kinetic )
      call EmDee_boost( md, zero, thermostat(1)%damping, dt_4 )

      call EmDee_boost( md, one, zero, dt_2 )

      call thermostat(1) % integrate( dt_4, two*md%Energy%Kinetic )
      call EmDee_boost( md, zero, thermostat(1)%damping, dt_4 )
    else
      call thermostat % integrate( dt_4, two*[md%Energy%Kinetic - md%Energy%Rotational, md%Energy%Rotational] )
      transOnly
      call EmDee_boost( md, zero, thermostat(1)%damping, dt_4 )
      rotOnly
      call EmDee_boost( md, zero, thermostat(2)%damping, dt_4 )

      trans_rot
      call EmDee_boost( md, one, zero, dt_2 )

      call thermostat % integrate( dt_4, two*[md%Energy%Kinetic - md%Energy%Rotational, md%Energy%Rotational] )
      transOnly
      call EmDee_boost( md, zero, thermostat(1)%damping, dt_4 )
      rotOnly
      call EmDee_boost( md, zero, thermostat(2)%damping, dt_4 )

      call EmDee_displace( md, one, zero, dt )

      call thermostat % integrate( dt_4, two*[md%Energy%Kinetic - md%Energy%Rotational, md%Energy%Rotational] )
      transOnly
      call EmDee_boost( md, zero, thermostat(1)%damping, dt_4 )
      rotOnly
      call EmDee_boost( md, zero, thermostat(2)%damping, dt_4 )

      trans_rot
      call EmDee_boost( md, one, zero, dt_2 )

      call thermostat % integrate( dt_4, two*[md%Energy%Kinetic - md%Energy%Rotational, md%Energy%Rotational] )
      transOnly
      call EmDee_boost( md, zero, thermostat(1)%damping, dt_4 )
      rotOnly
      call EmDee_boost( md, zero, thermostat(2)%damping, dt_4 )
      trans_rot
    end if
  end subroutine Hybrid_Step
  !-------------------------------------------------------------------------------------------------
  subroutine Boosting_Step
    real(rb) :: alpha1, alpha2
    if (single) then
      call thermostat(1) % integrate( dt_2, two*md%Energy%Kinetic )
      alpha1 = thermostat(1)%p(1)*thermostat(1)%InvQ(1)
      call EmDee_boost( md, one, alpha1, dt_2 )
      call EmDee_displace( md, one, zero, dt )
      thermostat(1)%eta(1) = thermostat(1)%eta(1) + alpha1*dt
      call EmDee_boost( md, one, alpha1, dt_2 )
      call thermostat(1) % integrate( dt_2, two*md%Energy%Kinetic )
    else
      call thermostat % integrate( dt_2, two*[md%Energy%Kinetic - md%Energy%Rotational, md%Energy%Rotational] )
      alpha1 = thermostat(1)%p(1)*thermostat(1)%InvQ(1)
      alpha2 = thermostat(2)%p(1)*thermostat(2)%InvQ(1)
      transOnly
      call EmDee_boost( md, one, alpha1, dt_2 )
      rotOnly
      call EmDee_boost( md, one, alpha2, dt_2 )
      trans_rot
      call EmDee_displace( md, one, zero, dt )
      thermostat(1)%eta(1) = thermostat(1)%eta(1) + alpha1*dt
      thermostat(2)%eta(1) = thermostat(2)%eta(1) + alpha2*dt
      transOnly
      call EmDee_boost( md, one, alpha1, dt_2 )
      rotOnly
      call EmDee_boost( md, one, alpha2, dt_2 )
      call thermostat % integrate( dt_2, two*[md%Energy%Kinetic - md%Energy%Rotational, md%Energy%Rotational] )
      trans_rot
    end if
  end subroutine Boosting_Step
  !-------------------------------------------------------------------------------------------------
  subroutine Kamberaj_Single_Thermostat_Step
    call EmDee_boost( md, one, zero, dt_2 )
    call EmDee_boost( md, zero, thermostat(1)%p(1)*thermostat(1)%InvQ(1), dt_2 )
    call EmDee_displace( md, one, zero, dt )
    call thermostat(1) % integrate( dt, two*md%Energy%Kinetic )
    call EmDee_boost( md, zero, thermostat(1)%p(1)*thermostat(1)%InvQ(1), dt_2 )
    call EmDee_boost( md, one, zero, dt_2 )
  end subroutine Kamberaj_Single_Thermostat_Step
  !-------------------------------------------------------------------------------------------------
  subroutine Kamberaj_Step
    call EmDee_boost( md, one, zero, dt_2 )
    if (single) then
      call EmDee_boost( md, zero, thermostat(1)%p(1)*thermostat(1)%InvQ(1), dt_2 )
      call EmDee_displace( md, one, zero, dt )
      call thermostat(1) % integrate( dt, two*md%Energy%Kinetic )
      call EmDee_boost( md, zero, thermostat(1)%p(1)*thermostat(1)%InvQ(1), dt_2 )
    else
      transOnly
      call EmDee_boost( md, zero, thermostat(1)%p(1)*thermostat(1)%InvQ(1), dt_2 )
      rotOnly
      call EmDee_boost( md, zero, thermostat(2)%p(1)*thermostat(2)%InvQ(1), dt_2 )
      trans_rot
      call EmDee_displace( md, one, zero, dt )
      call thermostat % integrate( dt, two*[md%Energy%Kinetic - md%Energy%Rotational, md%Energy%Rotational] )
      rotOnly
      call EmDee_boost( md, zero, thermostat(2)%p(1)*thermostat(2)%InvQ(1), dt_2 )
      transOnly
      call EmDee_boost( md, zero, thermostat(1)%p(1)*thermostat(1)%InvQ(1), dt_2 )
      trans_rot
    end if
    call EmDee_boost( md, one, zero, dt_2 )
  end subroutine Kamberaj_Step
  !-------------------------------------------------------------------------------------------------
end program lj_nvt
