#include "emdee.f03"

program lj_nvt

use mGlobal
use mConfig
use mRandom
use mThermostat
use iso_c_binding
use EmDee
use mCorrelate

implicit none

real(rb), parameter :: mvv2e = 2390.057364_rb         ! Da*A²/fs² to kcal/mol
real(rb), parameter :: kB = 8.31451E-7_rb             ! Boltzmann constant in Da*A²/(fs²*K)
real(rb), parameter :: Pconv = 1.6388244954E+8_rb     ! Da/(A*fs²) to atm
real(rb), parameter :: kCoul = 0.13893545755135628_rb ! Coulomb constant in Da*A³/(fs*e)²

! Simulation specifications:
character(sl) :: Base
integer       :: i, N, NB, seed, MDsteps, Nconf, thermo, Nequil, Nprod, rotationMode
real(rb)      :: T, Rc, Rm, dt, skin, alpha

! System properties:
integer  :: dof
real(rb) :: Volume

! Thermostat variables:
integer :: method, M, ndamp, nloops, nts
logical :: single
class(nhc), pointer :: thermostat(:)

! Mean Square Displacement variables: 
real(rb), pointer :: Rcm(:,:)
integer :: DoMsd, nevery, blocksize, nfreq

!Dipole moment variables:
integer :: DoDipole, Dnevery, Dblocksize, Dnfreq, NP
real(rb), pointer :: q(:,:), R(:,:)
real(rb), allocatable :: mu(:,:), theta(:,:) 

!Radial distribution function variables:
integer :: DoRdf, Gnevery, Rnfreq, bins, npairs, counter
real(rb), allocatable :: gr(:,:), rdf(:,:)
integer , allocatable :: itype(:), jtype(:)

! Other variables:
integer  :: step
real(rb) :: dt_2, dt_4, KE_sp, kT
character(256) :: filename, configFile

integer :: threads
type(tEmDee) :: md
type(c_ptr), allocatable :: model(:)
type(kiss) :: random
type(tMSD) :: MSD
type(tACF) :: ACF
integer :: out
#define trans_rot md%Options%translate = .true.;  md%Options%rotate = .true.
#define transOnly md%Options%translate = .true.;  md%Options%rotate = .false.
#define rotOnly   md%Options%translate = .false.; md%Options%rotate = .true.

character(*), parameter :: titles = "Step rank Temp Ts Press Ps KinEng KinEng_t KinEng_r "// &
                                    "KinEng_r1 KinEng_r2 KinEng_r3 DispEng CoulEng PotEng "// &
                                    "TotEng H_nhc Virial BodyVirial Ks Ks_t Ks_r Us Hs Hs_nhc"
! MPI variables:
include 'mpif.h'
integer :: ierr, nprocs, my_rank
character(sl) :: rank

! Executable code:
#ifdef coul
  call writeln( "md/lj/coul ("//__DATE__//")" )
#else
  call writeln( "md/lj ("//__DATE__//")" )
#endif

call Get_Command_Line_Args( threads, filename )

call MPI_Init( ierr )
call MPI_Comm_Size( MPI_COMM_WORLD, nprocs, ierr )
call MPI_Comm_Rank( MPI_COMM_WORLD, my_rank, ierr )
rank = adjustl(int2str(my_rank + 1))
filename = trim(filename)//"_"//trim(rank)//".inp"

call Read_Specifications( filename )

call Config % Read( configFile )
call Setup_Simulation

allocate( model(Config%ntypes) )
call Configure_System( md, Rm, Rc, alpha )
call Config % Save_XYZ( trim(Base)//".xyz" )

call writeln( titles )
step = 0
call writeln( properties() )
do step = 1, NEquil
  call execute_step
  if (mod(step,thermo) == 0) call writeln( properties() )
end do
call writeln( "Loop time of", real2str(md%Time%Total), "s." )
call Report( md )

if (my_rank == 0) call writeln( )
call writeln( "Memory usage" )
call writeln( titles )
step = NEquil
call writeln( properties() )

if (DoMsd == 1) then
  allocate(Rcm(3,NB)) 
  call EmDee_download( md, "centersOfMass"//c_null_char, c_loc(Rcm(1,1)) )  
  call MSD % setup( nevery, blocksize, Nprod, Rcm )
  open(newunit = out, file = trim(Base)//".msd", status = "replace")
  close(out)
end if

if (DoDipole == 1) then
  allocate(q(4,NB)) 
  allocate(R(3,N)) 
  allocate(mu(3,NB)) 
  allocate(theta(3,NB)) 
  NP = N/NB  
  call EmDee_download( md, "coordinates"//c_null_char, c_loc(R(1,1)) )
  call EmDee_download(md, "quaternions"//c_null_char, c_loc(q(1,1)))
  call ComputeTheta
  call ACF % setup( Dnevery, Dblocksize, Nprod, mu )
  open(newunit = out, file = trim(Base)//".dipole", status = "replace")
  close(out)
end if

if (DoRdf == 1) then 
  allocate(gr(bins,npairs), source = 0.0_rb) 
  allocate(rdf(bins,npairs), source = 0.0_rb)
  open(newunit = out, file = trim(Base)//".rdf", status = "replace")
  counter = 1   
  call EmDee_rdf(md, bins, npairs, itype, jtype, rdf)
  gr = rdf
  close(out)
end if

do step = NEquil+1, NEquil+NProd
  call execute_step
  if (mod(step,Nconf)==0) then
    call EmDee_download( md, "coordinates"//c_null_char, c_loc(Config%R(1,1)) )
    call Config % Save_XYZ( trim(Base)//".xyz", append = .true. )
  end if
 if (mod(step,thermo) == 0) call writeln( properties() ) 
 if (DoMSD == 1 .AND. mod(step,nevery) == 0) then
    call EmDee_download( md, "centersOfMass"//c_null_char, c_loc(Rcm(1,1)) )
    call MSD % sample( Rcm ) 
    if (mod(step,nfreq) == 0) then
      call MSD % save( trim(Base)//".msd", append = .true. )
    end if
  end if
  if (DoDipole == 1 .AND. mod(step,Dnevery) == 0)  then 
    call EmDee_download(md, "quaternions"//c_null_char, c_loc(q(1,1)))
    call ComputeMu
    call ACF % sample( mu ) 
    if (mod(step,Dnfreq) == 0) then
      call ACF % save( trim(Base)//".dipole", append = .true. )
    end if
  end if
  if (DoRdf == 1 .AND. mod(step,Gnevery) == 0)  then 
    call EmDee_rdf(md, bins, npairs, itype, jtype, rdf)
    gr = gr + rdf
    counter = counter + 1   
    if (mod(step,Rnfreq) == 0) then
      call rdf_save_to_file( trim(Base)//".rdf", append = .true. )
    end if
  end if
end do
call writeln( "Loop time of", real2str(md%Time%Total), "s." )
call Report( md )
call EmDee_download( md, "coordinates"//c_null_char, c_loc(Config%R(1,1)) )
call Config % Write( trim(Base)//"_out.lmp", velocities = .true. )
call stop_log
call MPI_Finalize( ierr )

contains
  !-------------------------------------------------------------------------------------------------
  subroutine execute_step
    select case (method)
      case (0); call EmDee_verlet_step( md, dt )
      case (1); call Pscaling_Step
      case (2); call Boosting_Step
      case (3); call Kamberaj_Step
      case (4); call Hybrid_Step
      case (5); call New_Hybrid_Step
      case (6); call Pscaling_Shadow_Step
    end select
  end subroutine execute_step
  !-------------------------------------------------------------------------------------------------
  subroutine Configure_System( md, Rm, Rc, alpha )
    type(tEmDee), intent(inout) :: md
    real(rb),     intent(in)    :: Rm, Rc, alpha

    md = EmDee_system( threads, 1, Rc, skin, Config%natoms, &
                       c_loc(Config%Type(1)), c_loc(Config%mass(1)), c_loc(Config%Mol(1)) )
    md%Options%rotationMode = rotationMode

    do i = 1, Config%ntypes
      if (abs(Config%epsilon(i)) < epsilon(1.0_rb)) then
        model(i) = EmDee_pair_none()
      else
        model(i) = EmDee_smoothed( &
                     EmDee_pair_lj_cut( Config%epsilon(i)/mvv2e, Config%sigma(i) ), Rc-Rm )
      end if
      call EmDee_set_pair_model( md, i, i, model(i), kCoul )
    end do

#   ifdef coul
      call EmDee_set_coul_model( md, EmDee_coul_damped_smoothed( alpha, Rc-Rm ) )
      call EmDee_upload( md, "charges"//c_null_char, c_loc(Config%Charge(1)) )
#   endif

    call EmDee_upload( md, "box"//c_null_char, c_loc(Config%Lx) )
    call EmDee_upload( md, "coordinates"//c_null_char, c_loc(Config%R(1,1)) )
    call EmDee_random_momenta( md, kT, .true._1, seed )

  end subroutine Configure_System
  !-------------------------------------------------------------------------------------------------
  subroutine Report( md )
    type(tEmDee), intent(in) :: md
    real(rb) :: other
    call writeln( repeat("-",40) )
    call writeln( "Pair time      =", real2str(md%Time%Pair), "s." )
    call writeln( "Neighbor time  =", real2str(md%Time%Neighbor), "s." )
    other = md%Time%Total - (md%Time%Pair + md%Time%Neighbor)
    call writeln( "Neighbor list builds =", int2str(md%builds) )
    call writeln( repeat("-",40) )
  end subroutine Report
  !-------------------------------------------------------------------------------------------------
  character(sl) function properties()
    real(rb) :: Temp, Ts, H, Hs, Hthermo
    Temp = (md%Energy%Kinetic/KE_sp)*T
    Ts = (md%Energy%ShadowKinetic/KE_sp)*T
    H = md%Energy%Potential + md%Energy%Kinetic
    Hs = md%Energy%ShadowPotential + md%Energy%ShadowKinetic
    Hthermo = sum(thermostat%energy())
    properties = trim(adjustl(int2str(step))) // " " // trim(rank) // " " // &
                 join(real2str([ Temp, &
                                 Ts, &
                                 Pconv*((NB-1)*kB*Temp + md%Virial/3.0_rb)/Volume, &
                                 Pconv*((NB-1)*kB*Ts + md%Virial/3.0_rb)/Volume, &
                                 mvv2e*[md%Energy%Kinetic, &
                                        md%Energy%Kinetic - md%Energy%Rotational, &
                                        md%Energy%Rotational, &
                                        md%Energy%RotPart, &
                                        md%Energy%Dispersion, &
                                        md%Energy%Coulomb, &
                                        md%Energy%Potential, &
                                        H, &
                                        H + Hthermo, &
                                        md%Virial, &
                                        md%BodyVirial, &
                                        md%Energy%ShadowKinetic, &
                                        md%Energy%ShadowKinetic - md%Energy%ShadowRotational, &
                                        md%Energy%ShadowRotational, &
                                        md%Energy%ShadowPotential, &
                                        Hs, &
                                        Hs + Hthermo]]))
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
    integer :: inp, i
    open( newunit=inp, file = file, status = "old" )
    read(inp,*); read(inp,*) Base
    read(inp,*); read(inp,*) configFile
    read(inp,*); read(inp,*) T
    read(inp,*); read(inp,*) Rc, Rm, alpha
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
    read(inp,*); read(inp,*) DoMSD, nevery, blocksize, nfreq
    read(inp,*); read(inp,*) DoDipole, Dnevery, Dblocksize, Dnfreq
    read(inp,*); read(inp,*) npairs
    allocate( itype(npairs), jtype(npairs) )
    read(inp,*); read(inp,*) DoRdf, Gnevery, Rnfreq, bins, (itype(i),jtype(i),i=1,npairs)
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
      case (0,1,4,6); allocate( nhc_pscaling :: thermostat(nts) )
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
  subroutine Pscaling_Step
    if (single) then
      call thermostat(1) % integrate( dt_2, two*md%Energy%Kinetic )
      call EmDee_boost( md, zero, thermostat(1)%damping, dt_2 )
      call EmDee_verlet_step( md, dt )
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
  subroutine Pscaling_Shadow_Step
    real(rb) :: alpha, factor
    if (single) then
      call thermostat(1) % integrate( dt_2, two*md%Energy%ShadowKinetic )
      call EmDee_boost( md, zero, thermostat(1)%damping, dt_2 )
      call EmDee_verlet_step( md, dt )
      call thermostat(1) % integrate( dt_2, two*md%Energy%ShadowKinetic )
      alpha = thermostat(1)%damping
      call EmDee_boost( md, zero, alpha, dt_2 )
      factor = exp(-alpha*dt)
      md%Energy%ShadowKinetic = factor*md%Energy%ShadowKinetic
      md%Energy%ShadowRotational = factor*md%Energy%ShadowRotational
    else
      stop "P-scaling shadow only admits single thermostat"
    end if
  end subroutine Pscaling_Shadow_Step
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
!-----------------------------------------------------------------------------------------------------
  subroutine Separate_Boost_Step
    transOnly
    call EmDee_boost( md, one, zero, dt_4 ) 
    rotOnly
    call EmDee_boost( md, one, zero, dt_4 )
    trans_rot
    call EmDee_displace( md, one, zero, dt_2 )    
    rotOnly
    call EmDee_boost( md, one, zero, dt_2 )
    transOnly
    call EmDee_boost( md, one, zero, dt_2 )
    trans_rot
    call EmDee_displace( md, one, zero, dt_2 )    
    rotOnly
    call EmDee_boost( md, one, zero, dt_4 )
    transOnly
    call EmDee_boost( md, one, zero, dt_4 )     
   end subroutine Separate_Boost_Step
!-----------------------------------------------------------------------------------------------------
  subroutine ComputeTheta
    integer :: ind, i, j
    real(rb) :: A(3,3)
    ind = 1
    do i = 1, NB
      do j = 1, NP       
        mu(:,i) = Config%Charge(ind)*R(:,ind) 
        ind = ind + 1 
      end do
      A =  matmul(matrix_Bt(q(:,i)),matrix_C(q(:,i)))
      theta(:,i) = matmul(A,mu(:,i))
    end do
  end subroutine ComputeTheta
!---------------------------------------------------------------------------------------------------
  subroutine ComputeMu
    integer ::  i
    real(rb) :: A(3,3)
    do i = 1, NB
      A =  matmul(matrix_Bt(q(:,i)),matrix_C(q(:,i)))      
      mu(:,i) = matmul( transpose(A),theta(:,i) )
    end do
  end subroutine ComputeMu
!---------------------------------------------------------------------------------------------------
  subroutine rdf_save_to_file(  file, append )
    character(*),         intent(in)           :: file
    logical,        intent(in), optional :: append
    integer, parameter :: out = 65
    logical :: app
    app = present(append)
    if (app) app = append
    if (app) then
      open(unit=out,file=file,status="old",position="append")
    else
      open(unit=out,file=file,status="replace")
    end if
    call rdf_save_to_unit( out )
    close(out)
  end subroutine rdf_save_to_file
!---------------------------------------------------------------------------------------------------
 subroutine rdf_save_to_unit( unit)
   integer,        intent(in)    :: unit
   integer :: i
   character(3) :: Ci, Cj
   character(sl) :: title
   title = "r"
   do i = 1, npairs
     write(Ci,'(I3)') itype(i)
     write(Cj,'(I3)') jtype(i)
     title = trim(title)//" g("//trim(adjustl(Ci))//","//trim(adjustl(Cj))//")"
   end do
   write(unit,'(A)') trim(title)
     do i = 1, bins 
       write(unit,*) (i-0.5)*Rc/bins, gr(i,:)/real(counter,rb)
     end do
 end subroutine rdf_save_to_unit
!---------------------------------------------------------------------------------------------------
  pure function matrix_B( q ) result( B )
    real(rb), intent(in) :: q(0:3)
    real(rb)             :: B(4,3)
    B = reshape( [-q(1),  q(0),  q(3), -q(2),  &
                  -q(2), -q(3),  q(0),  q(1),  &
                  -q(3),  q(2), -q(1),  q(0)], [4,3] )
  end function matrix_B
!---------------------------------------------------------------------------------------------------
  pure function matrix_C( q ) result( C )
    real(rb), intent(in) :: q(0:3)
    real(rb)             :: C(4,3)
    C = reshape( [-q(1),  q(0), -q(3),  q(2),  &
                  -q(2),  q(3),  q(0), -q(1),  &
                  -q(3), -q(2),  q(1),  q(0)], [4,3] )
  end function matrix_C
!---------------------------------------------------------------------------------------------------
  pure function matrix_Bt( q ) result( Bt )
    real(rb), intent(in) :: q(0:3)
    real(rb)             :: Bt(3,4)
    Bt = reshape( [-q(1), -q(2), -q(3), &
                    q(0), -q(3),  q(2), &
                    q(3),  q(0), -q(1), &
                   -q(2),  q(1),  q(0)  ], [3,4] )
  end function matrix_Bt
!---------------------------------------------------------------------------------------------------
  pure function matrix_Ct( q ) result( Ct )
    real(rb), intent(in) :: q(0:3)
    real(rb)             :: Ct(3,4)
    Ct = reshape( [-q(1), -q(2), -q(3), &
                    q(0),  q(3), -q(2), &
                   -q(3),  q(0),  q(1), &
                    q(2), -q(1),  q(0)  ], [3,4] )
  end function matrix_Ct
!===================================================================================================
end program lj_nvt
