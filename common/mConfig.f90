module mConfig

use mGlobal

implicit none

character, parameter :: element(5) = ["H","C","N","O","S"]
integer,   parameter :: atomic_mass(5) = [1,12,14,16,32]

type tConfig

  real(rb) :: xmin, xmax, ymin, ymax, zmin, zmax
  real(rb), pointer :: Lx, Ly, Lz

  integer  :: ntypes
  real(rb), allocatable :: epsilon(:), sigma(:)

  integer  :: natoms
  integer,  pointer :: Mol(:), Type(:)
  real(rb), pointer :: Charge(:)

  real(rb), pointer :: R(:,:), F(:,:), P(:,:)

  real(rb), pointer :: mass(:)
  real(rb), pointer :: Rx(:), Ry(:), Rz(:) ! Positions
  real(rb), pointer :: Fx(:), Fy(:), Fz(:) ! Forces
  real(rb), pointer :: Px(:), Py(:), Pz(:) ! Linear momenta
  real(rb), allocatable :: InvMass(:)

  logical :: velocity_input = .false.

  contains

    procedure :: tConfig_Read, tConfig_Read_from_File
    generic :: Read => tConfig_Read, tConfig_Read_from_File

    procedure :: tConfig_Write, tConfig_Write_to_File
    generic :: Write => tConfig_Write, tConfig_Write_to_File

    procedure :: tConfig_Save_XYZ_to_unit, tConfig_Save_XYZ_to_file
    generic :: Save_XYZ => tConfig_Save_XYZ_to_unit, tConfig_Save_XYZ_to_file

end type tConfig

type(tConfig) :: Config

contains

  !=================================================================================================

  subroutine tConfig_Read( me, unit )
    class(tConfig), intent(inout) :: me
    integer,        intent(in)    :: unit
    integer       :: narg, i, k
    character(sl) :: arg(10)
    real(rb) :: mass
    call next_command( unit, narg, arg )
    do while (narg > 0)
      call next_command( unit, narg, arg )

      if ((narg == 3).and.(join(arg(2:3)) == "atom types")) then
        me % ntypes = str2int( arg(1) )

      else if ((narg == 2).and.(arg(2) == "atoms")) then
        me % natoms = str2int( arg(1) )

      else if ((narg == 4).and.(join(arg(3:4)) == "xlo xhi")) then
        me % xmin = str2real(arg(1))
        me % xmax = str2real(arg(2))
        allocate( me % Lx )
        me % Lx = me % xmax - me % xmin

      else if ((narg == 4).and.(join(arg(3:4)) == "ylo yhi")) then
        me % ymin = str2real(arg(1))
        me % ymax = str2real(arg(2))
        allocate( me % Ly )
        me % Ly = me % ymax - me % ymin

      else if ((narg == 4).and.(join(arg(3:4)) == "zlo zhi")) then
        me % zmin = str2real(arg(1))
        me % zmax = str2real(arg(2))
        allocate( me % Lz )
        me % Lz = me % zmax - me % zmin

      else if ((narg == 1).and.(arg(1) == "Masses")) then
        allocate( me % mass(me % ntypes) )
        do i = 1, me % ntypes
          call next_command( unit, narg, arg )
          me % mass(i) = str2real(arg(2))
        end do

      else if ((narg == 2).and.(join(arg(1:2)) == "Pair Coeffs")) then
        allocate( me % epsilon(me % ntypes), me % sigma(me % ntypes) )
        do i = 1, me % ntypes
          call next_command( unit, narg, arg )
          me % epsilon(i) = str2real(arg(2))
          me % sigma(i) = str2real(arg(3))
        end do

      else if ((narg == 1).and.(arg(1) == "Atoms")) then
        associate( N => me % natoms )
          allocate( me%Mol(N), me%Type(N), me%Charge(N) )

          allocate( me%R(3,N), me%F(3,N), me%P(3,N) )
          me%Rx => me%R(1,:)
          me%Ry => me%R(2,:)
          me%Rz => me%R(3,:)
          me%Fx => me%F(1,:)
          me%Fy => me%F(2,:)
          me%Fz => me%F(3,:)
          me%Px => me%P(1,:)
          me%Py => me%P(2,:)
          me%Pz => me%P(3,:)

!          allocate( me%Rx(N), me%Ry(N), me%Rz(N) )
!          allocate( me%Fx(N), me%Fy(N), me%Fz(N) )
!          allocate( me%Px(N), me%Py(N), me%Pz(N) )
          allocate( me%InvMass(N) )
        end associate
        do k = 1, me % natoms
          call next_command( unit, narg, arg )
          i = str2int(arg(1))
          me % Mol(i) = str2int(arg(2))
          me % Type(i) = str2int(arg(3))
          me % Charge(i) = str2real(arg(4))
          me % Rx(i) = str2real(arg(5))
          me % Ry(i) = str2real(arg(6))
          me % Rz(i) = str2real(arg(7))
          me % Px(i) = zero
          me % Py(i) = zero
          me % Pz(i) = zero
          if (narg == 10) then
            me % Rx(i) = me % Rx(i) + str2real(arg( 8)) * me % Lx
            me % Ry(i) = me % Ry(i) + str2real(arg( 9)) * me % Ly
            me % Rz(i) = me % Rz(i) + str2real(arg(10)) * me % Lz
          end if
        end do
        me % InvMass = one/me%Mass(me%Type)

      else if ((narg == 1).and.(arg(1) == "Velocities")) then

        me % velocity_input = .true.
        do k = 1, me % natoms
          call next_command( unit, narg, arg )
          i = str2int(arg(1))
          mass = me%Mass(me%Type(i))
          me % Px(i) = mass*str2real(arg(2))
          me % Py(i) = mass*str2real(arg(3))
          me % Pz(i) = mass*str2real(arg(4))
        end do

      end if
    end do
  end subroutine tConfig_Read

  !=================================================================================================

  subroutine tConfig_Read_from_File( me, file )
    class(tConfig), intent(inout) :: me
    character(*),   intent(in)    :: file
    integer :: inp, stat
    inp = 69
    open( unit = inp, file = file, status = "old", iostat = stat )
    if (stat /= 0) call error( "Configuration file", trim(file), "was not found." )
    call me % Read( inp )
    close(inp)
  end subroutine tConfig_Read_from_File

  !=================================================================================================

  subroutine tConfig_Write( me, unit, velocities )
    class(tConfig), intent(inout) :: me
    integer,        intent(in)    :: unit
    logical,        intent(in), optional :: velocities
    integer :: i
    real(rb) :: xi, yi, zi, mass
    write(unit,'("LAMMPS data file",/)')
    write(unit,'(A," atom types",/)') trim(int2str(me % ntypes))
    write(unit,'(A," atoms",/)') trim(int2str(me % natoms))
    write(unit,'(A," xlo xhi")') trim(join(real2str([me%xmin,me%xmax])))
    write(unit,'(A," ylo yhi")') trim(join(real2str([me%ymin,me%ymax])))
    write(unit,'(A," zlo zhi")') trim(join(real2str([me%zmin,me%zmax])))
    write(unit,'(/,"Masses",/)')
    do i = 1, me % ntypes
      write(unit,'(A)') trim(join([int2str(i),real2str(me%mass(i))]))
    end do
    write(unit,'(/,"Pair Coeffs",/)')
    do i = 1, me % ntypes
      write(unit,'(A)') trim(join([int2str(i),real2str(me%epsilon(i)),real2str(me%sigma(i))]))
    end do
    write(unit,'(/,"Atoms",/)')
    do i = 1, me % natoms
!      xi = (me%Rx(i) - me%xmin)/me%Lx
!      yi = (me%Ry(i) - me%ymin)/me%Ly
!      zi = (me%Rz(i) - me%zmin)/me%Lz
!      xi = me%xmin + me%Lx*(xi - floor(xi))
!      yi = me%ymin + me%Ly*(yi - floor(yi))
!      zi = me%zmin + me%Lz*(zi - floor(zi))
      xi = me%Rx(i)
      yi = me%Ry(i)
      zi = me%Rz(i)
      write(unit,'(A,X,A)') trim(join(int2str([i,me%Mol(i),me%Type(i)]))),  &
                            trim(join(real2str([me%Charge(i),xi,yi,zi])))
    end do
    if (present(velocities)) then
      if (velocities) then
        write(unit,'(/,"Velocities",/)')
        do i = 1, me % natoms
          mass = me%mass(me%Type(i))
          xi = me%Px(i)/mass
          yi = me%Py(i)/mass
          zi = me%Pz(i)/mass
          write(unit,'(A,X,A)') trim(int2str(i)), trim(join(real2str([xi,yi,zi])))
        end do
      end if
    end if
  end subroutine tConfig_Write

  !=================================================================================================

  subroutine tConfig_Write_to_File( me, file, velocities )
    class(tConfig), intent(inout) :: me
    character(*),   intent(in)    :: file
    logical,        intent(in), optional :: velocities
    integer :: out, stat
    out = 69
    open( unit = out, file = file, status = "replace", iostat = stat )
    if (stat /= 0) call error( "Cannot open file", trim(file), "for writing." )
    if (present(velocities)) then
      call me % Write( out, velocities )
    else
      call me % Write( out )
    end if
    close(out)
  end subroutine tConfig_Write_to_File

  !=================================================================================================

  subroutine tConfig_Save_XYZ_to_file( me, file, append )
    class(tConfig), intent(inout)        :: me
    character(*),   intent(in)           :: file
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
    call tConfig_Save_XYZ_to_unit( me, out )
    close(out)
  end subroutine tConfig_Save_XYZ_to_file

  !=================================================================================================

  subroutine tConfig_Save_XYZ_to_unit( me, unit )
    class(tConfig), intent(inout) :: me
    integer,        intent(in)    :: unit
    integer  :: i, m, nmol, imol
    real(rb) :: imass, ri(3), L(3), Lmin(3)
    real(rb), allocatable :: rcm(:,:), molMass(:)
    character(sl) :: itype
    write(unit,*) me % natoms
    write(unit,*)
    nmol = maxval(me%mol)
    allocate( rcm(3,nmol), molMass(nmol) )
    rcm = 0.0_rb
    molMass = 0.0_rb
    do i = 1, me%natoms
      imol = me%mol(i)
      imass = me%mass(me%Type(i))
      molMass(imol) = molMass(imol) + imass
      rcm(:,imol) = rcm(:,imol) + imass*me%R(:,i)
    end do
    forall (imol=1:nmol) rcm(:,imol) = rcm(:,imol)/molMass(imol)
    L = [me%Lx, me%Ly, me%Lz]
    Lmin = [me%xmin, me%ymin, me%zmin]
    do imol = 1, nmol
      ri = rcm(:,imol)/L
      rcm(:,imol) = L*(ri - floor(ri))
    end do
    do i = 1, me % natoms
      imol = me%mol(i)
      ri = me%R(:,i) - L*anint((me%R(:,i) - rcm(:,imol))/L)
      m = nint(me%Mass(me%Type(i)))
      if (any(atomic_mass == m)) then
        itype = element(maxloc(transfer(atomic_mass == m, atomic_mass),dim=1))
      else
        itype = int2str(me%Type(i))
      end if
      write(unit,*) trim(itype), ri
    end do
  end subroutine tConfig_Save_XYZ_to_unit

  !=================================================================================================

end module mConfig
