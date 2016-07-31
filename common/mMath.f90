!> This module contains general mathematical functions and subroutines.
module mMath
use mGlobal
implicit none

integer, parameter :: ib = 16

contains
  !-------------------------------------------------------------------------------------------------
  function gauss_elim( a, b ) result( x )
    real(rb), intent(in) :: a(:,:), b(:)
    real(rb)             :: x(size(b))
    real(rb), allocatable :: U(:,:)
    integer :: n, i, k, ip
    n = size(b)
    allocate( U(n,n+1) )
    U(:,1:n) = a
    U(:,n+1) = b
    do k = 1, n-1
      ip = maxloc(abs(U(k:n,k)),1) + k - 1
      call swap( U(k,k:), U(ip,k:) )
      do i = k+1, n
        U(i,k+1:) = U(i,k+1:) - (U(i,k)/U(k,k))*U(k,k+1:)
      end do
    end do
    do i = n, 1, -1
      x(i) = (U(i,n+1) - sum(U(i,i+1:n)*x(i+1:n)))/U(i,i)
    end do
    deallocate( U )
    contains
      elemental subroutine swap( a, b )
        real(rb), intent(inout) :: a, b
        real(rb) :: c
        c = a
        a = b
        b = c
      end subroutine swap
  end function
  !-------------------------------------------------------------------------------------------------
  pure recursive function fac( n ) result( c )
    integer, intent(in) :: n
    integer(ib)         :: c
    if (n == 0) then
      c = 1
    else
      c = n*fac(n-1)
    end if
  end function fac
  !-------------------------------------------------------------------------------------------------
  pure recursive function fac_ratio( n, m ) result( c )
    integer, intent(in) :: n, m
    integer(ib)         :: c
    if (n == m) then
      c = 1
    else
      c = n*fac_ratio(n-1,m)
    end if
  end function fac_ratio
  !-------------------------------------------------------------------------------------------------
  pure recursive function binom( n, k ) result( c )
    integer, intent(in) :: n, k
    integer(ib)         :: c
    if ((k == 0).or.(k == n)) then
      c = 1
    else
      c = binom(n-1,k-1) + binom(n-1,k)
    end if
  end function binom
  !-------------------------------------------------------------------------------------------------
  pure recursive function stir1( n, k ) result( c )
    integer, intent(in) :: n, k
    integer(ib)         :: c
    if ((n == 0).and.(k == 0)) then
      c = 1
    elseif ((n == 0).or.(k == 0)) then
      c = 0
    else
      c = (n-1)*stir1(n-1,k) + stir1(n-1,k-1)
    end if
  end function stir1
  !-------------------------------------------------------------------------------------------------
  recursive function comb( n, k ) result( c )
    integer, intent(in) :: n, k
    integer             :: c
    if ((k == 0).or.(k == n)) then
      c = 1
    else
      c = comb(n-1,k-1) + comb(n-1,k)
    end if
  end function comb
  !-------------------------------------------------------------------------------------------------
  !> Computes weights for finite difference formulas using the Fornberg (1988) algorithm.
  !!
  !! @param[in] z (scalar, real) location where approximations are to be accurate
  !! @param[in] x (vector, real) grid point locations
  !! @param[in] m (scalar, integer) maximum order of derivative for which weights are sought
  function finite_diff_weights( z, x, m ) result( c )
    real(rb), intent(in)  :: z
    integer,  intent(in)  :: m
    real(rb), intent(in)  :: x(0:)
    real(rb)              :: c(0:size(x)-1, 0:m)
    real(rb), parameter :: tol = 1.e-15_rb
    real(rb) :: c1, c2, c3, c4, c5
    integer :: i, j, k, mn
    c1 = 1.0_rb
    c4 = x(0) - z
    c = 0.0_rb
    c(0,0) = 1.0_rb
    do i = 1, size(x)-1
      mn = min(i,m)
      c2 = 1.0_rb
      c5 = c4
      c4 = x(i) - z
      do j = 0, i-1
        c3 = x(i) - x(j)
        c2 = c2*c3
        if (j == i-1) then
          do k = mn, 1, -1
            c(i,k) = c1*(k*c(i-1,k-1) - c5*c(i-1,k))/c2
          end do
          c(i,0) = -c1*c5*c(i-1,0)/c2
        endif
        do k = mn, 1, -1
          c(j,k) = (c4*c(j,k) - k*c(j,k-1))/c3
        end do
        c(j,0) = c4*c(j,0)/c3
      end do
      c1 = c2
    end do
    where(abs(c) < tol) c = 0.0_rb
  end function finite_diff_weights
  !-------------------------------------------------------------------------------------------------
  pure function product_matrix_vector( A, b ) result( c )
    real(rb), intent(in) :: A(:,:), b(:)
    real(rb)             :: c(size(b))
    c = matmul(A,b)
  end function product_matrix_vector
  !-------------------------------------------------------------------------------------------------
  pure function product_transpose_matrix_vector( A, b ) result( c )
    real(rb), intent(in) :: A(:,:), b(:)
    real(rb)             :: c(size(b))
    c = matmul(transpose(A),b)
  end function product_transpose_matrix_vector
  !-------------------------------------------------------------------------------------------------
  function eigenvalues( Matrix ) result( Lambda )
    real(rb), intent(in) :: Matrix(3,3)
    real(rb)             :: Lambda(3)
    real(rb), parameter :: SQRT3 = 1.7320508075688773_rb
    real(rb) :: A(3,3), M, C1, C0
    real(rb) :: DE, DD, EE, FF
    real(rb) :: P, SQRTP, Q, C, S, PHI
    A = Matrix
    DE = A(1,2)*A(2,3)
    DD = A(1,2)**2
    EE = A(2,3)**2
    FF = A(1,3)**2
    M  = A(1,1) + A(2,2) + A(3,3)
    C1 = ( A(1,1)*A(2,2) + A(1,1)*A(3,3) + A(2,2)*A(3,3) ) - (DD + EE + FF)
    C0 = A(3,3)*DD + A(1,1)*EE + A(2,2)*FF - A(1,1)*A(2,2)*A(3,3) - 2.0_rb * A(1,3)*DE
    P = M**2 - 3.0_rb * C1
    Q = M*(P - 1.5_rb*C1) - 13.5_rb*C0
    SQRTP = sqrt(abs(P))
    PHI = 27.0_rb*(0.25_rb * C1**2 * (P - C1) + C0 * (Q + 6.75_rb*C0))
    PHI = atan2(sqrt(abs(PHI)),Q)/3.0_rb
    C = SQRTP*cos(PHI)
    S = SQRTP*sin(PHI)/SQRT3
    Lambda(2) = (M - C)/3.0_rb
    Lambda(3) = Lambda(2) + S
    Lambda(1) = Lambda(2) + C
    Lambda(2) = Lambda(2) - S
  end function eigenvalues
  !-------------------------------------------------------------------------------------------------
  function eigenvectors( Matrix, W ) result( Q )
    real(rb), intent(in) :: Matrix(3,3), W(3)
    real(rb)             :: Q(3,3)
    real(rb), parameter :: EPS = 2.2204460492503131e-16_rb
    real(rb) :: A(3,3), NORM, N1, N2, N1TMP, N2TMP
    real(rb) :: THRESH, ERROR, WMAX, F, T
    integer :: I, J
    logical :: SUCCESS
    A = Matrix
    WMAX   = max(abs(W(1)), abs(W(2)), abs(W(3)))
    THRESH = (8.0_rb * EPS * WMAX)**2
    N1TMP = A(1,2)**2 + A(1,3)**2
    N2TMP = A(1,2)**2 + A(2,3)**2
    Q(1,1) = A(1,2) * A(2,3) - A(1,3) * A(2,2)
    Q(1,2) = Q(1,1)
    Q(2,1) = A(1,3) * A(1,2) - A(2,3) * A(1,1)
    Q(2,2) = Q(2,1)
    Q(3,2) = A(1,2)**2
    A(1,1) = A(1,1) - W(1)
    A(2,2) = A(2,2) - W(1)
    Q(1,1) = Q(1,2) + A(1,3) * W(1)
    Q(2,1) = Q(2,2) + A(2,3) * W(1)
    Q(3,1) = A(1,1) * A(2,2) - Q(3,2)
    NORM = Q(1,1)**2 + Q(2,1)**2 + Q(3,1)**2
    N1 = N1TMP + A(1,1)**2
    N2 = N2TMP + A(2,2)**2
    ERROR = N1 * N2
    if (N1 <= THRESH) then
      Q(1,1) = 1.0_rb
      Q(2,1) = 0.0_rb
      Q(3,1) = 0.0_rb
    else if (N2 <= THRESH) then
      Q(1,1) = 0.0_rb
      Q(2,1) = 1.0_rb
      Q(3,1) = 0.0_rb
    else if (NORM < (64.0_rb * EPS)**2 * ERROR) then
      T = abs(A(1,2))
      F = -A(1,1) / A(1,2)
      if (abs(A(2,2)) > T) then
        T = abs(A(2,2))
        F = -A(1,2) / A(2,2)
      end if
      if (abs(A(2,3)) < T) then
        F = -A(1,3) / A(2,3)
      end if
      NORM = 1.0_rb / sqrt(1.0_rb + F**2)
      Q(1,1) = NORM
      Q(2,1) = F * NORM
      Q(3,1) = 0.0_rb
    else
      NORM = sqrt(1.0_rb / NORM)
      Q(:,1) = Q(:,1) * NORM
    end if
    T = W(1) - W(2)
    if (abs(T) < 8.0_rb * EPS * WMAX) then
      A(1,1) = A(1,1) + T
      A(2,2) = A(2,2) + T
      Q(1,2) = Q(1,2) + A(1,3) * W(2)
      Q(2,2) = Q(2,2) + A(2,3) * W(2)
      Q(3,2) = A(1,1) * A(2,2) - Q(3,2)
      NORM = Q(1,2)**2 + Q(2,2)**2 + Q(3,2)**2
      N1 = N1TMP + A(1,1)**2
      N2 = N2TMP + A(2,2)**2
      ERROR = N1 * N2
      if (N1 <= THRESH) then
        Q(1,2) = 1.0_rb
        Q(2, 2) = 0.0_rb
        Q(3,2) = 0.0_rb
      else if (N2 <= THRESH) then
        Q(1,2) = 0.0_rb
        Q(2, 2) = 1.0_rb
        Q(3,2) = 0.0_rb
      else if (NORM < (64.0_rb * EPS)**2 * ERROR) then
        T = abs(A(1,2))
        F = -A(1,1) / A(1,2)
        if (abs(A(2,2)) > T) then
          T = abs(A(2,2))
          F = -A(1,2) / A(2,2)
        end if
        if (abs(A(2,3)) > T) then
          F = -A(1,3) / A(2,3)
        end if
        NORM = 1.0_rb / sqrt(1.0_rb + F**2)
        Q(1,2) = NORM
        Q(2,2) = F * NORM
        Q(3,2) = 0.0_rb
      else
        NORM = sqrt(1.0_rb / NORM)
        Q(:,2) = Q(:,2) * NORM
      end if
    else
      A(2,1) = A(1,2)
      A(3,1) = A(1,3)
      A(3,2) = A(2,3)
      A(1,1) = A(1,1) + W(1)
      A(2,2) = A(2,2) + W(1)
      do I = 1, 3
        A(I,I) = A(I,I) - W(2)
        N1 = A(1,I)**2 + A(2,I)**2 + A(3,I)**2
        if (N1 > THRESH) then
          Q(1,2) = Q(2,1) * A(3,I) - Q(3,1) * A(2,I)
          Q(2,2) = Q(3,1) * A(1,I) - Q(1,1) * A(3,I)
          Q(3,2) = Q(1,1) * A(2,I) - Q(2,1) * A(1,I)
          NORM = Q(1,2)**2 + Q(2,2)**2 + Q(3,2)**2
          SUCCESS = NORM <= (256.0_rb * EPS)**2 * N1
          if (.not.SUCCESS) then
            NORM = sqrt(1.0_rb / NORM)
            Q(:, 2) = Q(:, 2) * NORM
            exit
          end if
        end if
      end do
      if (SUCCESS) then
        do J = 1, 3
          if (Q(J,1) <= 0.0_rb) then
            I = 1 + mod(J,3)
            NORM = 1.0_rb / sqrt(Q(J,1)**2 + Q(I,1)**2)
            Q(J,2) = Q(I,1) * NORM
            Q(I,2) = -Q(J,1) * NORM
            Q(I,2) = 0.0_rb
            exit
          end if
        end do
      end if
    end if
    Q(:,3) = cross_product(Q(:,1),Q(:,2))
  end function eigenvectors
  !-------------------------------------------------------------------------------------------------
  pure function cross_product(a, b) result( c )
    real(rb), intent(in) :: a(3), b(3)
    real(rb)             :: c(3)
    c = [ a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1) ]
  end function cross_product
  !-------------------------------------------------------------------------------------------------
  function inv3x3( A ) result( B )
    real(rb), intent(in) :: A(3,3)
    real(rb)             :: B(3,3)
    logical :: OK
    call M33INV( A, B, OK )
    if (.not.OK) call error( "could not invert 3x3 matrix" )
  end function inv3x3
  !-------------------------------------------------------------------------------------------------
  !  M33INV  -  Compute the inverse of a 3x3 matrix.
  !
  !  A       = input 3x3 matrix to be inverted
  !  AINV    = output 3x3 inverse of matrix A
  !  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.

      SUBROUTINE M33INV (A, AINV, OK_FLAG)

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV
      LOGICAL, INTENT(OUT), optional :: OK_FLAG

      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
      DOUBLE PRECISION :: DET
      DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR


      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)

      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         if (present(OK_FLAG)) OK_FLAG = .FALSE.
         RETURN
      END IF

      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      AINV = TRANSPOSE(COFACTOR) / DET

      if (present(OK_FLAG)) OK_FLAG = .TRUE.

      RETURN

      END SUBROUTINE M33INV
  !-------------------------------------------------------------------------------------------------
end module mMath

