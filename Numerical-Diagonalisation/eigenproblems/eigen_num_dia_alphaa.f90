! Computes the next-to-leading eigenvalue and eigenvector for the (alpha,a)-family
! of jump distributions.

! Additional functions that it employs:
!   - w_alphaa.f90
!   - U_fun.f90
!   - R_fun_alphaa.f90

subroutine eigen_num_dia_alphaa(lambda,phi,alpha,a,L,N)
  
  use omp_lib                ! Import OpenMP library for parallelization
  implicit none              ! Ensure all variables are explicitly declared

  ! Input parameters
  integer, intent(in) :: N       ! Number of grid points
  real*8, intent(in) :: alpha    ! Parameter for the weight function
  real*8, intent(in) :: a        ! Additional parameter for the weight function
  real*8, intent(in) :: L        ! Half-length of the domain (spanning [-L, L])

  ! Output parameters
  real*8, intent(out) :: lambda  ! Next-to-leading eigenvalue
  double precision, intent(out) :: phi(N)  ! Corresponding eigenvector (normalized)

  ! Local variables
  integer :: i, j                ! Loop counters
  integer :: lwork, info         ! Workspace size and info for LAPACK routine
  real*8 :: x, y                 ! Coordinates in the domain
  real*8 :: dl                   ! Grid spacing
  real*8 :: diffU                ! Difference in potential U(x) - U(y)
  real*8 :: pi                   ! π
  real*8 :: w_alphaa, U_fun      ! Weight function and potential function
  real*8 :: R_fun_alphaa, sum    ! Boundary condition function, normalization sum
  double precision :: wr(N)      ! Array to store eigenvalues
  double precision :: work2(1)   ! Workspace query for LAPACK routine
  double precision, allocatable :: work(:), K(:,:) ! Workspace and kernel matrix

  ! Parameters
  dl = 2.d0 * L / N             ! Compute grid spacing
  pi = 4.d0 * datan(1.d0)       ! Compute π using arctangent
  allocate(K(N, N))             ! Allocate memory for the kernel matrix

  ! Parallelized construction of the kernel matrix K
  !$omp parallel private(i, j, x, y, diffU)
  !$omp do
  do i = 1, N
     x = dl * i - L             ! Compute the x-coordinate for grid point i
     do j = 1, i
        y = dl * j - L          ! Compute the y-coordinate for grid point j
        diffU = U_fun(x) - U_fun(y)  ! Compute potential difference

        ! Populate the kernel matrix with the weight function
        K(i, j) = w_alphaa(x - y, alpha, a, L, N) * exp(-abs(diffU) / 2.d0) * dl
        if (i .eq. j) then
           K(i, j) = K(i, j) + R_fun_alphaa(y, alpha, a)  ! Add diagonal correction
        endif
        K(j, i) = K(i, j)       ! Enforce symmetry of the kernel matrix
     enddo
  enddo
  !$omp enddo
  !$omp end parallel

  ! Perform eigenvalue decomposition of the symmetric matrix K using LAPACK
  call dsyev('V', 'U', N, K, N, wr, work2, -1, info)  ! Query workspace size
  lwork = int(work2(1))        ! Determine optimal workspace size
  allocate(work(lwork))        ! Allocate workspace
  call dsyev('V', 'U', N, K, N, wr, work, lwork, info) ! Compute eigenvalues/vectors

  ! Extract the next-to-leading eigenvalue
  lambda = wr(N-1)

  ! Normalize the eigenvector corresponding to the next-to-leading eigenvalue
  sum = 0.d0
  do i = 1, N
      sum = sum + K(i, N-1)**2 * dl
  enddo
  do i = 1, N
      phi(i) = K(i, N-1) / sqrt(sum)  ! Normalize eigenvector
  enddo

end subroutine eigen_num_dia_alphaa