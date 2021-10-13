!! /codes/lapack/solve1.f90
!! This routine solves Ax=b using DGESV routine in LAPACK.
!!
!! Recall that the naming convention follows (see http://www.netlib.org/lapack/lug/node26.html)
!!   D:  double precision
!!   GE: general type of matrix
!!   SV: simple driver by factoring A and overwriting B with the solution x
!!
!! See also:
!! a. https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/dgesv.htm
!! b. http://www.netlib.org/clapack/old/double/dgesv.c
!!


program solve1
  implicit none
  integer, parameter :: n = 3
  real, dimension(n) :: x,b
  real, dimension(n,n) :: a
  integer :: i, info, lda, ldb, nrhs
  integer, dimension(n) :: ipiv

  a = reshape((/1., 2., 2., 0., -4., -6., 0., 0., -1./), (/n,n/))
  b = (/3.,-6.,1./)
  x = b

  nrhs = 1
  lda = n
  ldb = n

  call dgesv(n, nrhs, a, lda, ipiv, x, ldb, info)

  print*, 'solving Ax = b'
  do i = 1,n
     print *, i, x(i)
  end do
  
end program solve1
