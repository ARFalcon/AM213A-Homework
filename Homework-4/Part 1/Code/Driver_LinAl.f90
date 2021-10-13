Program Driver_LinAl
  use LinAl
  implicit none

  real, allocatable, dimension(:,:) :: A, B, eign, C, EigVec, mu, D, mu2
  real, allocatable, dimension(:) :: eig
  real :: norm
  integer :: n

  allocate(A(4,4))
  allocate(B(3,3))
  allocate(eig(3))
  allocate(eign(3,1))
  allocate(C(4,4))
  allocate(EigVec(4,3))
  allocate(mu(4,1))
  allocate(D(2,2))
  allocate(mu2(2,1))



  !Problem 1
    print *, '########################################################'
    print *, '#                       Problem 1                      #'
    print *, '########################################################'
    A = reshape((/5.0,4.0,1.0,1.0,4.0,5.0,1.0,1.0,1.0,1.0,4.0,2.0,1.0,1.0,2.0,4.0/),shape(A))
    ! Calls the subroutine to tridiagonalize the matrix
    call TriDiag(A)
    ! Prints the tridiagonal matrix
    call printMat(A,4,4)
  
  !Problem 2
    print *, '########################################################'
    print *, '#                       Problem 2                      #'
    print *, '########################################################'
    !Without Shifts
    B = reshape((/3.0,1.0,0.0,1.0,2.0,1.0,0.0,1.0,1.0/),shape(B))
    ! Calls the subroutine to find eigen values without shifts
    call QRNoShiftsEig(B,eig)

    eign(:,1) = eig
    ! Prints the vector of eigenvalues 
    call printMat(eign,3,1)

    !With Shifts and deflation

    B = reshape((/3.0,1.0,0.0,1.0,2.0,1.0,0.0,1.0,1.0/),shape(B))

    ! Calls the subroutine to find eigen values with shifts and deflation
    call QRShiftsEig(B,eig)
  
    eign(:,1) = eig
    ! Prints the vector of eigenvalues 
    call printMat(eign,3,1)

  !Problem 3
    print *, '########################################################'
    print *, '#                       Problem 3                      #'
    print *, '########################################################'
    ! Input the matrix from problem 3 with the associated eigenvalues
    C = reshape((/2.0,1.0,3.0,4.0,1.0,-3.0,1.0,5.0,3.0,1.0,6.0,-2.0,4.0,5.0,-2.0,-1.0/),shape(C))
    mu = reshape((/-8.0286,7.9329,5.6689,-1.57319/),shape(mu))
    
    ! Apply the inverse iteration method
    call InvIter(C,mu,EigVec)
    ! Prints the Eigenvector Matrix 
    call printMat(EigVec,4,4)


End Program Driver_LinAl