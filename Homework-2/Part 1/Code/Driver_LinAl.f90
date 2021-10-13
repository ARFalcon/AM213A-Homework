Program Driver_LinAl
  !!!! To complie code run make -f Makefile and to see results open output.txt!!!!
  use LinAl, only: readMat, traceMat, printMat, twoNorm, gEWPP, gEWPP_B, LU, LUB

  implicit none
  character(len=100) :: myFileName
  integer :: i, j, nsizeA, msizeA, msizeB, nsizeB
  real :: traceA, norm, PI, ex, s2
  real, dimension(:), allocatable :: numbers
  integer, dimension(:), allocatable :: s
  real, allocatable :: As(:,:), Bs(:,:), A(:,:), B(:,:), E(:,:), Bss(:,:), Ass(:,:), L(:,:), U(:,:), X(:,:), Plane(:,:)
  logical :: logic

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Problem 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! Reads in matrix A for the homework.
  myFileName = 'Amat.dat'
  call readMat(myFileName,A,nsizeA,msizeA)

  ! Allocate all matricies that are used in problem 2,3,4 
  allocate(numbers(nsizeA))
  allocate(As(msizeA,nsizeA))
  allocate(Ass(msizeA,nsizeA))
  allocate(L(msizeA,nsizeA))
  allocate(U(msizeA,nsizeA))
  As=A
  Ass=A

  ! Reads in martix B for the homework
  myFileName = 'Bmat.dat'
  call readMat(myFileName,B,msizeB,nsizeB)

  ! Allocate matricies that were based off the dimenstions of B
  allocate(Bs(msizeB,nsizeB))
  allocate(Bss(msizeB,nsizeB))
  allocate(E(msizeB,nsizeB))
  Bs=B
  Bss=B

  !Printing the trace of A
  call traceMat(A,msizeA,traceA)
  print *, 'The trace of the matrix is',traceA
  

  !Printing the 2-norm of the columns of the matrix A
  do i = 1,nsizeA


    numbers = A(:,i)
    call twoNorm(numbers,nsizeA,norm)
    print *, 'The norm of column',i,'is',norm
    print *, ''


  enddo
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Problem 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  !Prints A and B before GE with Partial Pivioting takes place
  print*,'Matrix A before Gaussian elimination'
  call printMat(A,msizeA,nsizeA)
  print*,'Matrix B before Gaussian elimination'
  call printMat(B,msizeB,nsizeB)
  
  ! Apply GE with PP
  call gEWPP(A,B,msizeA,nsizeB,logic)

  !Prints A and B after GE with Partial Pivioting takes place 
  print *,'Matrix A after Gaussian elimination'
  call printMat(A,msizeA,nsizeA)
  print *,'Matrix B after Gaussian elimination'
  call printMat(B,msizeB,nsizeB)

  ! Backsubs  A and B and prints A and B
  call gEWPP_B(A,B)
  print *, ''
  print *, 'Soluition Matrix X'
  call printMat(B,msizeB,nsizeB)
  E = matmul(As,B) - Bs
  print *, ''
  print *, 'The error matrix for GE with PP'
  call printMat(E,msizeB,nsizeB)

  ! Prints the 2-Norm of the columns of AsX-Bs
  ! This should be close to machine accuracy.
  print *, ''
  print *, 'The norm of the columns of the error matrix.'
  do i = 1,nsizeB


    numbers = E(:,i)
    call twoNorm(numbers,nsizeB,norm)
    print *, 'The norm of column',i,'is',norm
    print *, ''


  enddo

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Problem 4 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  !Decomposes A into LU form and prints the matirx as well
  call LU(As,msizeA,logic,s)
  call printMat(As,msizeA,nsizeA)


  print *, ''

  ! Preforms the backsubstitution 
  call LUB(As,msizeA,Bs,s)
  call printMat(Bs,msizeB,nsizeB)


  print *, ''

  !Prints L from A
  print *, 'The lower triangular matrix from the LU decomposition'
  print *, ''

  do i = 1, nsizeA
    do j = 1, msizeA

      if ( j < i  ) then
        L(j,i) = 0
      elseif (j == i) then
        L(j,i) = 1.0
      else
        L(j,i) = As(j,i)
      endif

    enddo

  enddo

  call printMat(L,msizeA,nsizeA)


  print *, ''

  !Prints U from A
  print *, 'The Upper triangular matrix from the LU decomposition'
  print *, ''


  do i = 1, nsizeA
    do j = 1, msizeA

      if ( j > i  ) then
        U(j,i) = 0
      elseif (j == i) then
        U(j,i) = As(j,i)
      else
        U(j,i) = As(j,i)
      endif

    enddo
    
  enddo

  call printMat(U,msizeA,nsizeA)


  print *, ''
  print *, 'The Soluition Matrix'
  call printMat(Bs,msizeB,nsizeB)
  ! Prints the 2-Norm of the columns of Ass X- Bs
  ! This should be close to machine accuracy.
  E = matmul(Ass,Bs) - Bss
  do i = 1,nsizeB

    numbers = E(:,i)
    call twoNorm(numbers,nsizeB,norm)
    print *, 'The norm of column',i,'is',norm
    print *, ''


  enddo

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Problem 5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!


  allocate(Plane(3,3))
  allocate(X(3,1))

  PI = 4.D0*DATAN(1.D0)
  ex = exp(1.0)
  s2 = sqrt(2.0)

  ! In this problem we want to solve ax + by + cz + d = 0
  ! We need to get 3 independent eqations so we can solve for (a,b,c) 
  ! The 3 equations we are going to use is 
  ! x1 - x2 + y1 - y2 + z1 - z2 = 0
  ! x1 - x3 + y1 - y3 + z1 - z3 = 0
  !    x1   +    y1   +    z1   = 1

  Plane = transpose(reshape((/4.0,0.0,-2.0,1.0-PI,2.0-ex,3.0+s2,1.0,2.0,3.0/),shape(Plane)))
  X = reshape((/0.0,0.0,1.0/),shape(X))

  call gEWPP(Plane,X,3,1,logic)
  call gEWPP_B(Plane,X)

 !Prints the soluition vector

  call printMat(X,3,1)
deallocate(As,Bs,A,B,E,Bss,Ass,L,U,X,Plane,numbers,s)

End Program Driver_LinAl