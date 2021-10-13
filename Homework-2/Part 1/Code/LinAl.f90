module LinAl

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Read in Matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  subroutine readMat(filename,A,msize,nsize)

    implicit none
    character(len=*) :: filename
    real,dimension(:,:),intent(out),allocatable :: A
    integer :: i,j
    integer, intent(out) :: msize, nsize

    ! Reads a file containing the matrix A 
    ! Sample file:
    !
    ! 4 4 
    ! 2.0 1.0 1.0 0.0
    ! 4.0 3.0 3.0 1.0
    ! 8.0 7.0 9.0 5.0
    ! 6.0 7.0 9.0 8.0
    !
    ! Note that the first 2 numbers in the first line are the matrix dimensions, i.e., 4x4,
    ! then the next msize lines are the matrix entries. This matrix is found in Eq. 2.18 of the lecture note.
    ! Note that entries must be separated by a tab.


    open(10,file=filename)

    ! Read the matrix dimensions
    read(10,*) msize,nsize
    allocate(A(msize,nsize))
    ! Read matrix
    do i=1,msize
       read(10,*) ( A(i,j), j=1,nsize )
    enddo

    close(10)
    
  end subroutine readMat

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Trace of Matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  subroutine traceMat(mat,dim,trace)

    ! Takes in a square matrix and its dimension
    ! Outputs the trace of the matrix


    implicit none
    real, dimension(:,:), allocatable, intent(in) :: mat
    integer, parameter :: dp = kind(1.d0)
    integer, intent(in) :: dim
    real(kind=dp), intent(out) :: trace
    integer :: i
    trace = 0
    do i=1,dim
      trace = trace + mat(i,i)
    enddo

  end subroutine traceMat

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 2-Norm of a vector ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!  

  subroutine twoNorm(vector,dim,norm)

    ! Takes in a vector with its lengths
    ! Outputs the 2-norm

    implicit none
    real, dimension(:), allocatable, intent(in) :: vector
    integer, parameter :: dp = kind(1.d0)
    integer, intent(in) :: dim
    real(kind=dp), intent(out) :: norm
    integer :: i

    norm = 0
    do i=1,dim
      norm = norm + vector(i)**2
    enddo
    norm = sqrt(norm)

  end subroutine twoNorm

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Print a Matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  subroutine printMat(mat,dimM,dimN)

    ! Takes in a m x n matrix with its dimension
    ! Writes the matrix in human readable format

    implicit none
    real, dimension(:,:), allocatable, intent(in) :: mat
    integer, intent(in) :: dimM, dimN
    integer :: i
    print *, "The size of the matrix is",  dimM,'x',dimN
    do i = 1, dimM
      print *, mat(i,:)
    enddo



  end subroutine printMat

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Gaussian Elimination with Partial Pivioting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  subroutine gEWPP(A,B,mA,nB,logic)
    implicit none
    real, dimension(:,:), allocatable, intent(inout) :: A, B
    real, allocatable :: u(:,:) , v(:)
    logical, intent(out) :: logic
    integer, intent(in) :: mA, nB
    integer :: i,j,k
    real :: factor

    !Takes in 2 matrics to do the forward steps to solve the equation AX=B using Gaussian elimination
    !with partial pivioting with the dimension of A and the number of columns of B
    !with a logic variable to determine if the matrix is singular or not.
    !Note: If the B is a single vector B must be defined as a m x 1 matrix.

    ! allocate tempoary matrix and vector for swapping
    allocate(u(mA,mA+nB))
    allocate(v(mA+nB))
    !Define agmented matrix
    u(:,1:mA) = A
    u(:,mA+1:mA+nB) = B

    !Start the algorithm
    do j = 1,mA-1
      
      k = maxloc(abs(u(j:mA,j)),1) + j-1 !Determine location of largest element

      if ( k .NE. j ) then

        !Sawp the rows if necessary
        v = u(k,:)
        u(k,:) = u(j,:)
        u(j,:) = v

      endif

      !Tell us if the matrix is singular or not
      print *,''

      if (u(j,j) == 0) then
        logic = .TRUE.
        stop "Matrix is singular"
      endif

      !This gets the matrix in REF form
      do i = j+1, mA 

        factor = u(i,j)/u(j,j)
        u(i,:) = u(i,:) - factor*u(j,:)

      enddo      

    enddo

    !Extracts A and B from the agumented matrix
    B = u(:,mA+1:mA+nB)
    A = u(:,1:mA)
    logic = .False.
    
  end subroutine gEWPP

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Gaussian Elimination with Partial Pivioting Backsubstitution ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  subroutine gEWPP_B(A,B)
    implicit none
    real, allocatable, intent(inout) :: A(:,:)
    real, allocatable, intent(inout) :: B(:,:)
    real, allocatable :: u(:,:)
    integer :: i, j, mA, nB, k
    real :: factor
    mA = size(B,1)
    nB = size(B,2)
    allocate(u(mA,mA+nB))

    !Takes in A and B from the output of Gaussian Elimination with Partial Pivioting (gEWPP)
    !And returns A and B after back substiution where B is the solution matrix

    u(:,1:mA) = A
    u(:,mA+1:mA+nB) = B
    
     do j = mA, 2,-1            ! j is column
        do i = j-1, 1, -1        ! i is row

          if (U(i,i)==0) then
            stop "Matrix is singular"
          endif
          
          factor = u(i,j)/u(j,j)
          u(i,:) = u(i,:) - factor*u(j,:)
        
        enddo
      
      enddo

      do i=1,mA
        do j=mA+1,mA+nB
            u(i,j)=u(i,j)/u(i,i)
        enddo

      enddo

    B = u(:,mA+1:mA+nB)
    A = u(:,1:mA)
  end subroutine gEWPP_B

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LU decomposition with Partial Pivioting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  subroutine LU(A,m,logic,s)
    implicit none
    real, allocatable, intent(inout) :: A(:,:)
    integer, allocatable, intent(out) :: s(:)
    integer, intent(in) :: m
    logical, intent(out) :: logic
    real, allocatable :: v(:)
    integer :: i, j, k, temp, l
    allocate(v(m))
    allocate(s(m))

  ! Solves the eqaution AX=B by LU - decomposition
  ! this subroutine decomposes A as a single LU matrix
  ! Input is A square matrix, dimension of A, logic variable to determine if A is singular or not
  ! and s the permutation vector. 

    do j = 1,m
      s(j) = j
    enddo 

    do j = 1,m
    l = maxloc(abs(A(j:m,j)),1) + j-1  
      if ( l .NE. j ) then
        !Sawp the rows if necessary
        v = A(l,:)
        A(l,:) = A(j,:)
        A(j,:) = v

        temp = s(l)
        s(l) = s(j)
        s(j) = temp
      endif

      if (A(j,j) == 0) then
        logic = .TRUE.
        stop "Matrix is singular"
      endif

      do i = j+1,m
        A(i,j) = A(i,j)/A(j,j)

        do k=j+1,m
          A(i,k)=A(i,k)-(A(i,j)*A(j,k))
        enddo

      enddo

      logic = .FALSE.

    enddo

  end subroutine LU

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LU decomposition with Partial Pivioting Backsubstitution ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  subroutine LUB(A,m,B,s)
    real, allocatable, intent(inout) :: B(:,:)
    real, allocatable, intent(in) :: A(:,:)
    real, allocatable :: Y(:,:) , sum(:)
    integer, allocatable, intent(inout) :: s(:)
    integer, intent(in) :: m
    integer :: i, j, k, l, w, x

    !Finishes the solutition of AX=B by doing LU backsubstitution.
    !Input is A in LU form, m is dimenstion of A, B is matrix you want to solve and s is permutation vector.

    w = size(B,1)
    x = size(B,2)
    allocate(Y(w,x))
    allocate(sum(x))
    
    do j = 1,m
    
      l = s(j)
      Y(j,:) = B(l,:)
    
    enddo
    
    do j = 1,m-1

      do i = j+1,m

        Y(i,:) = Y(i,:) - Y(j,:)*A(i,j)

      enddo

    enddo

    do i = m,1,-1
      if ( A(i,i) == 0 ) then

        stop "Matrix is singular"

      endif

      sum = 0.0

      do k = i+1,m

        sum = sum + A(i,k)*B(k,:)

      enddo

      B(i,:) = (Y(i,:)-sum)/A(i,i)
    enddo

  end subroutine LUB




end module LinAl