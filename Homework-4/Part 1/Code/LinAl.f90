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
    integer, intent(in) :: dim
    real, intent(out) :: trace
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
    integer, intent(in) :: dim
    real, intent(out) :: norm
    integer :: i, j

    

    norm = 0
    do i=1,dim
      norm = norm + vector(i)**2
    enddo
    norm = sqrt(norm)

  end subroutine twoNorm

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 2-Norm of a Mat ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!  

  subroutine twoNormMat(Mat,norm)

    ! Takes in a matrix 
    ! Outputs the 2-norm element wise

    implicit none
    real, dimension(:,:), allocatable, intent(in) :: Mat
    real, intent(out) :: norm
    integer :: i, j, m, n

    m = size(Mat,1)
    n = size(Mat,2)

    norm = 0.0
    do i=1,m
      do j = 1, n
        norm = norm + mat(i,j)**2
      enddo
    enddo
    norm = sqrt(norm)

  end subroutine twoNormMat

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
    implicit none
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

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Cholesky Decomposition ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  subroutine CD(A,logic)

    ! Takes in matrix A
    ! Outputs matrix A in cholesky decomposition form with a PD and singualr logical.


    implicit none
    real, allocatable, intent(inout) :: A(:,:)
    logical :: logic
    integer :: n, m, i, j, k
    m = size(A,1)
    do j = 1, m

      do k = 1, j-1
        A(j,j) = A(j,j) - A(j,k)*A(j,k)
        
      enddo

      if ( A(j,j) < 1.0e-15 ) then 
        
        logic = .False.
        stop "Matrix is not positive definite or matrix is singular."

      else

        A(j,j) = sqrt(A(j,j))
      
      endif

      do i = j+1, m

        do k = 1, j-1

          A(i,j) = A(i,j) - A(i,k)*A(j,k)

        enddo

        A(i,j) = A(i,j)/A(j,j)

      enddo
    enddo
    logic = .TRUE.

  end subroutine CD

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Cholesky Decomposition Backsumstitution ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  subroutine CDB(L,B)

  ! Takes in lower triaglular matrix from Cholesky decomposition with B matrix from AX=B
  ! Returns B matrix after applying back substitution to it.


  implicit none
  real, allocatable, intent(inout) :: L(:,:), B(:,:)
  integer :: i, j, k, n, m
  real, allocatable :: sum(:)
  n = size(B,2)
  m = size(L,1)
  allocate(sum(n))

  sum = 0.0
  
  !foward substitution
  do i = 1, m
    sum = B(i,:)

    do j = 1, i-1

      sum = sum - B(j,:)*L(i,j)

    enddo

    B(i,:) = sum/L(i,i)
  
  enddo

  !backwards substitution
  do i = m, 1, -1
    if ( L(i,i) < 1.0e-16 ) then

      stop "Matrix is singular!"

    endif

    do k = i+1 , m
    
      B(i,:) = B(i,:) - L(k,i)*B(k,:)

    enddo

    B(i,:) = B(i,:)/L(i,i)

  enddo


  end subroutine CDB

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ QR Factorization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  subroutine QR(A,Q,R)

    ! Takes in matrix A and does a full QR decomposition 
    ! Returns Q and R.

    real, allocatable, intent(in) :: A(:,:)
    real, allocatable, intent(out) :: Q(:,:), R(:,:)
    real, allocatable :: s(:), v(:,:), temp(:,:), Iden(:,:), Z(:,:)
    integer :: i, j, k, l, m, n
    real :: norm        

    m = size(A,1)
    n = size(A,2)

    allocate(s(n))
    allocate(v(m,n))
    allocate(Iden(m,m))
    allocate(Q(m,m))
    allocate(R(m,n))
    allocate(Z(m,m))

    ! Creat an Identity matrix for Q = I - 2vv^T
    Iden = 0.0
    do i = 1,m
      Iden(i,i) = 1
    enddo

    v = 0.0
    Q = Iden
    R = A

    !Start the Loop

    do j = 1, n

      norm = 0.0

      !Calculate the norm of A(i,:) from i to n
      do i = j, m

        norm = norm + R(i,j)**2

      enddo

      norm = sqrt(norm)


      s(j) = sign(norm,R(j,j))
      !Creats the vector that is used to creat a householder matrix
      v(j,j) = R(j,j) + s(j)

      do i = j+1, m

        v(i,j) = R(i,j)

      enddo

      norm = 0.0

      do i=j, m

        norm = norm + v(i,j)**2

      enddo

      norm = sqrt(norm)
      
      v(:,j) = v(:,j)/norm

      !Create an outer product matrix for R = R - 2*(vv^T)*R
      ! and for Q = I - 2vv^T

      do k = 1, m
        do l = 1, m

          Z(l,k) =  v(l,j)*v(k,j)

        enddo
      enddo

      ! the creation of Q and R

      R = R - 2.0*matmul(Z,R)
      Q = matmul(Q,Iden - 2*Z)

    enddo
    
    R = -R
    Q = -Q

  end subroutine QR

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Vandermonde Matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  subroutine VanMat(A,Van,P)

    ! Takes in a vector A in Matrix form i.e. (mx1) and what power of polynomial P
    ! output the vandermonde matrix

    implicit none
    real, allocatable, intent(in) :: A(:,:)
    real, allocatable, intent(out) :: Van(:,:)
    integer, intent(in) :: P
    integer :: i, j, msize
    
    msize = size(A,1)
    allocate(Van(msize,p+1))


    do i = 1, msize
      do j = 1, p+1

        van(i,j) = A(i,1)**(j-1)

      enddo

    enddo
    
  end subroutine VanMat

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Mat to File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  subroutine matToFile(Mat,filename)

    !Takes in a matrix and a filename
    !Writes the matrix to the file

    implicit none
    real, allocatable, intent(in) :: Mat(:,:)
    character(len=*),intent(in) :: filename
    integer :: i, msize, nsize, j

    msize = size(Mat,1)
    nsize = size(Mat,2)

    open(1, file = filename)

    do i = 1, msize

      write(1,*) (Mat(i,j) , j = 1, nsize)

    enddo

    close(1)

  end subroutine matToFile
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fourbenius-Norm Matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  subroutine ForNormMat(mat,norm)

    ! Takes in a matrix
    ! Outputs the 2-norm

    implicit none
    real, dimension(:,:), allocatable, intent(in) :: mat
    real, intent(out) :: norm
    integer :: i, j, n, m

    n = size(Mat,2)
    m = size(Mat,1)

    norm = 0
    do i = 1, m
      do j = 1, n

        norm = norm + mat(i,j)**2

      enddo

    enddo
    
    norm = sqrt(norm)

  end subroutine ForNormMat

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Tridiagonalize Matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  subroutine TriDiag(A)
    
    !Takes in a Matrix A
    !Outputs a Tridiagonal matrix using householder matrices
    
    implicit none
    real, allocatable, dimension(:,:), intent(inout) :: A
    integer :: m, n, i, j, k, l
    real :: s, norm
    real, allocatable :: v(:), Z(:,:)


    m = size(A,1)
    
    allocate(v(m))
    allocate(Z(m,m))

    do j = 1, m-1
      
      v = 0.0
      ! Computes the norm for the Signed norm in the next step.
      norm = 0
      do i = j+1, m

        norm = norm + A(i,j)**2

      enddo
      norm = sqrt(norm)
      !Computes the signed norm
      s = sign(norm,A(j+1,j))
      ! Extracts the vector from A for the application of A - 2*vv^TA and A - 2*Avv^T
      v(j+1) = A(j+1,j) + s
      do i = j+2, m

        v(i) = A(i,j)

      enddo
      !Finds the norm of v 
      norm = 0
      do i = j+1, m

        norm = norm + v(i)**2

      enddo
      norm = sqrt(norm)
      ! Normalizes v
      v = v/norm
      ! Calculate the outerproduct vv^T
      do k = 1, m
        do l = 1, m

          Z(l,k) =  v(l)*v(k)

        enddo
      enddo
      
      ! Applys the householder transformation to make A tridiagonal
      A = A - 2*matmul(Z,A)
      A = A - 2*matmul(A,Z)

    enddo
  end subroutine TriDiag

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Elementwise 2-Norm of a the Sub-diagonal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  subroutine SubDiagNorm(A,norm)
    implicit none
    real, allocatable, dimension(:,:), intent(in) :: A
    integer :: i, m
    real, intent(out):: norm
    m = size(A,1)
    norm = 0.0

    do i = 1, m-1

      norm = norm + A(i+1,i)**2

    enddo

    norm = sqrt(norm)
    
  end subroutine SubDiagNorm

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ QR No Shifts Eigenvalues ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  subroutine QRNoShiftsEig(A,eig)
    real, allocatable, dimension(:,:), intent(inout) :: A
    real, allocatable, dimension(:), intent(out) :: eig
    real :: norm, tol
    real, allocatable, dimension(:,:) :: Q, R
    integer :: i, m


    m = size(A,1)
    tol = 10.0**(-14)

    allocate(eig(m))

    call SubDiagNorm(A,norm)

    do while (norm > tol)

      call QR(A,Q,R)
      A = matmul(R,Q)

      call SubDiagNorm(A,norm)

    end do

    do i = 1, m
      eig(i) = A(i,i)
    enddo

    
  end subroutine QRNoShiftsEig

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ QR Shifts Eigenvalues ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  subroutine QRShiftsEig(A,eig)
    implicit none
    real, allocatable, dimension(:,:), intent(inout) :: A
    real, allocatable, dimension(:), intent(out) :: eig
    real :: norm, tol, u
    real, allocatable, dimension(:,:) :: Q, R, Iden, B
    integer :: i, j, m, k


    m = size(A,1)
    tol = 10.0**(-4)
    allocate(eig(m))

    
    do j = m, 2, -1

      allocate(B(j,j))
      B = A(1:j,1:j)
  
      allocate(Iden(j,j))
      Iden = 0.0
      do i = 1,j
        Iden(i,i) = 1.0
      enddo
      
      norm = 1.0

      do while (norm > tol)

        u = B(j,j)
        B = B - u*Iden
        call QR(B,Q,R)
        B = matmul(R,Q) + u*Iden
        call SubDiagNorm(B,norm)

      end do

      A(1:j,1:j) = B
      deallocate(B,Iden)

    enddo

    

    do i = 1, m
      eig(i) = A(i,i)
    enddo
    
  end subroutine QRShiftsEig

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Inverse Iteration Algorithm ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  
  subroutine InvIter(A,mu,EigVec)
    implicit none
    real, allocatable, dimension(:,:), intent(inout) :: A
    real, allocatable, dimension(:,:), intent(in) :: mu
    real, allocatable, dimension(:,:), intent(out) :: EigVec
    real, allocatable, dimension(:,:) :: Iden, B, Y, R, Q, X, ASave, XSave, S
    real :: norm, u, tol, norm_R, norm_S
    integer :: i, j, k, m, n
    logical :: logic

    m = size(A,1)
    n = size(mu,1)

    allocate(Iden(m,m))
    allocate(B(m,m))
    allocate(Y(m,1))
    allocate(R(m,1))
    allocate(X(m,1))
    allocate(EigVec(m,n))
    allocate(ASave(m,m))
    allocate(XSave(m,1))
    allocate(S(m,1))

    tol = 10.0**(-4)

    ASave = A
    Iden = 0.0
    do i = 1,m
        Iden(i,i) = 1.0
    enddo

    do i = 1, n

      u = mu(i,1)

      A = ASave - u*Iden

      norm_R = 1.0
      norm_S = 1.0
      X = 0.0
      X(1,1) = 1.0
      
      do while (norm_R > tol .and. norm_S > tol)
        B = A
        XSave = X
        call gEWPP(B,X,m,1,logic)
        call gEWPP_B(B,X) 
        call twoNormMat(X,norm)

        X = X/norm

        R = X - XSave
        S = X + XSave


        call twoNormMat(R,norm_R)
        call twoNormMat(S,norm_S) 

      enddo

      EigVec(:,i) = X(:,1)
    
    enddo
    
  end subroutine InvIter  



end module LinAl