Program Driver_LinAl
    use LinAl

    Implicit none
    real, allocatable :: A(:,:), B(:,:), R(:,:), Q(:,:), Y(:,:), G(:,:), S(:,:), E(:,:), FitVals(:,:)
    real, allocatable :: Iden(:,:), QB(:,:), SQR(:,:), SQB(:,:), Err(:,:)
    logical :: logic
    character(len=100) :: filename
    integer :: msize, nsize, input, j, i, GramSize, Qsize
    real, allocatable :: error(:)
    real :: norm, QRNorm, QQNorm

    !Reads in the data file
    filename = 'atkinson.dat'
    call readMat(filename,A,msize,nsize)

!~~~~~~~~~~~~~~~~~~~~~~Problem 1~~~~~~~~~~~~~~~~~~~~~~~~!

  

    !For a 3rd Degree Polynomial

    input = 3
    allocate(Y(msize,1))
    allocate(G(input+1,input+1))
    allocate(S(input+1,1))


    !Creats Vandermonde Matrix
    call VanMat(A,B,input)

    Y(:,1) = A(:,2)

    
    !Creats the normal equation A^TAx=A^Tb
    S = matmul(transpose(B),Y)
    G = matmul(transpose(B),B)
    GramSize = size(G,1)

    !Solves the normal equation Using Cholesky decomposition and backsubstitution.
    call CD(G,logic)
    call CDB(G,S)
  
    !Prints the solution vector
    call printMat(S,GramSize,1)


    !Calculates the 2-norm error between the fitted curve and the data
    allocate(E(msize,1))
    allocate(FitVals(msize,1))
    FitVals = matmul(B,S)
    allocate(error(msize))

    E = FitVals-Y
    error = E(:,1)

    !Prints the fitted curve to a file
    call matToFile(FitVals,'solnVecP3.dat')
    call twoNorm(error,msize,norm)
    print *, 'The norm of error vector is',norm
    print *, ''


    deallocate(error,E,FitVals,S,G,Y,B)

    !For a 5th Degree Polynomial

    input = 5
    allocate(Y(msize,1))
    allocate(G(input+1,input+1))
    allocate(S(input+1,1))



    !Creats Vandermonde Matrix
    call VanMat(A,B,input)

    Y(:,1) = A(:,2)

    !Creats the normal equation A^TAx=A^Tb
    S = matmul(transpose(B),Y)
    G = matmul(transpose(B),B)
    GramSize = size(G,1)


    !Solves the normal equation Using Cholesky decomposition and backsubstitution.
    call CD(G,logic)
    call CDB(G,S)


    !Prints the solution vector
    call printMat(S,GramSize,1)


    !Calculates the 2-norm error between the fitted curve and the data
    allocate(E(msize,1))
    allocate(FitVals(msize,1))
    FitVals = matmul(B,S)
    allocate(error(msize))

    E = FitVals-Y
    error = E(:,1)


    !Prints the fitted curve to a file
    call matToFile(FitVals,'solnVecP5.dat')
    call twoNorm(error,msize,norm)
    print *, 'The norm of error vector is',norm
    print *, ''

    deallocate(error,E,FitVals,S,G,Y,B)



!~~~~~~~~~~~~~~~~~~~~~~Problem 2~~~~~~~~~~~~~~~~~~~~~~~~!


    !For a 3rd Degree Polynomial
    input = 3

    !Creating vandermonde matrix
    allocate(Y(msize,1))
    call VanMat(A,B,input)

    Y(:,1) = A(:,2)    

    GramSize = size(B,2)

    !Does a QR decomposition on the vandermonde matrix
    call QR(B,Q,R)

    allocate(E(msize,GramSize))

    call printMat(R,msize,GramSize)

    !Computes A-QR and prints it to the screen
    E = B - matmul(Q,R)

    call printMat(E,msize,GramSize)

    !Computes the Forbenius norm and prints it to the screen
    call ForNormMat(E,QRNorm)

    print *, 'The norm of A-QR matrix is',QRNorm
    print *, ''

    !Creates an Idenity matrix 
    Qsize = size(Q,1)
    allocate(Iden(Qsize,Qsize))
    Iden = 0.0
    
    do i = 1,Qsize

      Iden(i,i) = 1

    enddo

    allocate(S(Qsize,Qsize))

    !Computes Q^TQ -I
    S = matmul(transpose(Q),Q) - Iden

    call printMat(S,Qsize,Qsize)

    !Prints norm to screen
    call ForNormMat(S,QQNorm)
    print *, 'The norm of Q^TQ-I matrix is',QQNorm
    print *, ''


    !Solves the least square eq
    allocate(QB(Qsize,1))
    allocate(SQR(GramSize,GramSize))
    allocate(SQB(GramSize,1))
    QB = matmul(transpose(Q),Y)

    SQB(:,1) = QB(1:GramSize,1)

    SQR = R(1:GramSize,1:GramSize)

    call gEWPP_B(SQR,SQB)

    !Prints the solution vector
    call printMat(SQB,GramSize,1)

    allocate(Err(Qsize,1))

    !Calculates the 2-norm of Fitted data and Original data
    Err = matmul(B,SQB) - Y

    call ForNormMat(Err,norm)
    print *, 'The norm of A*x-b vector is',norm
    print *, ''

    !Writes the fitted data to a file
    allocate(FitVals(msize,1))
    FitVals = matmul(B,SQB)
    call matToFile(FitVals,'solnVecP3QR.dat')

    deallocate(Err,SQB,SQR,QB,S,Iden,E,Y,FitVals)

    !For a 5th Degree Polynomial

    input = 5

    !Creating vandermonde matrix
    allocate(Y(msize,1))
    call VanMat(A,B,input)

    Y(:,1) = A(:,2)    

    GramSize = size(B,2)

    !Does a QR decomposition on the vandermonde matrix
    call QR(B,Q,R)

    allocate(E(msize,GramSize))

    call printMat(R,msize,GramSize)

    !Computes A-QR and prints it to the screen
    E = B - matmul(Q,R)

    call printMat(E,msize,GramSize)

    !Computes the Forbenius norm and prints it to the screen
    call ForNormMat(E,QRNorm)

    print *, 'The norm of A-QR matrix is',QRNorm
    print *, ''

    !Creates an Idenity matrix 
    Qsize = size(Q,1)
    allocate(Iden(Qsize,Qsize))
    Iden = 0.0
    
    do i = 1,Qsize

      Iden(i,i) = 1

    enddo

    allocate(S(Qsize,Qsize))

    !Computes Q^TQ -I
    S = matmul(transpose(Q),Q) - Iden

    call printMat(S,Qsize,Qsize)

    !Prints norm to screen
    call ForNormMat(S,QQNorm)
    print *, 'The norm of Q^TQ-I matrix is',QQNorm
    print *, ''


    !Solves the least square eq
    allocate(QB(Qsize,1))
    allocate(SQR(GramSize,GramSize))
    allocate(SQB(GramSize,1))
    QB = matmul(transpose(Q),Y)

    SQB(:,1) = QB(1:GramSize,1)

    SQR = R(1:GramSize,1:GramSize)

    call gEWPP_B(SQR,SQB)

    !Prints the solution vector
    call printMat(SQB,GramSize,1)

    allocate(Err(Qsize,1))

    !Calculates the 2-norm of Fitted data and Original data
    Err = matmul(B,SQB) - Y

    !Writes the fitted data to a file
    call ForNormMat(Err,norm)
    print *, 'The norm of A*x-b vector is',norm
    print *, ''
    allocate(FitVals(msize,1))
    FitVals = matmul(B,SQB)
    call matToFile(FitVals,'solnVecP5QR.dat')

    deallocate(Err,SQB,SQR,QB,S,Iden,E,Y,FitVals)

End Program Driver_LinAl