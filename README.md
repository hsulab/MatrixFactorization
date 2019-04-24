# Fast Polar Decomposition Specialized for 3x3 Real Matrix

## Reference
**[1] An algorithm to compute the polar decomposition of a 3x3 matrix**

## BasicQH.m
    **Algorithm 3.1**

## CompleteQH.m
    **Matrix Generation** 
        *FormMatrixB.m*
        ! Generate 4x4 matrix B from 3x3 matrix A.
        
        *FormMatrixQ.m*
        ! Generate 3x3 matrix Q from eigen vector v.

    **Lambda Calculation**
        *EstimateLambda1.m*
        ! Estimate the dominant eigenvalue lambda1.
            --- *PolyMethod.m*
                ! Use the characteristic polynomial formula.
            --- *NewtonMethod.m*
                ! Use Newton's method.
                --- *HornerPoly.m*
                    ! Use Horner's method to calculate polynomial formula.

    **Eigvec Calculation**
        *CalcBpEigvec.m*
        ! Calculate 2x2 Matrix Bp's eigen vector.

    **Matrix Decomposition**
        *LU.m*
        ! complete or partial pivoting
        ! A=LU, PA=LU, P1AP2=LU

        *LDL.m*
        ! Bunch-Parlett pivoting
        ! A=LDL', P'AP=LDL'

        *QR.m*
        ! Gram-Schmidt or Householder Transformation
        ! A=QR

        *SVD.m*
        ! singular value decomposition
        ! A=USV
## Usage
    ```
    MATLAB or Octave command line
    > A = [0.1,0.2,0.3;0.1,-0.1,0;0.3,0.2,0.1;]; % 3x3 matrix 
    > [Q,H] = Complete(A); % polar decomposition of A with Higham's Method 
    ```
