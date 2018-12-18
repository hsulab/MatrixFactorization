# Fast Polar Decomposition Specialized for 3x3 Real Matrix

##Reference
[1] An algorithm to compute the polar decomposition of a 3x3 matrix

   / - Algorithm 3.1
- |                     / - Algorithm 3.2 - Algorithm 3.3 - Algorithm 3.4
   \ - Algorithm 3.5 - |
                        \ - ...

###BasicQH.m
    ** Algorithm 3.1 **

###FastQH.m
    ** Matrix Generation ** 
        *FormMatrixB.m*
        ! Generate 4x4 matrix B from 3x3 matrix A.
        
        *FormMatrixQ.m*
        ! Generate 3x3 matrix Q from eigen vector v.

    ** Lambda Calculation **
        *EstimateLambda1.m*
        ! Estimate the dominant eigenvalue lambda1.
            --- *PolyMethod.m*
                ! Use the characteristic polynomial formula.
            --- *NewtonMethod.m*
                ! Use Newton's method.
                --- *HornerPoly.m*
                    ! Use Horner's method to calculate polynomial formula.

    ** Eigvec Calculation **
        *CalcBpEigvec.m*
        ! Calculate 2x2 Matrix Bp's eigen vector.

    ** Matrix Decomposition **
        *LUpivot.m*
        ! A=LU, PA=LU, P1AP2=LU

        *LDLpivot.m*
        ! A=LDL', P'AP=LDL'

        *QRfactorization.m*
        ! A=QR
