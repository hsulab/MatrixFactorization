#!/usr/local/bin/octave
#####################################################
### Octave is a free software similiar to MATLAB. ###
#####################################################
function [Q,H,A_fro] = BasicQH(A)
    % This program computes the polar decomposition A = QH of A R(3x3).
    % This program is modified from algorithm 3.1.

    % Initial Check.
    [m,n] = size(A);
    if (m!=n) || (m!=3) || (n!=3)
        error('Matrix A is not 3x3.\n');
    end
   
    %0. transform input matrix A by A / ||A||F.
    A_fro = norm(A,'fro');
    A = A/A_fro;
    
    %1. Form B R(4x4) in from A.
    B = FormMatrixB(A);
    
    %2. Compute an eigendecomposition B = Z^Z' using the QR algorithm.
    [eigvecs, eigvals] = eig(B);

    %3. Select a dominant eigenvector v.
    maxeigval = 0;
    for i=1:length(eigvals)
        if abs(eigvals(i)) >= maxeigval
            maxeig = eigvals(i);
            maxeigvec = eigvecs(:,i);
        end
    end

    %4. Form the matrix Q from v using (2.7)
    Q = FormMatrixQ(maxeigvec);
    
    %5. Compute the upper triangle of H = Q'*A, 
    % and set the lower triangle equal to the upper triangle.
    H = Q'*A;
    H = triu(H,0) + triu(H,1)';

end
