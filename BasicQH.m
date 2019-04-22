#!/usr/local/bin/octave
function [Q,H] = BasicQH(A)
% BasicQH computes the polar decomposition A = QH of A R(3x3).
% IN:
%   A - 3x3 matrix
% OUT:
%   Q - 3x3 orthogonal matrix
%   H - 3x3 symmetric positive semidefinite matrix
% Usage:
%   [Q,H] = BasicQH(A)

%0. normalize A
[A,A_fro] = ScaleA(A);

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
H = A_fro*Q'*A;
H = triu(H,0) + triu(H,1)';

end
