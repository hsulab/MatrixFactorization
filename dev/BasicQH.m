#!/usr/local/bin/octave
#####################################################
### Octave is a free software similiar to MATLAB. ###
#####################################################
% This program computes the polar decomposition A = QH of A R(3x3).
% This algorithm is modified from algorithm 3.1.

%0. transform input matrix A by A / ||A||F, frobenius norm
% ||A||F = (A'A)^(1/2)
A = [1,2,3; 4,5,6; 7,8,9];
A = A / norm(A, 'fro');

%1. Form B R(4x4) in from A.
B = FormMatrixB(A);

%2. Compute an eigendecomposition B = Z^Z' using the QR algorithm.
% B is symmetric
Bk = B;
for i=1:1000
    [Qk,Rk] = qr(Bk);
    Bk = Rk*Qk;
end

[eigvecs, eigvals] = eig(B);
%3. Select a dominant eigenvector v.
eigvals = diag(Bk,0);

maxeigval = 0;
for i=1:length(eigvals)
    if abs(eigvals(i)) >= maxeigval
        maxeig = eigvals(i);
    end
end

maxeigvec = eigvecs(:,1);
%4. Form the matrix Q from v using (2.7)

v1 = maxeigvec(1); v2 = maxeigvec(2);
v3 = maxeigvec(3); v4 = maxeigvec(4);

Q = [1-2*(v3^2+v4^2), 2*(v2*v3+v1*v4), 2*(v2*v4-v1*v3); ...
        2*(v2*v3-v1*v4), 1-2*(v2^2+v4^2), 2*(v3*v4+v1*v2); ...
        2*(v2*v4+v1*v3), 2*(v3*v4-v1*v2), 1-2*(v2^2+v3^2);];

%5. Computer the upper triangle of H = Q' and set the lower triangle equal to the upper triangle.
H = Q'*A;
H = triu(H,0) + triu(H,1)';
A;
Q*H;
