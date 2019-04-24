function [U,S,V] = SVD(A,tol)
%SVDSIM computes SVD of a matrix.
% IN:
%   A - matrix
%   tol - tolerance
% Usage: 
%   [U,S,V] = SVD(A)
%   or S = SVD(A)
% Notes:
%   with A = U*S*V' , S>=0 , U'*U = Iu  , and V'*V = Iv
%   The idea is to use the QR decomposition on A to gradually "pull" U out from
%   the left and then use QR on A transposed to "pull" V out from the right.
%   This process makes A lower triangular and then upper triangular alternately.
%   Eventually, A becomes both upper and lower triangular at the same time,
%   (i.e. Diagonal) with the singular values on the diagonal.

% check arguments
if ~exist('tol','var')
   tol = eps*1024;
end
   
%reserve space in advance
sizea = size(A);
loopmax = 100*max(sizea);
loopcount = 0;

% or use Bidiag(A) to initialize U, S, and V
U = eye(sizea(1));
S = A';
V = eye(sizea(2));

% iteration
Err = realmax;
while Err > tol && loopcount < loopmax 
    % QR iteration
    [Q,S] = QR(S'); U = U*Q; % flops: 60 + 45
    [Q,S] = QR(S'); V = V*Q;

    % exit when we get "close"
    e = triu(S,1);
    E = norm(e(:));
    F = norm(diag(S));
    if F == 0 
        F = 1;
    end
    Err = E/F;
    loopcount = loopcount + 1;
end
loopcount;

%fix the signs in S
ss = diag(S);
S = zeros(sizea);
for n=1:length(ss)
    ssn = ss(n);
    S(n,n) = abs(ssn);
    if ssn < 0
       U(:,n) = -U(:,n);
    end
end

% check out
if nargout<=1
   U=diag(S);
end

end
