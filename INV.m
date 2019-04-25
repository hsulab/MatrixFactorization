#!/usr/local/bin/octave
function A_inv = INV(A)
% INV computes the inverse of a triangle-like square matrix.
% IN:
%   A - triangle square matrix
% OUT:
%   A_inv - the inverse of A
% Usage:
%   A_inv = INV(A);
% Notes:
%   The inverse of a triangle square matrix equals a diagonal matrix,
%   whose each diagonal entry equals the reciprocal of the entry in 
%   the same place of the matrix.

% check square
[m,n] = size(A);
if m ~= n
    error('Invalid Input Matrix: Not Square.')
end

% init inv
A_inv = zeros(m,n);

% calc diagonal
for i=1:n
    A_inv(i,i) = 1.0/A(i,i);
end

% check upper or lower
t = 0;
if norm(triu(A,1)) > eps
    % upper 
    t = 1;
else 
    if norm(tril(A,-1)) > eps 
        % lower
        A = A';
        t = 2;
    else 
        % diagonal 
        t = 3;
    end
end

% calc other entries 
if t ~= 3
    for i=n-1:-1:1
        for j=i+1:n
            temp = 0;
            for k=i+1:j
                temp = temp - A_inv(k,j)*A(i,k);
            end
            A_inv(i,j) = temp/A_inv(i,i);
        end
    end
end

% check lower
if t == 2 
    A_inv = A_inv';
end

end
