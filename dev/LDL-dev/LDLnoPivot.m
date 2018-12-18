#!/usr/local/bin/octave
########################################################################
######               The code is the same for MATLAB.             ######
########################################################################
function [L,D] = LDLnoPivot(A)
    % LDL' Decomposition without Pivoting

    % A must be a symmetric matrix.
    n = length(A);

    % LDL' Decomposition
    L = tril(A,-1)+eye(n);
    D = diag(A);

    L(2:end,1) = L(2:end,1)/D(1);
    for k=2:n
        v = L(k,1:k-1).*D(1:k-1)';
        D(k) = D(k) - L(k,1:k-1)*v';
        L(k+1:end,k) = (L(k+1:end,k)-L(k+1:end,1:k-1)*v')/D(k);
    end

    D = diag(D);

end
