#!/usr/local/bin/octave
########################################################################
######               The code is the same for MATLAB.             ######
########################################################################
function [L,D,P] = LDLdiagonalPivot(A)
    % LDL' Decomposition with Complete Pivoting
    % P'AP = LDL'
    % SIAM J. Numer. Anal., Vol. 8, 1971, 639-655.

    % A must be a symmetric matrix.
    n = length(A);

    % Init
    P = eye(n);
    A_diag = diag(A);
    for k=1:n-1
        pivot = max(abs(A_diag(k:n)));
        for j=k:n
            if abs(A_diag(j)) == pivot
                ind = j;
            end
        end

        if j != k
            P(:,[k,ind]) = P(:,[ind,k]);
            A([k,ind],:) = A([ind,k],:);
            A(:,[k,ind]) = A(:,[ind,k]);
        end
    end

    % LDLt Decomposition
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
