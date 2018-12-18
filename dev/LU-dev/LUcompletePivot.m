#!/usr/local/bin/octave
########################################################################
######               The code is the same for MATLAB.             ######
########################################################################
function [L,U,P1,P2] = LUcompletePivot(A)
    % LU-Factorization with Complete Pivoting
    % P1*A*P2 = LU

    % Init
    n = length(A);

    P1 = eye(n); P2 = eye(n);

    % Calc
    for k=1:n-1
        % Find the max element as the pivot
        pivot = max(max(abs(A(k:n,k:n))));

        for i=k:n % row
            for j=k:n % column
                if abs(A(i,j)) == pivot
                    row_ind = i;
                    col_ind = j;
                    break
                end
            end
        end

        % Interchange rows and cols
        if (k != row_ind) || (k != col_ind)
            P1([k,row_ind],:) = P1([row_ind,k],:);
            P2(:,[k,col_ind]) = P2(:,[col_ind,k]);
            A([k,row_ind],:) = A([row_ind,k],:);
            A(:,[k,col_ind]) = A(:,[col_ind,k]);
        end
    end

    % Gaussian Elimination
    L = eye(n); U = A;
    for k=1:n-1
        for j=k+1:n
            L(j,k) = U(j,k)/U(k,k);
            U(j,k:n) = U(j,k:n) - L(j,k)*U(k,k:n);
        end
    end

end
