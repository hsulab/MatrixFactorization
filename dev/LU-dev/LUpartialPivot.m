#!/usr/local/bin/octave
########################################################################
######               The code is the same for MATLAB.             ######
########################################################################
function [L,U,P] = LUpartialPivot(A)
    % LU-Factorization with Partial Pivoting
    % P*A = LU

    % Init
    n = length(A);

    P = eye(n);

    % Calc
    for k=1:n-1
        % Find the max entry as the pivot
        pivot = max(abs(A(k:n,k)));

        for i=k:n % row
            if abs(A(i,k)) == pivot
                row_ind = i;
                    break
            end
        end

        % Interchange rows and cols
        if k != row_ind
            P([k,row_ind],:) = P([row_ind,k],:);
            A([k,row_ind],:) = A([row_ind,k],:);
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
