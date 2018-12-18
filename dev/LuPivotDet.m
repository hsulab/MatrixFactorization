#!/usr/local/bin/octave
########################################################################
######               The code is the same for MATLAB.             ######
########################################################################
function [A_det] = LuPivotDet(A)
    % LU-Factorization with Partial Pivoting
    % PA = LU

    % Init
    n = length(A);
    L=eye(n); P=eye(n); U=A;
    n_trans = 0;
    % Calc
    for k=1:n-1 % row
        % Find the max entry as the pivot
        pivot = max(abs(U(k:n,k)));
        for j=k:n % row
            if abs(U(j,k)) == pivot
                ind = j;
                break
            end
        end
        % Interchange rows
        if ind != k
            U([k,ind],k:n) = U([ind,k],k:n);
            L([k,ind],1:k-1) = L([ind,k],1:k-1);
            P([k,ind],:) = P([ind,k],:);
            n_trans = n_trans + 1;
        end
        % Gaussian Elimination
        for j=k+1:n
            L(j,k) = U(j,k)/U(k,k);
            U(j,k:n) = U(j,k:n) - L(j,k)*U(k,k:n);
        end
    end
    % Calc Det
    L_det = 1;
    for i=1:length(diag(L))
        L_det *= L(i,i);
    end
    U_det = 1;
    for i=1:length(diag(U))
        U_det *= U(i,i);
    end
    A_det = (L_det*U_det)/((-1)^n_trans);
    if abs(A_det) < eps
        A_det = 0;
    end
end
