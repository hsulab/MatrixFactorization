#!/usr/local/bin/octave
########################################################################
######               The code is the same for MATLAB.             ######
########################################################################
function [L,D,P] = LDLpivot(A)
    % LDL' Decomposition with Block Complete Pivoting
    % P'AP = LDL'
    % SIAM J. Numer. Anal., Vol. 8, 1971, 639-655.

    % A must be a symmetric matrix.
    n = length(A);

    % Block Criterion.
    alpha = (1+sqrt(17))/8;

    mu0 = 0;
    for i=1:n
        for j=1:i
            if abs(A(i,j))>mu0
                mu0=abs(A(i,j));
            end
        end
    end
    
    mu1 = 0;
    for i=1:n
        if abs(A(i,i))>mu1
            mu1=abs(A(i,i));
        end
    end

    % Block Pivoting.
    if mu1>=alpha*mu0
        % Use 1x1 Pivot.
        disp('1x1')
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

    else
        % Use 2x2 Pivot.
        disp('2x2')
        P = eye(n); P1 = eye(n); P2 = eye(n);
        for k=1:n-1
            A_sub = A(k:n,k:n);
            pivot = 0;
            for i=1:n
                for j=1:i
                    if abs(A(i,j))>pivot
                        pivot = abs(A(i,j));
                        row_ind = i;
                        col_ind = j;
                    end
                end
            end

            if (k!=row_ind) && (k+1!=col_ind)
                % Interchange row/col k and row_ind.
                P1([k,row_ind],:) = P1([row_ind,k],:);
                A([k,row_ind],:) = A([row_ind,k],:);
                A(:,[k,row_ind]) = A(:,[row_ind,k]);

                % Interchange row/col k+1 and col_ind.
                P2([k+1,col_ind],:) = P2([col_ind,k+1],:);
                A([k+1,col_ind],:) = A([col_ind,k+1],:);
                A(:,[k+1,col_ind]) = A(:,[col_ind,k+1]);

                %
                P = P2*P1*P;
            end
        end

        %% LDL' Decomposition.
        %E = A(1:2,1:2);
        %E_inv = E^(-1);
        %c = A(3:n,1:2);
        %B = B(3:n,3:n);

        %L = eye(n);
        %L(3:n,1:2) = c*E_inv;

        %D = zeros(n);
        %D(1:2,1:2) = E;
        %D(3:n,3:n) = B-c*E_inv*c';
    end

    % LDL' Decomposition
    %             E | d'      1 |        a11| d'
    % 1x1 pivot ( - - - ) = ( - - - ) * ( - - - )
    %             c | B     c/a11|I         | A2

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
