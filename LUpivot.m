#!/usr/local/bin/octave
#####################################################
### Octave is a free software similiar to MATLAB. ###
#####################################################
function varargout = LUpivot(varargin)
    % LU-Factorization.

        % No Pivoting.
        % A = LU, [L,U] = LUpivot(A)

        % Partial Pivoting.
        % PA = LU, [L,U,P] = LUpivot(A,'partial')

        % Complete Pivoting.
        % P1*A*P2 = LU, [L,U,P1,P2] = LUpivot(A,'complete')

    % Input Arguments.
    if nargin == 1
        A = varargin{1};
        pivotstrategy = 'no';
    elseif nargin == 2
        A = varargin{1};
        if strcmp(varargin{2},'partial') || strcmp(varargin{2},'complete')
            pivotstrategy = varargin{2};
        else
            error('Wrong Pivoting Strategy.\n');
        end
    else
        error('Wrong Input Arguments.\n');
    end

    % Get the dimension of Matrix A.
    n = length(A);

    % Different Pivot Strategies.
    if strcmp(pivotstrategy,'partial')
        % Set Initial P.
        P = eye(n);

        % Calculate P.
        for k=1:n-1
            % Find the max entry as the pivot.
            pivot = max(abs(A(k:n,k)));

            for i=k:n % row
                if abs(A(i,k)) == pivot
                    row_ind = i;
                        break
                end
            end

            % Interchange rows.
            if k != row_ind
                P([k,row_ind],:) = P([row_ind,k],:);
                A([k,row_ind],:) = A([row_ind,k],:);
            end
        end

    elseif strcmp(pivotstrategy,'complete')
        % Set Initial P1 and P2.
        P1 = eye(n); P2 = eye(n);

        % Calculate P1 and P2.
        for k=1:n-1
            % Find the max element as the pivot.
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

            % Interchange rows and cols.
            if (k != row_ind) || (k != col_ind)
                P1([k,row_ind],:) = P1([row_ind,k],:);
                P2(:,[k,col_ind]) = P2(:,[col_ind,k]);
                A([k,row_ind],:) = A([row_ind,k],:);
                A(:,[k,col_ind]) = A(:,[col_ind,k]);
            end
        end
    end

    % Gaussian Elimination.
    L = eye(n); U = A;
    for k=1:n-1
        for j=k+1:n
            L(j,k) = U(j,k)/U(k,k);
            U(j,k:n) = U(j,k:n) - L(j,k)*U(k,k:n);
        end
    end

    % Output Arguments.
    varargout{1} = L;
    varargout{2} = U;

    if nargin == 2
        if strcmp(pivotstrategy,'partial')
            varargout{3} = P;
        elseif strcmp(pivotstrategy,'complete')
            varargout{3} = P1;
            varargout{4} = P2;
        end
    end

end
