#!/usr/local/bin/octave
function varargout = LU(varargin)
% LU does LU factorization on a matrix.
% IN:
%   A - matrix 
%   pivotstrategy - optional, partial or complete
% OUT:
%   L - lower triangle matrix
%   U - upper triangle matrix
%   P1, P2 - optional, permutation matrix 
% Usage:
%   No Pivot: A = LU, [L,U] = LU(A)
%   Partial Pivot: PA = LU, [L,U,P] = LU(A,'partial')
%   Complete Pivot: P1*A*P2 = LU, [L,U,P1,P2] = LU(A,'complete')

% check input Arguments
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

% get the dimension of Matrix A
n = length(A);

% different pivot strategies.
if strcmp(pivotstrategy,'partial')
    % set Initial p
    P = eye(n);

    % calculate P
    for k=1:n-1
        % find the max entry as the pivot
        pivot = max(abs(A(k:n,k)));

        % find pivot index
        for i=k:n % row
            if abs(A(i,k)) == pivot
                row_ind = i;
                    break
            end
        end

        % interchange rows
        if k != row_ind
            P([k,row_ind],:) = P([row_ind,k],:);
            A([k,row_ind],:) = A([row_ind,k],:);
        end
    end
elseif strcmp(pivotstrategy,'complete')
    % set initial P1 and P2
    P1 = eye(n); P2 = eye(n);

    % calculate P1 and P2
    for k=1:n-1
        % find the max element as the pivot
        pivot = max(max(abs(A(k:n,k:n))));

        % find pivot index
        for i=k:n % row
            for j=k:n % column
                if abs(A(i,j)) == pivot
                    row_ind = i;
                    col_ind = j;
                    break
                end
            end
        end

        % interchange rows and cols
        if (k != row_ind) || (k != col_ind)
            P1([k,row_ind],:) = P1([row_ind,k],:);
            P2(:,[k,col_ind]) = P2(:,[col_ind,k]);
            A([k,row_ind],:) = A([row_ind,k],:);
            A(:,[k,col_ind]) = A(:,[col_ind,k]);
        end
    end
end

% gaussian elimination
L = eye(n); U = A;
for k=1:n-1
    for j=k+1:n
        L(j,k) = U(j,k)/U(k,k);
        U(j,k:n) = U(j,k:n) - L(j,k)*U(k,k:n);
    end
end

% output arguments
varargout{1} = L;
varargout{2} = U;

% return permutation matrix 
if nargin == 2
    if strcmp(pivotstrategy,'partial')
        varargout{3} = P;
    elseif strcmp(pivotstrategy,'complete')
        varargout{3} = P1;
        varargout{4} = P2;
    end
end
end
