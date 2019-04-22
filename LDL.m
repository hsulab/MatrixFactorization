#!/usr/local/bin/octave
function varargout = LDL(varargin)
% LDL does LDL' Decomposition on a square matrix.
% IN:
%   A - matrix 
%   pivotstrategy - optional, no or block
% OUT:
%   L - lower triangle matrix
%   D - diagonal matrix
%   P - optional, permutation matrix 
% Usage:
%   No Pivoting: A = LDL', [L,D] = LDL(A)
%   Complete Diagonal Pivoting: P'AP = LDL', [L,D,P] = LDL(A,'block')
% Reference:
%   SIAM J. Numer. Anal., Vol. 8, 1971, 639-655.

% input arguments
if nargin == 1
    A = varargin{1};
    pivotstrategy = 'no';
elseif nargin == 2
    A = varargin{1};
    pivotstrategy = varargin{2};
    if ~strcmp(pivotstrategy,'block')
        error('Wrong Pivoting Strategy.\n');
    end
else
    error('Wrong Input Arguments.\n');
end

% A must be a symmetric matrix
[m,n] = size(A);
if m != n
    error('Input Matrix is not Square.\n');
end
for i=1:n
    for j=1:i-1
        if A(i,j) != A(j,i)
            error('Input Matrix is not Symmetric.\n');
            break;
        end
    end
end

% block criterion
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

% block pivoting
if strcmp(pivotstrategy,'block')
    if mu1 >= alpha*mu0
        % use 1x1 pivot
        P = eye(n); A_diag = diag(A);
        for k=1:n-1
            % find pivot index
            pivot = max(abs(A_diag(k:n)));
            for j=k:n
                if abs(A_diag(j)) == pivot
                    ind = j;
                end
            end

            % permutate A
            if j != k
                P(:,[k,ind]) = P(:,[ind,k]);
                A([k,ind],:) = A([ind,k],:);
                A(:,[k,ind]) = A(:,[ind,k]);
            end
        end
    else
        % use 2x2 pivot
        P = eye(n); P1 = eye(n); P2 = eye(n);
        for k=1:n-1
            % get part of A
            A_sub = A(k:n,k:n);

            % find pivot index
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

            % permutate A
            if (k!=row_ind) && (k+1!=col_ind)
                % interchange row/col k and row_ind
                P1(:,[k,row_ind]) = P1(:,[row_ind,k]);
                A([k,row_ind],:) = A([row_ind,k],:);
                A(:,[k,row_ind]) = A(:,[row_ind,k]);

                % interchange row/col k+1 and col_ind
                P2(:,[k+1,col_ind]) = P2(:,[col_ind,k+1]);
                A([k+1,col_ind],:) = A([col_ind,k+1],:);
                A(:,[k+1,col_ind]) = A(:,[col_ind,k+1]);

                % calculate transformation matrix P
                P = P*P1*P2;
            end
        end
    end
end

% LDL' decomposition
L = tril(A,-1)+eye(n); D = diag(A);
L(2:end,1) = L(2:end,1)/D(1);
for k=2:n
    v = L(k,1:k-1).*D(1:k-1)';
    D(k) = D(k) - L(k,1:k-1)*v';
    L(k+1:end,k) = (L(k+1:end,k)-L(k+1:end,1:k-1)*v')/D(k);
end
D = diag(D);

% Output Arguments.
varargout{1} = L;
varargout{2} = D;

if nargin == 2
    varargout{3} = P;
end

end
