#!/usr/local/bin/octave
function [Q,R] = QR(A, strategy)
% QR computes A=QR with Gram-Schmidt Orthogonalization or Householder transformation.
% IN:
%   A - square matrix
% OUT:
%   Q - orthogonal matrix
%   R - upper triangle matrix
% Usage:
%   [Q,R] = QR(A)

% get strategy
if ~exist('strategy', 'var')
    strategy = 'gs'; % default: Gram-Schmidt
end

% get dimension of A
[m,n] = size(A);

% get the dimension of A
if strcmp(strategy,'gs')
    % init
    Q = zeros(m,n);
    R = zeros(n);
    
    % Gram-Schmidt
    for k=1:n
        % check full rank
        R(k,k) = norm(A(:,k)); % flops: nM+(n-1)A+1sqrt
        if R(k,k) == 0
            error('Input matrix must be full rank.')
        end
    
        % check 
        Q(:,k) = A(:,k)/R(k,k); % flops: nM
        for j=k+1:n
            R(k,j) = Q(:,k)'*A(:,j); % flops: nM+(n-1)A
            A(:,j) = A(:,j) - R(k,j)*Q(:,k); % flops: nM + nA
        end
    end
elseif strcmp(strategy,'ht')
    % init 
    E = eye(n);
    X = zeros(n,1);
    R = zeros(n);
    P1 = E;
    for k=1:n-1
        % get w, Pk=1-2ww'
        s = -sign(A(k,k))*norm(A(k:n,k));
        R(k,k) = -s;
        if k == 1
            w = [A(1,1)+s,A(2:n,k)']';
        else 
            w = [zeros(1,k-1),A(k,k)+s,A(k+1:n,k)']';
            R(1:k-1,k) = A(1:k-1,k);
        end
        if norm(w) ~= 0
            w = w/norm(w);
        end 
        P = E - 2*w*w';
        A = P*A;
        P1 = P*P1;
        R(1:n,n) = A(1:n,n);
    end
    Q = P1';
else 
    error('Wrong strategy')

end
