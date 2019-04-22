#!/usr/local/bin/octave
function [Q,H] = CompleteQH(A)
% CompleteQH compute the polar decomposition of A, which is A=QH.
% IN:
%   A - 3x3 matrix
% OUT:
%   Q - 3x3 orthogonal matrix
%   H - 3x3 symmetric positive semidefinite matrix
% Usage:
%   [Q,H] = CompleteQH(A)
% Reference:
%   Numer. Algor., 2016, 73, 349-369.

% normalize A
[A, A_fro] = ScaleA(A);

% tau2 Tolerance
tau2 = 1e-4;

% check sign of A
if det(A) >= eps
    eta = 1;
else
    eta = -1;
end
    
% form B from A
B = FormMatrixB(A);
    
% compute b=det(B) from an LU factorization with partial pivoting
[L,U,P] = LU(B,'partial');
b = det(L)*det(U)/det(P);
    
if b < 1-tau2
    % dominant eigenvalue of B is well separated
    %fprintf('Eigenvalues of B are well separated.\n');

    % compute d=detA using an LU factorization with partial pivoting
    [L,U,P] = LU(A,'partial');
    d = det(L)*det(U)/det(P);

    % change sign of B
    if d < 0
        B = -B; d = -d;
    end
    
    % estimate lambda1, a dominant eigenvalue of B
    lambda1 = EstDomiEigval(b,d);

    % get shifted matrix Bs
    Bs = lambda1*eye(length(B))-B;
        
    % LDL decomposition of Bs
    [L,D,P] = LDL(Bs,'block');
    
    % get approximate eigvec v
    v = P*(L^(-1))'*eye(4)(:,4);
    v = v/norm(v);
    
else
    % dominant eigenvalue of B is not well separated
    %fprintf('Eigenvalues of B are not well separated.\n');

    % compute d=detA using an LU factorization with partial pivoting
    [L,U,P] = LU(A,'complete');
    d = det(L)*det(U)/det(P);

    % change sign of B
    if d < 0
        B = -B;
    end
    
    % estimate lambda1
    lambda1 = EstDomiEigval(b,d);
    Bs = lambda1*eye(length(B)) - B;
    
    % Access to w by complete pivoting in P1AP2=LU, 
    % and evaluate w = -log10|u22|
    [L,U,P1,P2] = LU(A,'complete');
    u22 = U(2,2);
    
    % check iterative method
    if log10(abs(u22)) > -7.18
        % inverse iteration is fast
        %fprintf('Inverse Iteration.\n');
    
        % calculate number of iterations
        n_iterations = ceil(15/(16.86+2*log10(abs(u22))));

        % compute Bs = LDL' by block LDL' factorization
        [L,D,P] = LDL(Bs,'block');
    
        % initial guess for v
        v = P*(L^(-1))'*eye(4)(:,4);
        v = v/norm(v);
    
        % iterations
        for i=1:n_iterations
            % update v using one step of inverse iteration with LDL'
            v = (P^(-1))'*(L^(-1))*D*(L^(-1))'*(P^(-1))*v;
            v = v/norm(v);
        end
    
    else
        % subspace iteration is fast
        %fprintf('Inverse Subspace Iteration.\n');
    
        % compute Bs = LDL' by block LDL' factorization
        [L,D,P] = LDL(Bs,'block');
    
        % Initial guess for v.
        V = P*(L^(-1))'*eye(4)(:,3:4);

        % iteration
        for i=1:2
            % orthonormalize V via QR factorization
            [Q,R] = QR(V);
            V = Q(:,[1,2]);
    
            % update V using one step of inverse subspace iteration 
            % with LDL' factorization.
            V = (P^(-1))'*(L^(-1))*D*(L^(-1))'*(P^(-1))*V;
        end
    
        % Orthonormalize V via QR factorization.
        [Q,R] = QR(V);
        V = Q(:,[1,2]);
        Bp = V'*Bs*V;
    
        % Find w, eigenvector of smallest eigenvalue of Bp, by analytic formula.
        w = CalcBpEigvec(Bp);
        v = V*w;
    end
end

% form the matrix Q
Q = FormMatrixQ(v);

% check if sign of Q needs changing
Q = eta*Q;
    
% compute the upper triangle of H = Q'*A, 
% and set the lower triangle equal to the upper triangle
H = A_fro*Q'*A;
H = triu(H,0) + triu(H,1)';

end
