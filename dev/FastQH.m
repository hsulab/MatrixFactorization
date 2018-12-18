#!/usr/local/bin/octave
#####################################################
### Octave is a free software similiar to MATLAB. ###
#####################################################
% This program compute the polar decomposition of A, which is A=QH.
% This algorithem is modified from algorithm 3.5.

% tau2 Tolerance
tau2 = 1e-4;

% Norm A.
%A=[4,-1,1;-1,4.25,2.75;1,2.75,3.5];
A = [0.1 0.2 0.3;0.1,-0.1,0;0.3,0.2,0.1];
A_fro = norm(A,'fro');
A = A/A_fro;

% Form B from A
B = FormMatrixB(A);

% Compute b=detB from an LU factorization with partial pivoting
b = LuPivotDet(B);

if b < 1-tau2
    disp('b is smaller than 0. -> Dominant eigenvalue of B is well separated.')
    % Dominant eigenvalue of B is well separated.
    % Compute d=detA using an LU factorization with partial pivoting.
    d = LuPivotDet(A);

    if d < 0
        B = -B;
        d = -d;
    end

    % Estimate lambda1, a dominant eigenvalue of B.
    lambda1 = EstimateLambda1(b,d);

    Bs = lambda1*eye(length(B))-B;
    
    [L,D,P] = LDLdiagonalPivot(Bs);

    v = (L^(-1))'*D(:,4);
    v = P*v/norm(v);

    % Form the matrix Q.
    Q = FormMatrixQ(v);

    % Compute the upper triangle of H = Q', 
    % and set the lower triangle equal to the upper triangle.
    H = Q'*A;
    H = tril(H)+tril(H,-1)';

else
    % Compute d=detA using an LU factorization with complete pivoting.
    d = LuPivotDet(A);

    if d < 0
        B = -B;
    end

    % Estimate lambda1 using algo33
    lambda1 = EstimateLambda1(b,d);
    Bs = lambda1*eye(length(B)) - B;

    % Access to w by complete pivoting in P1AP2=LU,
    % and evaluate w = -log10|u22|
    [L,U,P1,P2] = LUcompletePivot(A);
    u22 = U(2,2);

    if log10(abs(u22)) > -7.18
        disp('Inverse Iteration.')

        n_iterations = ceil(15/(16.86+2*log10(abs(u22))));
        % Compute Bs = LDL' by block LDL' factorization with Bunch-Parlett.
        [L,D] = LDLnoPivot(Bs);

        % Initial guess for v.
        v = (L^(-1))'*D(:,4);
        v = v/norm(v);

        for i=1:n_iterations
            % Update v using one step of inverse iteration with LDL'.
            v = (L^(-1))*D*(L^(-1))'*v;
            v = v/norm(v);
        end

    else
        disp('Inverse Subspace Iteration.')

        % Compute Bs = LDL' by block LDL' factorization with Bunch-Parlett.
        [L,D] = LDLnoPivot(Bs);

        % Initial guess for v.
        V = (L^(-1))'*[D(:,3) D(:,4)];
        for i=1:2
            % Orthonormalize V via QR factorization
            [Q,R] = QRfactorization(V);
            V = Q(:,[1,2]);

            % Update V using one step of inverse subspace iteration 
            % with LDL' Factorization.
            V = (L^(-1))*D*(L^(-1))'*V;
        end

        % Orthonormalize V via QR factorization.
        [Q,R] = QRfactorization(V);
        V = Q(:,[1,2]);
        Bp = V'*Bs*V;

        % Find w, eigenvector of smallest eigenvalue of Bp, by analytic formula.
        w = CalcBpEigvec(Bp);
        v = V*w;
    end

    % Form the matrix Q.
    Q = FormMatrixQ(v);

    % Compute the upper triangle of H = Q'*A, 
    % and set the lower triangle equal to the upper triangle.
    H = Q'*A;
    H = tril(H)+tril(H,-1)';

end
