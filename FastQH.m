#!/usr/local/bin/octave
#####################################################
### Octave is a free software similiar to MATLAB. ###
#####################################################
function [Q,H,A_fro] = FastQH(A)
    % This program compute the polar decomposition of A, which is A=QH.
    % This algorithem is modified from algorithm 3.5.
    
    % Initial Check.
    [m,n] = size(A);
    if (m!=n) || (m!=3) || (n!=3)
        error('Matrix A is not 3x3.\n');
    end

    % tau2 Tolerance
    tau2 = 1e-4;
    
    % Norm A.
    A_fro = norm(A,'fro');
    A = A/A_fro;
    
    % Form B from A
    B = FormMatrixB(A);
    
    % Compute b=det(B) from an LU factorization with partial pivoting
    [L,U,P] = LUpivot(B,'partial');
    b = det(L)*det(U)/det(P);
    
    if b < 1-tau2
        % Dominant eigenvalue of B is well separated.
        fprintf('Algo 3.2\n');

        % Compute d=detA using an LU factorization with partial pivoting.
        [L,U,P] = LUpivot(A,'partial');
        d = det(L)*det(U)/det(P);

        if d < 0
            B = -B;
            d = -d;
        end
    
        % Estimate lambda1, a dominant eigenvalue of B.
        lambda1 = EstimateLambda1(b,d);
        Bs = lambda1*eye(length(B))-B;
        
        [L,D,P] = LDLpivot(Bs,'block');
    
        v = (L^(-1))'*D(:,4);
        v = P*v/norm(v);
    
    else
        % Dominant eigenvalue of B is not well separated.
        fprintf('Algo 3.5\n');

        % Compute d=detA using an LU factorization with partial pivoting.
        [L,U,P] = LUpivot(A,'complete');
        d = det(L)*det(U)/det(P);

        if d < 0
            B = -B;
        end
    
        % Estimate lambda1 using algo33.
        lambda1 = EstimateLambda1(b,d);
        Bs = lambda1*eye(length(B)) - B;
    
        % Access to w by complete pivoting in P1AP2=LU,
        % and evaluate w = -log10|u22|
        [L,U,P1,P2] = LUpivot(A,'complete');
        u22 = U(2,2);
    
        if log10(abs(u22)) > -7.18
            fprintf('Inverse Iteration.\n');
    
            n_iterations = ceil(15/(16.86+2*log10(abs(u22))));

            % Compute Bs = LDL' by block LDL' factorization.
            [L,D] = LDLpivot(Bs);
    
            % Initial guess for v.
            v = (L^(-1))'*D(:,4);
            v = v/norm(v);
    
            for i=1:n_iterations
                % Update v using one step of inverse iteration with LDL'.
                v = (L^(-1))*D*(L^(-1))'*v;
                v = v/norm(v);
            end
    
        else
            fprintf('Inverse Subspace Iteration.\n');
    
            % Compute Bs = LDL' by block LDL' factorization.
            [L,D] = LDLpivot(Bs);
    
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
    
    end

    % Form the matrix Q.
    Q = FormMatrixQ(v);
    
    % Compute the upper triangle of H = Q'*A, 
    % and set the lower triangle equal to the upper triangle.
    H = Q'*A;
    H = triu(H,0) + triu(H,1)';

end
