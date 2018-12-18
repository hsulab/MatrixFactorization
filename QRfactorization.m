#!/usr/local/bin/octave
#####################################################
### Octave is a free software similiar to MATLAB. ###
#####################################################
function [Q,R] = QRfactorization(A)
    % QR Factorization with Gram-Schmidt process.

    % Get th dimension of A.
    [m,n] = size(A);
    Q = zeros([m,m]);

    % Orthonormalize A to Generate Q.
    Q(:,1) = A(:,1)/norm(A(:,1));
    for j=2:n
        for k=1:j-1
            Q(:,j) = A(:,j);
            Q(:,j) = Q(:,j)-(A(:,j)'*Q(:,k))*Q(:,k);
        end
        Q(:,j) = Q(:,j)/norm(Q(:,j));
    end

    % Calculate R.
    R = zeros([m,n]);
    for j=1:n
        R(1:j,j) = Q(:,1:j)'*A(:,j);
    end

end
