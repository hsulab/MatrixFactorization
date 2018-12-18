#!/usr/local/bin/octave
#####################################################
### Octave is a free software similiar to MATLAB. ###
#####################################################
function [lambda1] = NewtonMethod(b,d)
    % This program calculates the dominant eigenvalue by Newton's Method.
    % b=det(B), d=det(A)

    % coefs :=
    % Coefficients of the charateristic polynomial formula.
    coefs = [1,0,-2,-8*d,b]

    % coefs_d :=
    % Coefficients of the derivative of charateristic polynomial formula.
    coefs_d = [4,0,-4,-8*d]

    % Set the initial guess x=sqrt(3).
    x = sqrt(3); x_old = 3;

    % Use Newton's Method to find lambda.
    while x_old-x>10^-15
        x_old = x;

        % Horner's Method for Calculating the value of polynomial formula.
        p = HornerPoly(coefs,x);
        p_d = HornerPoly(coefs_d,x);

        % Update x.
        x = x - p/p_d;
    end

    lambda1 = x;

end
