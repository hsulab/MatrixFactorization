#!/usr/local/bin/octave
function [lambda1] = NewtonMethod(b,d)
% NewtonMethod calculates the dominant eigenvalue by Newton's Method.
% IN:
%   b - det(B)
%   d - det(A)
% OUT:
%   lambda1 - dominant eigenvalue
% Usage:
%   lambda1 = NewtonMethod(b,d)

% Coefficients of the charateristic polynomial formula.
coefs = [1,0,-2,-8*d,b];

% Coefficients of the derivative of charateristic polynomial formula.
coefs_d = [4,0,-4,-8*d];

% Set the initial guess x=sqrt(3).
x = sqrt(3); x_old = 3;

% Use Newton's Method to find lambda.
while x_old-x>10^-15
    % save old x
    x_old = x;

    % Horner's Method for calculating the value of polynomial formula.
    p = HornerPoly(coefs,x);
    p_d = HornerPoly(coefs_d,x);

    % update x
    x = x - p/p_d;
end

% get lambda1
lambda1 = x;

end
