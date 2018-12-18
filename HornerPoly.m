#!/usr/local/bin/octave
#####################################################
### Octave is a free software similiar to MATLAB. ###
#####################################################
function polyvalue = HornerPoly(A,x)
    % Horner's Method is fast to calculate polynomial value.

    % Set the initial value as 0.
    polyvalue = 0;

    % Determine the degree of coefficient matrix A.
    % a_n*x^n + a_(n-1)*x^(n-1) + ... + a0.
    % A = [a_n, a_(n-1), ..., a0].
    d = length(A);

    % Calculate the value of the polynomial.
    for i=1:d
        polyvalue = A(i)+polyvalue*x;
    end

end

