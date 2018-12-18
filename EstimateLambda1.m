#!/usr/local/bin/octave
#####################################################
### Octave is a free software similiar to MATLAB. ###
#####################################################
function [lambda1] = EstimateLambda1(b,d)
    % This algorithm computes an estimate of lambda1,
    % a dominant eogenvalue lambda1 of B given b=detB and d=detA.

    % Tolerance.
    tau1 = 1e-4;

    if b + 1/3 > tau1
        % Analytic Charateristic Polynomial.
        lambda1 = PolyMethod(b,d);
    else
        % Newton's Method.
        lambda1 = NewtonMethod(b,d);
    end

end
