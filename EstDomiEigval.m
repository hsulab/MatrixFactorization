#!/usr/local/bin/octave
function [lambda1] = EstDomiEigval(b,d)
% EstimateDominantEigenvalue computes an estimate of lambda1.
% IN:
%   b - detB
%   d - detA
% OUT:
%   lambda1 - dominant eigenvalue
% Usage:
%   lambda1 = EstimateDominantEigenvalue(b,d)

% Tolerance.
tau1 = 1e-4;

% choose method
if b + 1/3 > tau1
    % Analytic Charateristic Polynomial.
    lambda1 = PolyMethod(b,d);
else
    % Newton's Method.
    lambda1 = NewtonMethod(b,d);
end
end
