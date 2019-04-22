#!/usr/local/bin/octave
function [A_scaled,A_fro] = ScaleA(A)
% ScaleA normalize A with F-norm.
% IN:
%   A - input 3x3 matrix
% OUT:
%   A_scaled - 3x3 matrix normalized with F-norm
%   A_fro - F-norm of A
% Usage:
%   [A_scaled,A_fro] = ScaleA(A)

% Initial Check.
[m,n] = size(A);
if (m!=n) || (m!=3) || (n!=3)
    error('Matrix A is not 3x3.\n');
end

% normalize input matrix A with built-in function
A_fro = norm(A,'fro');
A_scaled = A/A_fro;

end
