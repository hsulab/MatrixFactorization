#!/usr/local/bin/octave
function [a] = DET(A)
% DET computes the det of a triangle square matrix.
% IN:
%   A - triangle square matrix
% OUT:
%   a - determinant value
% Usage:
%   a = DET(A);
% Notes:
%   The det of a triangle matrix equals its diagonal entries' multiplication.

% check square
[m,n] = size(A);
if m ~= n
    error('Invalid Input Matrix: Not Square.')
end

% det
a = 1;
for i=1:n
    a = a*A(i,i);
end

end
