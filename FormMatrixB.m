#!/usr/local/bin/octave
function [B] = FormMatrixB(A)
% FormMatrixB generates a 4x4 matrix from a 3x3 matrix A
% IN:
%   A: 3x3 matrix
% OUT:
%   B: 4x4 matrix
% Usage:
%   B = FormMatrix(A)

% get each element in A
a11 = A(1,1); a12 = A(1,2); a13 = A(1,3);
a21 = A(2,1); a22 = A(2,2); a23 = A(2,3);
a31 = A(3,1); a32 = A(3,2); a33 = A(3,3);

% direct but complex
%B = [a11+a22+a33, a23-a32, a31-a13, a12-a21; ...
%     a23-a32, a11-a22-a33, a12+a21, a13+a31; ...
%     a31-a13, a12+a21, a22-a11-a33, a23+a32; ...
%     a12-a21, a13+a31, a23+a32, a33-a11-a22];

% use the property of symmetric matrix
B = [a11+a22+a33, a23-a32,     a31-a13,     a12-a21; ...
     0.0,         a11-a22-a33, a12+a21,     a13+a31; ...
     0.0,         0.0,         a22-a11-a33, a23+a32; ...
     0.0,         0.0,         0.0,         a33-a11-a22];
B = triu(H,0) + triu(H,1)';
end
