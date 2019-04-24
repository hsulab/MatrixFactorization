#!/usr/local/bin/octave
function [Q,H] = SVDQH(A)
% SVDQH computes polar decomposition A with SVD.
% IN:
%   A - 3x3 matrix
% OUT:
%   Q - 3x3 orthogonal matrix
%   H - 3x3 symmetric positive semidefinite matrix
% Usage:
%   [Q,H] = SVDQH()

[U,S,V] = SVD(A);
Q = U*V';
H = V*S*V';

end
