#!/usr/local/bin/octave
function [Q,H] = DirectQH(A)
% DirectQH computes polar decomposition of A with definition.
% IN:
%   A - 3x3 matrix
% OUT:
%   Q - 3x3 orthogonal matrix
%   H - 3x3 symmetric positive semidefinite matrix
% Usage:
%   [Q,H] = DirectQH(A)

H = sqrtm(A'*A);
Q = A*inv(H);

end
