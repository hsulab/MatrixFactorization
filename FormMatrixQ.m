#!/usr/local/bin/octave
function [Q] = FormMatrixQ(v)
% FormMatrixQ generates 3x3 Q from 4x1 v.
% IN:
%   v - 4x1 vector
% OUT:
%   Q - 3x3 matrix
% Usage:
%   Q = FormMatrixQ(v)

% Get each element in v.
v1 = v(1); v2 = v(2);
v3 = v(3); v4 = v(4);

% Form B.
Q = [1-2*(v3^2+v4^2), 2*(v2*v3+v1*v4), 2*(v2*v4-v1*v3); ...
    2*(v2*v3-v1*v4), 1-2*(v2^2+v4^2), 2*(v3*v4+v1*v2); ...
    2*(v2*v4+v1*v3), 2*(v3*v4-v1*v2), 1-2*(v2^2+v3^2);];

end
