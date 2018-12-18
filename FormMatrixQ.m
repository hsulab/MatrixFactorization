#!/usr/local/bin/octave
#####################################################
### Octave is a free software similiar to MATLAB. ###
#####################################################
function [Q] = FormMatrixQ(v)
    % Form 3x3 Q from 4x1 v.

    % Get each element in v.
    v1 = v(1); v2 = v(2);
    v3 = v(3); v4 = v(4);

    % Form B.
    Q = [1-2*(v3^2+v4^2), 2*(v2*v3+v1*v4), 2*(v2*v4-v1*v3); ...
        2*(v2*v3-v1*v4), 1-2*(v2^2+v4^2), 2*(v3*v4+v1*v2); ...
        2*(v2*v4+v1*v3), 2*(v3*v4-v1*v2), 1-2*(v2^2+v3^2);];

end
