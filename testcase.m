#!/usr/local/bin/octave
#####################################################
### Octave is a free software similiar to MATLAB. ###
#####################################################

% 3x3 Matrix A
A = [4,-1,1;-1,4.25,2.75;1,2.75,3.5];

% Basic QH Algo
[Q,H,f] = BasicQH(A);
A/f
Q*H

% Fast QH Algo
[Q,H,f] = FastQH(A);
A/f
Q*H
