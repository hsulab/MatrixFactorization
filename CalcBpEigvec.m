#!/usr/local/bin/octave
#####################################################
### Octave is a free software similiar to MATLAB. ###
#####################################################
function [w] = CalcBpEigvec(Bp)
    % Bp
    a = Bp(1,1); b = Bp(1,2);
    c = Bp(2,1); d = Bp(2,2);

    % lambda^2-(a+d)*lambda+ad-bc=0
    delta = (a+d)^2-4*(a*d-b*c);

    % Find w, eigenvector of smallest eigenvalue of Bp
    lambda = ((a+d)-sqrt(delta))/2;

    % w = [x1;x2] set x1=1
    x2 = (a+c-lambda)/(lambda-2*b);
    w = [1;x2];
    w = w/norm(w);

end

