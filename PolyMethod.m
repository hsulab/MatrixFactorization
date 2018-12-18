#!/usr/local/bin/octave
########################################################################
######               The code is the same for MATLAB.             ######
########################################################################
function [lambda1] = PolyMethod(b,d)
    %%% Documentation for Algorithm 3.3 %%%
    % This algorithm computes an estimate of lambda1.
    % lambda1 is the dominant eigenvalue of B.

    % The characteristic polynomial is 'x^4-2x^2-8dx+b=0'.
    % LU partial pivoting is used to calculate b = detB and d = detA.

    %%% Ferrari-Lagrange Method %%%
    % To solve the quartic equation, Ferrari-Lagrange Method is used.

    % The standard quartic equation is 
    % 'x^4+ax^3+bx^2+cx+d=(x^2+g1x+h1x)(x^2+g2x+h2)=0' where 
    % g1+g2=a, g1g2+h1+h2=b, g1h2+g2h1=c, h1h2=d.

    % The standard transformed cubic equation is 
    % 'y^3-by^2+(ac-4d)y-a^2d-c^2+4bd=0' where
    % g^2-ag+b-y=0, h^2-yh+d=0.

    % The transformed cubic equation is 'y^3+2y^2-4by-64d^2-8b=0'.

    %%% Trigonometric Function %%%
    % To solve the cubic equation, trigonometric function is used.

    % For standard fomula ax^3+bx^2+cx+d=0 (a!=0),
    % set delta=(-b^3/27a^3+-d/2a+bc/6a^2)^2+(c/3a-b^2/9a^2)^3
    %          =alpha^2+beta^3<0

    % root x=-b/3a+2sqrt(-beta)cos(arccos(alpha/(-beta)^3/2+k)/3)
    %      k=0,2pi,-2pi

    %%% Compute the lambda1 %%%
    % for simple
    c = 8*d;

    % beta=-4/9*(1+3b)
    delta0 = 1+3*b;

    % alpha=8/27*(-1+108d^2+9b)
    delta1 = -1+(27/16)*c^2+9*b;

    % (8/27)/(--4/9)^(3/2)=1
    alpha = (delta1)/(delta0^(3/2));

    % y=z-2
    z = (4/3)*(1+delta0^(1/2)*cos(acos(alpha)/3));

    % g1=sqrt(z), g2=-sqrt(z)
    s = z^0.5/2;

    % h1=(y+sqrt(y^2-4b))/2, h2=(y-sqrt(y^2-4b))/2.
    % The maximum lambda is (-g+sqrt(g^2-4h))/2.
    % 4-z+c/s=sqrt(g^2-4h).
    lambda1 = s+(max(0,4-z+c/s))^(1/2)/2;
end
