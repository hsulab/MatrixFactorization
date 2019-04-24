#!/usr/local/bin/octave
% warning
warning('off')

% 3x3 Matrix A
A1 = [720,-650,710;396,-145,178;972,610,-529;];
A2 = [-25,300,300;70,-840,-840;-10,120,120;];

% generate A
y = 1;
A = 1/1275*(A1*y+A2);

% test
%A = [0.1,0.2,0.3;0.1,-0.1,0;0.3,0.2,0.1;];

%
s = 10
if s == 0
    tic
    for i=1:1
        [Q,H] = SVDQH(A);
    end
    toc

elseif s == 1
    tic
    for i=1:1
        [Q,H] = CompleteQH(A);
    end
    toc
end

cond(A)

% reference
[Q_ref,H_ref] = DirectQH(A);
S = svd(A)

% Higham
disp('Higham')
[Q,H] = CompleteQH(A);
fQ = norm(Q-Q_ref,'fro')/norm(Q_ref) % forward error in Q
fH = norm(H-H_ref,'fro')/norm(H_ref) % forward error in H
be = norm(A-Q*H,'fro')/norm(A,'fro') % back error 

% SVD
disp('SVD')
[Q,H] = SVDQH(A);
fQ = norm(Q-Q_ref,'fro')/norm(Q_ref) % forward error in Q
fH = norm(H-H_ref,'fro')/norm(H_ref) % forward error in H
be = norm(A-Q*H,'fro')/norm(A,'fro') % back error 
