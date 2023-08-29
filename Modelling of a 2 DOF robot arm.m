clear all;
clc;
close all;
syms t x y;
c = 8;
% Parameters-------------------------------
l = 1;
m = 1;

% State space ------------------------
A = [0 0 1 0;0 0 0 1;-0.4568 -0.6196 0 0;0.2485 -6.6174 0 0];
B = [0 0;0 0;0.7870 -0.0426;0.0426 0.1349];
C = [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
robotarm = ss(A,B,C,0);

%Jordan canonical form matrix-------------
canon(robotarm,"Modal");
canon_o = canon(robotarm,"Modal");

%Eigenvalues-----------------------------
eigA_o = eig(A);

%Controllability-------------------------
phi_c = ctrb(A,B);
rank(phi_c);

% Observability--------------------------
phi_o = obsv(A,C);
rank(phi_o);

% Realization----------------------------
canon_cont = canon(robotarm,"Companion");
transpose(canon_cont);

% Transfer function----------------------
tf(robotarm);
tf_o = tf(robotarm);

tf1 = tf_o(1,1);
tf2 = tf_o(1,2);
tf3 = tf_o(2,1);
tf4 = tf_o(2,2);
tf5 = tf_o(3,1);
tf6 = tf_o(3,2);
tf7 = tf_o(4,1);
tf8 = tf_o(4,2);

%Poles and zeros-------------------------
pole(tf_o);
tzero(tf(robotarm));

%Step response of the open loop system---
%figure(1)
step(robotarm)

%Step response of the closed loop system---
%Pole assignment
p1 = [-1,-2,-3,-4]; %near to zero
k1 = place(A,B,p1);
A1 = A-B*k1;
syscl = ss(A1,B,C,0);
%figure(2);
step(syscl);

p2 = 1000*[-1,-2,-3,-4]; %far from zero
k2 = place(A,B,p2);
A2 = A-B*k2;
syscl = ss(A2,B,C,0);
%figure(3);
step(syscl);

%Rootlocu----------------------------
%figure(4);
rlocus(tf1)
%figure(5);
rlocus(tf2)
%figure(6);
rlocus(tf3)
%figure(7);
rlocus(tf4)
%figure(8);
rlocus(tf5)
%figure(9);
rlocus(tf6)
%figure(10);
rlocus(tf7)
%figure(11);
rlocus(tf8)

%Observer designing
p3 = [-1, -2, -3, -4];
L = place(A',C',p3)';
A_obs = A-L*C;
B_obs = [B, L];
C_obs = [C;eye(4)];
sysObserver = ss(A_obs,B_obs,C_obs,0);

%Time response
%openloop
[y,t] = step(tf_o,8);
opts = timeoptions;
opts.Grid = 'on';
impulseplot(tf_o,opts);
%closedloop
[y,t] = step(syscl,8);
opts = timeoptions;
opts.Grid = 'on';
impulseplot(syscl,opts);
%observer
[y,t] = step(sysObserver,8);
opts = timeoptions;
opts.Grid = 'on';
impulseplot(sysObserver,opts);


%simulation of the robot
Q = [0; 0];
dQ = [0; 0];
T = [0; 0];
t = 0;
dt = 0.01;
while(t<10)

   c1 = cos(Q(1)); 
   s1 = sin(Q(1));
   c2 = cos(Q(2));
   s2 = sin(Q(2));

% Freefall
   T = [-0*dQ(1);-0*dQ(2)];
   ddQ = TwoLinkDynamics(Q, dQ, T);
   dQ = dQ + ddQ * dt;
   Q = Q + dQ*dt;
   t = t + dt;
% Plot the Robot
   x0 = 0;
   y0 = 0;
   x1 = cos(Q(1));
   y1 = sin(Q(1));
   x2 = x1 + cos(Q(1) + Q(2));
   y2 = y1 + sin(Q(1) + Q(2));
   clf;
   plot([-2,2],[-2,2],'x');
   hold on;
   plot([-2,2],[0,0], 'r');
   plot([0,0],[-2,2], 'r');
   plot([y0, y1, y2], [x0, x1, x2], 'b.-');
   pause(0.01);
end
function [ddQ] = TwoLinkDynamics( Q, dQ, T)
    q1 = Q(1);
    q2 = Q(2);
    dq1 = dQ(1);
    dq2 = dQ(2);
    c1 = cos(q1);
    c2 = cos(q2);
    s1 = sin(q1);
    s2 = sin(q2);
    s12 = sin(q1 + q2);
    c12 = cos(q1 + q2);
    M = [3+2*c2,1+c2;1+c2,1;]
    C = [2*s2*dq1*dq2 + s2*dq2*dq2; -s2*dq1*dq1];
    G = 10*[2*s1+s12;s12];
    ddQ = inv(M)*(T + C + G);
end
 

