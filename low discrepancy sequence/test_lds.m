clear all,close all, clc;


p = 2;     % base
Np = 100;  % Nb projections
Ns = 1;    % Partition of the circle into Ns segments
Nz = 50;   % Add a projection at a reference angle each Nz projections

Seq = Low_discrepency_squence(Np,Ns,p,Nz);


theta = 0:pi/1000:2*pi;
xcircle = cos(theta);
ycircle = sin(theta);


x = cos(Seq*pi/180);
y = sin(Seq*pi/180);

h = plot(xcircle, ycircle);
hold on;
h0 = plot(x, y, 'ro');
axis off;
axis square;