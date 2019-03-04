function [g] = bungeRotationSample2Crystal(angles)
% rotation from sample to crystal axes

phi1 = angles(1);
Phi = angles(2);
phi2 = angles(3);


Z1 = [cos(phi1) sin(phi1) 0;
     -sin(phi1) cos(phi1) 0;
         0         0      1];
     
X  = [   1         0      0
         0   cos(Phi) sin(Phi);
         0  -sin(Phi) cos(Phi);];
     
Z2 = [cos(phi2) sin(phi2) 0;
     -sin(phi2) cos(phi2) 0;
         0         0      1];

g = Z2*X*Z1;