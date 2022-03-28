function ddq = fcn(u,q,dq,d)
%Variables de junta y sus derivadas.
q1 = q(1); q2 = q(2);
dq1 = dq(1); dq2 = dq(2);

%Parámetros físicos eslabon 1.
m1 = 152;
I1x = 276.66;
I1y = 147;
I1z = 151.1;
L1 = 0.382;
Lcm1 = L1/2;
%Parámetros físicos eslabon 2.
m2 = 382;
I2x = 463.64;
I2y = 1352.36;
I2z = 902.58;
L2 = 0.3;
Lcm2 = 0.297;
%Parámetros generales.
g = 9.81;
L0 = 3 - L1;
%Parámetros del motor.
R = 1.31;
Bm = 1.4683e-3;
Ki = 0.4225;
Kb = 0.4808;
Jm1 = 0.05305;
Jm2 = 0.07135;
r = 100*48;

%Variables auxiliares para simplificar.
Ix_y = I2x - I2y;
Ixy = I2x + I2y;
m2Lcm2 = m2*Lcm2^2;

%Matriz de inercia.
d11 = 0.5*(  m2Lcm2 + 2*I1y + Ixy + (m2Lcm2 - Ix_y)*cos(2*q2)  );
d12 = Ix_y*cos(q2)*sin(q1)*sin(q2);
d21 = d12;
d22 = I2z*cos(q1)^2 + 0.5*( 2*m2Lcm2 + I2y - I2y*cos(2*q1) )*sin(q2)^2 +...
    ( m2Lcm2 + I2x*sin(q1)^2 )*cos(q2)^2;
D = [ d11 d12 ; d21 d22 ];

%Matriz de coriolis.
c11 = ( -m2Lcm2 + Ix_y )*cos(q2)*sin(q2)*dq2;
c12 = ( -m2Lcm2 + Ix_y )*cos(q2)*sin(q2)*dq1 - 0.25*(  (Ixy - 2*I2z)*sin(2*q1)...
    + Ix_y*cos(2*q2)*( sin(2*q1) - 4*sin(q1) )  )*dq2;
c21 = 0.5*( m2Lcm2 - Ix_y + Ix_y*cos(q1) )*sin(2*q2)*dq1 + 0.25*( Ixy - 2*I2z...
    + Ix_y*cos(2*q2) )*sin(2*q1)*dq2;
c22 = -Ix_y*cos(q2)*sin(q2)*dq2*sin(q1)^2 + 0.25*( Ixy - 2*I2z +...
    Ix_y*cos(2*q2) )*sin(2*q1)*dq1;
C = [ c11 c12 ; c21 c22 ];

%Matriz de gravedad.
g1 = 0;
g2 = g*Lcm2*m2*cos(q2);
G = [ g1 ; g2 ];

%Salida del vector de aceleración.
Jeff = [Jm1 0 ; 0 Jm2]*r + (1/r)*D;
Beff = (Bm + Kb*Ki/R)*r*eye(2) + (1/r)*C;
% Jeff*ddq + Beff*dq + (1/r)*G = u;
ddq = Jeff\( u - Beff*dq - (1/r)*G - (d/r)*[1 1]' );
