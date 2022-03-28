clear; close all; clc;

materiales = ["Acero A36","Acero A529 Gr50","Acero A572 Gr50","Acero A992",...
    "Acero A500 GrB","Acero A500 GrC","Acero A501 GrA","Acero A53 GrB",...
    "Aluminio 6061 T4","Aluminio 6061 T6","Aluminio 6063 T4",...
    "Aluminio 6063 T5","Aluminio 6063 T6"]';
RcKSI = [36,50,50,50,42,46,36,35,18.855,39.162,11.893,20.306,30.46]';
RuKSI = [58,70,65,65,58,62,58,60,33.360,44.964,23.207,26.108,34.81]';
RcMPa = round(RcKSI*6.8944);
RuMPa = round(RuKSI*6.8944);
L = 2.1;       %[m].
F = 2000;    %[N].

% b = 101.60;  b = b/1000;
% h = 101.60;  h = h/1000;
% e = 2.39;   e = e/1000;
% Ac = b*h - (b-e)*(h-e);
% I = (1/12)*( b*h^3 - (b-2*e)*(h-2*e)^3 );
% c = h/2;

D = 88.90;   D = D/1000;
e = 6.35;     e = e/1000;
d = D - 2*e;
Ac = (pi/4)*(D^2 - d^2);
I = (pi/64)*( D^4 - d^4 );
c = D/2;

t = 18.7;
Fa = F*cosd(t);
Mf = (F*sind(t))*L;
Sa = Fa/Ac;
Sf = Mf*c/I;
SmaxMPa = ((Sa + Sf)/1000000)*ones(size(materiales));
FS = round(RcMPa./SmaxMPa,2);
disp(table(materiales,RuMPa,RcMPa,SmaxMPa,FS));




