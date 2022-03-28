clear; close all; clc;
syms D;

% Parámetros físicos.
p = 1.22;   v = 31.5;   Lc = 3;   h = 0.3;
M = 500;    g = 9.81;   cd = 2;
Rced = 248e6;   FS = 2;
% Parámetros del tubo.
L1 = 1.642;   L2 = 0.994;   u = 1/1.1;  %d/D.
c = D/2;   d = u*D;   I = (pi/64)*(D^4 - d^4);   
% Cálculo de áreas de ataque.
Af = (3*L1)*(4*L2);
Ac = (pi/4)*(D^2 - d^2);
% Esfuerzos cortante y flexionante.
Mf = (1/3)*(0.5*p*cd*(v^2)*Af)*Lc;
Sf = Mf*c/I;
Sc = -(M*g)/Ac;
% Cálculo del diámetro mínimo.
Snmax = Sf + Sc;
res = vpasolve(Snmax==(Rced/FS),D);
res = double(res);
sol = zeros(1,length(res));
for k=1:length(res)
    if ( isreal(res(k)) )
        sol(k) = res(k);
    end
end
% Se muestran los valores para el diámetro exterior e interior.
sol = sol(sol>0);
disp('Diametro mayor en [cm]: '); disp(sol'*100);
disp('Diametro menor en [cm]: '); disp(sol'*u*100);
Dfinal = sol*100/2.54;
aux = mod(Dfinal,0.25);
if( aux>0.125 )
    Dfinal = Dfinal + (0.25-aux);
else
    Dfinal = Dfinal - aux;
end
disp('Diametro mayor comercial en [cm]: '); disp(Dfinal*2.54);
disp('Diametro mayor comercial en [in]: '); disp(Dfinal);

% Se calcula el esfuerzo y factor de seguridad obtenidos con los
% diámetros calculados anteriormente.
clear D;
D = 8*2.54/100;
e = 6.35/1000;
d = D - 2*e;
c = D/2;   I = (pi/64)*(D^4 - d^4); 

Af = (3*L1)*(4*L2);
Ac = (pi/4)*(D^2 - d^2);

Mf = (1/3)*(0.5*p*cd*(v^2)*Af)*Lc;
Sf = Mf*c/I;
Sc = -(M*g)/Ac;

Snmax = Sf + Sc;
disp('Esfuerzo máximo en [MPa]: '); disp(Snmax);
disp('Factor de seguridad alcanzado: '); disp(Rced/Snmax);
