% Cálculo de los parámetros del motor de CD  MTPM-P50-1L18.
clear; close all; clc;

%Información obtenida de la hoja de datos de AUTOMATIONDIRECT.
%Primer valor: SIN CARGA.
%Segundo valor: ESTIMADO.
%Tercer valor: MÁXIMA EFICIENCIA.
%Cuarto valor: MÁXIMA POTENCIA DE SALIDA.
%Quinto valor: MÁXIMO TORQUE.
%Sexto valor: ÚLTIMO VALOR.

R = 1.31;       %Resistencia de armadura [ohms].
Ja = 0.02365;   %Inercia del rotor [kg/m^2].
M = [0.077 2.115 2.092 3.684 3.684 3.684];      %Torque del motor [N*m].
Ia = [0.690 5.146 5.067 8.576 8.576 8.576];     %Corriente de armadura [A].
Va = [90.67 90.40 90.41 90.30 90.30 90.30];     %Voltaje de armadura [V].
w = (2*pi/60)*[1896 1693 1696 1551 1551 1551];  %Velocidad angular [rad/s].

%Cálculo de Ki: Tm = Ki*Ia -> Constante de torque [N*m/A].
Ki = mean(M(2:end)./Ia(2:end));

%Cálculo de Kb: Vb = Kb*w -> Constante de fuerza contraelectromotriz [V*s/rad].
%Vb: Fuerza contraelectromotriz
%Vb = Va - R*Ia (en estado estacionario).
Vb = Va(2:end) - R*Ia(2:end);
Kb = mean(Vb./w(2:end));

%Cálculo de Bm: Coeficiente de fricción viscosa [N*m*s].
%En estado estacionario sin carga: Bm*Wnom = Tm = Ki*Ia.
Bm = Ki*Ia(1)/w(1);

%Resultados.
fprintf("\nResistencia de armadura \nR = %g [ohms]\n",R);
fprintf("\nInercia del rotor \nJa = %g [kg/m^2]\n",Ja);
fprintf("\nConstante de torque \nKi = %g [V*s/rad]\n",Ki);
fprintf("\nConstante de fuerza electromotriz \nKb = %g [N*m/A]\n",Kb);
fprintf("\nCoeficiente de fricción viscosa \nBm = %e [N*m*s]\n",Bm);