%Calcula errores máximos de los controladores PD, PID y GPI.
clear; close all; clc;
% Ejecuta los modelos de simulink para obtener los datos.
sim('robotSeguidor_PD.slx');
sim('robotSeguidor_PID.slx');
sim('robotSeguidor_GPI.slx');
% Conversiones de unidades.
eq1_PD = eq1_PD*180/pi;
eq2_PD = eq2_PD*180/pi;
eq1_PID = eq1_PID*180/pi;
eq2_PID = eq2_PID*180/pi;
eq1_GPI = eq1_GPI*180/pi;
eq2_GPI = eq2_GPI*180/pi;
w_PD = max(w_PD*4800*60/(2*pi));
w_PID = max(w_PID*4800*60/(2*pi));
w_GPI = max(w_GPI*4800*60/(2*pi));
V_PD = max(V_PD);
V_PID = max(V_PID);
V_GPI = max(V_GPI);
Tm_PD = max(Tm_PD);
Tm_PID = max(Tm_PID);
Tm_GPI = max(Tm_GPI);

% Errores finales de q1.
eq1_PD = eq1_PD(end);
eq1_PID = eq1_PID(end);
eq1_GPI = eq1_GPI(end);
% Errores finales de q2.
eq2_PD = eq2_PD(end);
eq2_PID = eq2_PID(end);
eq2_GPI = eq2_GPI(end);
% Voltaje de armadura para q1.
Vq1_PD = V_PD(end,1);
Vq1_PID = V_PID(end,1);
Vq1_GPI = V_GPI(end,1);
% Voltaje de armadura para q2.
Vq2_PD = V_PD(end,2);
Vq2_PID = V_PID(end,2);
Vq2_GPI = V_GPI(end,2);
% Torque demandado para el motor de q1.
Tmq1_PD = Tm_PD(end,1);
Tmq1_PID = Tm_PID(end,1);
Tmq1_GPI = Tm_GPI(end,1);
% Torque demandado para el motor de q2.
Tmq2_PD = Tm_PD(end,2);
Tmq2_PID = Tm_PID(end,2);
Tmq2_GPI = Tm_GPI(end,2);
% Velocidad angular de q1.
wq1_PD = w_PD(end,1);
wq1_PID = w_PID(end,1);
wq1_GPI = w_GPI(end,1);
% Velocidad angular de q2.
wq2_PD = w_PD(end,2);
wq2_PID = w_PID(end,2);
wq2_GPI = w_GPI(end,2);

parametro = ["errorFinal_q1 [°]","errorFinal_q2 [°]","w_q1 [rpm]",...
    "w_q2 [rpm]","V_q1 [V]","V_q2 [V]","Tm_q1 [Nm]","Tm_q2 [Nm]"]';
PD = [eq1_PD,eq2_PD,wq1_PD,wq2_PD,Vq1_PD,Vq2_PD,Tmq1_PD,Tmq2_PD]';
PID = [eq1_PID,eq2_PID,wq1_PID,wq2_PID,Vq1_PID,Vq2_PID,Tmq1_PID,Tmq2_PID]';
GPI = [eq1_GPI,eq2_GPI,wq1_GPI,wq2_GPI,Vq1_GPI,Vq2_GPI,Tmq1_GPI,Tmq2_GPI]';
disp(table(parametro,PD,PID,GPI));
