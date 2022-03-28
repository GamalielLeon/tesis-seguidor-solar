clear; close all; clc;
% Adquisición de información generada por el ASS GRENA 5 (2012).
M = xlsread('SV.xlsx'); % Lectura de archivo excel.
N = M(1);
M = M(2:end,:);
h = M(:,1) - 1.6;
q2 = M(:,2);
q1 = M(:,3);
delta = M(:,4);
q1(q1<=0) = 360 + q1(q1<=0);
h = h(q2>=0);           % Hora del día.
delta = delta(q2>=0);   % Declinación solar.
q1 = -q1(q2>=0);        % Ángulo azimutal.
q2 = q2(q2>=0);         % Ángulo de elevación.

% Gráficas de las coordenadas (x,y,z) obtenidas con los datos anteriores.
figure(1);
% Gráfica de las coordenadas espaciales del seguidor solar.
subplot(2,2,1); 
L2 = 0.3; L1 = 0.382; L0 = 3 - L1;
x = L2*cosd(q1).*cosd(q2);
y = L2*sind(q1).*cosd(q2);
z = L0 + L1 + L2*sind(q2);
plot3(x,y,z,'.','MarkerSize',12); grid on;
title('Desplazamiento espacial teórica del seguidor');
xlabel('X'); ylabel('Y'); zlabel('Z');
% Gráfica de las coordenadas espaciales del seguidor normalizadas (unitarias).
subplot(2,2,2);
L2 = 0.3; L1 = 0; L0 = 0; 
x2 = L2*cosd(q1).*cosd(q2);
y2 = L2*sind(q1).*cosd(q2);
z2 = L0 + L1 + L2*sind(q2);
norma = sqrt(x2.^2 + y2.^2 + z2.^2);
plot3(x2./norma,y2./norma,z2./norma,'.r','MarkerSize',12);
title('Coordenadas del seguidor normalizadas'); grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
% Comparación entre las curvas de las coordenadas normalizadas del seguidor
% y las ecuaciones vectoriales del movimiento solar (vector unitario).
subplot(2,2,3);
lambda = 0.3405426746;      % Latitud.
w = (h - 12)*15*pi/180;
e = 23.44*pi/180;           % Declinación natural.
sN = -sin(lambda)*cos(w).*cos(delta) + cos(lambda)*sin(delta);
sE = -sin(w).*cos(delta);
sZ = cos(lambda)*cos(w).*cos(delta) + +sin(lambda)*sin(delta);
xx2 = -double(sN);
yy2 = double(sE);
zz2 = double(sZ);
plot3(xx2,yy2,zz2,'k','linewidth',1.5); grid on; hold on;
plot3(x2./norma,y2./norma,z2./norma,'.r','MarkerSize',12);
title('Coordenadas normalizadas (seguidor vs ecuaciones vectoriales)');
xlabel('X'); ylabel('Y'); zlabel('Z'); hold off;
% Comparación entre las curvas de las coordenadas reales del seguidor
% y las ecuaciones vectoriales del movimiento solar (vector igual al del seguidor).
subplot(2,2,4);
L2 = 0.3; L1 = 0.382; L0 = 3 - L1;
xx = -L2*double(sN);
yy = L2*double(sE);
zz = L2*double(sZ) + L0 + L1;
plot3(xx,yy,zz,'k','linewidth',1.5); grid on; hold on;
plot3(x,y,z,'.','MarkerSize',12);
title('Coordenadas espaciales (seguidor vs ecuaciones vectoriales)');
xlabel('X'); ylabel('Y'); zlabel('Z'); hold off;

%% Configuracion del robot
%theta d a alpha sigma offset
L2 = 0.3; L1 = 0.382; L0 = 3 - L1;
L(1)=Link([0 L0+L1 0 pi/2 0 0]);   %parametros de la primer junta y eslabón.
L(2)=Link([0 0 L2 0 0 0]);         %parametros de la segunda junta y eslabón.
qli=[-pi pi ; 0 pi/2];             %definición de los limites de las juntas 1 y 2.
qz=[0 0];                          %matriz con vectores articulares renglon.
antropom=SerialLink(L,'name','antrop','qlim',qli); %construccion del objeto polar, 
antropom.plotopt = {'workspace', [-3 4 -3 4 0 5]}; %opcion del metodo plot,
%antropom.teach(qz)
figure(3);
%plot3(xx(1),yy(1),zz(1),'*k','linewidth',1); hold on;
plot3(xx,yy,zz,'r','linewidth',3);
antropom.plot([q1 q2]*pi/180,'fps',10,'trail',{'.','MarkerSize', 4});
% % %T=antropom.fkine(qfinal)
