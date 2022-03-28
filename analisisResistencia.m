function varargout = analisisResistencia(varargin)
% ANALISISRESISTENCIA MATLAB code for analisisResistencia.fig
%      ANALISISRESISTENCIA, by itself, creates a new ANALISISRESISTENCIA or raises the existing
%      singleton*.
%
%      H = ANALISISRESISTENCIA returns the handle to a new ANALISISRESISTENCIA or the handle to
%      the existing singleton*.
%
%      ANALISISRESISTENCIA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANALISISRESISTENCIA.M with the given input arguments.
%
%      ANALISISRESISTENCIA('Property','Value',...) creates a new ANALISISRESISTENCIA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before analisisResistencia_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to analisisResistencia_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help analisisResistencia

% Last Modified by GUIDE v2.5 20-Nov-2019 13:33:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @analisisResistencia_OpeningFcn, ...
                   'gui_OutputFcn',  @analisisResistencia_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before analisisResistencia is made visible.
function analisisResistencia_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to analisisResistencia (see VARARGIN)

% Choose default command line output for analisisResistencia
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes analisisResistencia wait for user response (see UIRESUME)
% uiwait(handles.figure1);
set(handles.panelA1,'visible','off');
set(handles.panelA2,'visible','off');
set(handles.panelAnalisis,'visible','on');
set(handles.salir,'visible','on');
set(handles.tabla,'visible','off');
handles.analisis1.Value = 0;
handles.analisis2.Value = 0;

set(handles.textFzas,'visible','off');
set(handles.textDist,'visible','off');
set(handles.distancias,'visible','off');
set(handles.fuerzas,'visible','off');

texto = ["Acero A36","Acero A529 Gr50","Acero A572 Gr50","Acero A992",...
"Acero A500 GrB","Acero A500 GrC","Acero A501 GrA","Acero A53 GrB",...
"Aluminio 6061 T4","Aluminio 6061 T6","Aluminio 6063 T4",...
"Aluminio 6063 T5","Aluminio 6063 T6"];
rc = [36,50,50,50,42,46,36,35,18.855,39.162,11.893,20.306,30.46];
ru = [58,70,65,65,58,62,58,60,33.360,44.964,23.207,26.108,34.81];
handles.text = texto;
handles.rcksi = rc;
handles.ruksi = ru;
set(handles.listMat,'String',texto);
set(handles.tablaA2,'visible','off');
guidata(hObject,handles);


% --- Outputs from this function are returned to the command line.
function varargout = analisisResistencia_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% --- Executes on button press in salir.
function salir_Callback(hObject, eventdata, handles)
delete(handles.output);

function dimPanel_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function dimPanel_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FS_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function FS_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Ce_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function Ce_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mp_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function mp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function vmax_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function vmax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in calcA1.
function calcA1_Callback(hObject, eventdata, handles)
materiales = handles.text';
%Resistencia de cedencia [ksi] y [MPa].
RcKSI = handles.rcksi';  RcMPa = RcKSI*6.8944;
%Resistencia de cedencia [ksi] y [MPa].
RuKSI = handles.ruksi';  RuMPa = RuKSI*6.8944;

obj = allchild(handles.panelMat);
aux = flipud(cell2mat(get(obj,'Value')));
materiales = materiales(aux==1);
RcMPa = RcMPa(aux==1);
RuMPa = RuMPa(aux==1);

g = 9.81;           %constante de aceleraci?n de la gravedad [m/s^2].
L1p = 1.642;        %Largo del panel (dimension mas larga) [m].
L2p = 0.994;        %Ancho del panel (dimension mas corta) [m].
ep = 0.04;          %Grosor o espesor del panel [m].
bp = 0.03;          %Ancho de la superficie de contacto del panel.
mp = str2double(get(handles.mp,'String'));      %Masa de cada panel [Kg].
FS = str2double(get(handles.FS,'String'));      %Factor de seguridad.
vmax = str2double(get(handles.vmax,'String'));  %Velocidad maxima del viento [m/s].
Ce = str2double(get(handles.Ce,'String'));      %Coeficiente eolico.
Fp = mp*g;                                      %Peso de cada panel [N].

%pAire = P/(RT)     %P:presi?n, R:constanteAireSeco=287.05, T:temperatura.
pAire = 1.22;       %Densidad te?rica del aire.
Aefp = L1p*L2p;         %?rea efectiva de contacto del panel.
Aefv = (3*L1p)*(4*L2p); %?rea efectiva de contacto del viento.
Fv = Ce*(0.5*pAire*(vmax^2)*Aefv);  %Fviento = Ce[(1/2)*p*v^2]*A.
bandFp = handles.cargaPanel.Value;
bandFv = handles.cargaViento.Value;
band2 = handles.bandFuerzas.Value;
dist = []; ff = [];
syms x;
La = str2double(handles.distLa.String); %Distancia desde el borde izquierdo de la viga al primer apoyo.
Lb = str2double(handles.distLb.String);	%Distancia desde el borde derecho de la viga al segundo apoyo.
L1 = str2double(handles.Lviga.String);	%Longitud de la viga.
L2 = 4*L2p;

if( sum(isnan([FS Ce mp vmax])) > 0 )
    msgbox('Un campo de datos está vacío o posee datos invalidos','Advertencia','warn');
    return;
end
if( sum(isnan([La Lb]) > 0) )
    msgbox('El campo de La, Lb o L está vacío o posee datos invalidos','Advertencia','warn');
    return;
end
if( (La<0)||(Lb>L1) )
    msgbox('El valor de La o de Lb está fuera de la longitud de la viga','Advertencia','warn');
    return;
end
if(band2 == 1)
    dist = str2num(handles.distancias.String);
    ff = str2num(handles.fuerzas.String);
end
if(max(dist) > L1)
    msgbox('Una distancia excede la longitud de la viga','Advertencia','warn');
    return;
end
if(length(dist) ~= length(ff))
    msgbox('El número de fuerzas y distancias debe ser igual','Advertencia','warn');
    return;
end
if( (band2==1)&&(isempty(dist)||isempty(ff)) )
    msgbox('El campo de distancias o fuerzas está vacío o posee datos invalidos','Advertencia','warn');
    return;
end

% Obtención de las reacciones de los dos apoyos de la viga, considerando
% todas las fuerzas presentes en la viga (cargas puntuales y lineales).
resa = abs(dot(dist-La,ff));
resb = abs(dot(dist-(L1-Lb),ff));
wv = bandFv*(0.2*Fv/Aefv)*L2;
wp = bandFp*(Fp/Aefp)*L2p;
Ra1 = ((wv*L1)*(0.5*L1-Lb))/(L1-La-Lb);
Rb1 = ((wv*L1)*(0.5*L1-La))/(L1-La-Lb);
Ra2 = ((wp*L1)*(0.5*L1-Lb) + resb)/(L1-La-Lb);
Rb2 = ((wp*L1)*(0.5*L1-La) + resa)/(L1-La-Lb);
ResF = @(x) sum(double((ff.*(dist - x)).*(x>=dist)));
M1 = @(x) double( Ra1*(x - La).*(x>=La) + Rb1*(x - L1 + Lb).*(x>=(L1-Lb)) - wv*(0.5*x.^2) );
M2 = @(x) double( Ra2*(x - La).*(x>=La) + Rb2*(x - L1 + Lb).*(x>=(L1-Lb)) - wp*(0.5*x.^2) + ResF(x) );
xdat = 0:0.001:L1;
datM2 = zeros(1,length(xdat));
for i=1:length(xdat)
    datM2(i) = M2(xdat(i));
end
% Obtención del momento flexionante máximo y calculo de espesor teórico
% con la ecuación de esfurzo flexionante en vigas.
Mfmax = max(sqrt(M1(xdat).^2 + datM2.^2));  Mfmax = Mfmax(1);
b = power(6*Mfmax*FS./(RcMPa*1000000),1/3);
espesor_mm = ceil(b*1000);
esfFlex_MPa = ((6*Mfmax./(espesor_mm/1000).^3))/1000000;

handles.datosFp.String = sprintf(strcat('Carga por paneles\n','Fp = ',...
    num2str(wp*L1),'[N]\n','wp = ',num2str(wp),'[N/m]'));
handles.datosFv.String = sprintf(strcat('Carga por viento\n','Fv = ',...
    num2str(wv*L1),'[N]\n','wv = ',num2str(wv),'[N/m]'));
axes(handles.grafViga);
plot([0,L1],[10,10],'k','linewidth',3); grid on; hold on;
plot([La,L1-Lb],[10,10],'ob','linewidth',3);
plot(dist,ones(1,length(dist))*10,'*r','linewidth',3);
title('Apoyos y fuerzas en la viga'); xlim([-0.5 L1+0.5]);
legend('Longitud de la viga','apoyos','fuerzas'); hold off;

axes(handles.grafMomento);
plot(xdat,M1(xdat),'r','linewidth',2); grid on; hold on;
plot(xdat,datM2,'b','linewidth',2); xlim([-0.3 L1+0.3]);
legend('Momento flector por viento','Momento flector por paneles/cargas externas');
title('Gráfica de momento flector'); hold off;

if(length(materiales) < 1)
    msgbox('No se ha seleccionado ningún material','Mensaje','help');
    set(handles.tabla,'Data',[]);
    return;
end
datos = [cellstr(materiales) num2cell(round(RcMPa)) num2cell(round(RuMPa)) num2cell(esfFlex_MPa) num2cell(espesor_mm)];
set(handles.tabla,'Data',datos);
set(handles.tabla,'visible','on');

% --- Executes on button press in mat1.
function mat1_Callback(hObject, eventdata, handles)

% --- Executes on button press in mat2.
function mat2_Callback(hObject, eventdata, handles)

% --- Executes on button press in mat3.
function mat3_Callback(hObject, eventdata, handles)

% --- Executes on button press in mat4.
function mat4_Callback(hObject, eventdata, handles)

% --- Executes on button press in mat6.
function mat6_Callback(hObject, eventdata, handles)

% --- Executes on button press in mat5.
function mat5_Callback(hObject, eventdata, handles)

% --- Executes on button press in mat8.
function mat8_Callback(hObject, eventdata, handles)

% --- Executes on button press in mat7.
function mat7_Callback(hObject, eventdata, handles)

% --- Executes on button press in mat9.
function mat9_Callback(hObject, eventdata, handles)

% --- Executes on button press in mat10.
function mat10_Callback(hObject, eventdata, handles)

% --- Executes on button press in mat11.
function mat11_Callback(hObject, eventdata, handles)

% --- Executes on button press in mat12.
function mat12_Callback(hObject, eventdata, handles)

% --- Executes on button press in mat13.
function mat13_Callback(hObject, eventdata, handles)

% --- Executes when selected object is changed in panelAnalisis.
function panelAnalisis_SelectionChangedFcn(hObject, eventdata, handles)
switch(hObject.String)
    case 'Primer Analisis'
        set(handles.panelAnalisis,'visible','off');
        set(handles.salir,'visible','off');
        set(handles.panelA1,'visible','on');
        set(handles.panelA2,'visible','off');
    case 'Segundo Analisis'
        set(handles.panelAnalisis,'visible','off');
        set(handles.salir,'visible','off');
        set(handles.panelA1,'visible','off');
        set(handles.panelA2,'visible','on');
        
        set(handles.t1_sec,'visible','off');
        set(handles.t2_sec,'visible','off');
        set(handles.t3_sec,'visible','off');
        set(handles.t4_sec,'visible','off');
        set(handles.t5_sec,'visible','off');
        set(handles.v1_sec,'visible','off');
        set(handles.v2_sec,'visible','off');
        set(handles.v3_sec,'visible','off');
        set(handles.v4_sec,'visible','off');
        set(handles.v5_sec,'visible','off');
end


% --- Executes on button press in volvA1.
function volvA1_Callback(hObject, eventdata, handles)
set(handles.panelA1,'visible','off');
set(handles.panelAnalisis,'visible','on');
set(handles.salir,'visible','on');


function distancias_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function distancias_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fuerzas_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function fuerzas_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in cargaViento.
function cargaViento_Callback(hObject, eventdata, handles)

% --- Executes on button press in bandFuerzas.
function bandFuerzas_Callback(hObject, eventdata, handles)
if(hObject.Value == 1)
    set(handles.textFzas,'visible','on');
    set(handles.textDist,'visible','on');
    set(handles.fuerzas,'visible','on');
    set(handles.distancias,'visible','on');
else
    set(handles.textFzas,'visible','off');
    set(handles.textDist,'visible','off');
    set(handles.fuerzas,'visible','off');
    set(handles.distancias,'visible','off');
end


% --- Executes on button press in calcA2.
% Cálculo de los momentos de inercia para perfiles comeriales.
function calcA2_Callback(hObject, eventdata, handles)
sec = handles.seccion.Value;
switch(sec)
    case 1  %Rectangular sólida.
        b = str2double(handles.v1_sec.String)/1000;
        h = str2double(handles.v2_sec.String)/1000;
        Ac = b*h;
        I = (b*h^3)/12;
        c = h/2;
        if( sum(isnan([b h])) > 0 )
            msgbox('Un campo de la sección está vacío o posee datos invalidos','Advertencia','warn');
            return;
        end
    case 2  %Rectangular hueca.
        b = str2double(handles.v1_sec.String)/1000;
        h = str2double(handles.v2_sec.String)/1000;
        e = str2double(handles.v3_sec.String)/1000;
        Ac = b*h - (b-e)*(h-e);
        I = (1/12)*( b*h^3 - (b-2*e)*(h-2*e)^3 );
        c = h/2;
        if( sum(isnan([b h e])) > 0 )
            msgbox('Un campo de la sección está vacío o posee datos invalidos','Advertencia','warn');
            return;
        end
    case 3  %Circular sólida.
        d = str2double(handles.v1_sec.String)/1000;
        Ac = (pi*d^2)/4;
        I = (pi*d^4)/64;
        c = d/2;
        if( isnan(d) )
            msgbox('Un campo de la sección está vacío o posee datos invalidos','Advertencia','warn');
            return;
        end
    case 4  %Circular hueca.
        D = str2double(handles.v1_sec.String)/1000;
        e = str2double(handles.v2_sec.String)/1000;
        d = 2*(0.5*D - e);
        Ac = (pi/4)*(D^2 - d^2);
        I = (pi/64)*( D^4 - d^4 );
        c = D/2;
        if( sum(isnan([D e])) > 0 )
            msgbox('Un campo de la sección está vacío o posee datos invalidos','Advertencia','warn');
            return;
        end
    case 5  %Tipo C.
        b = str2double(handles.v1_sec.String)/1000;
        h = str2double(handles.v2_sec.String)/1000;
        e = str2double(handles.v3_sec.String)/1000;
        Ic = [(h*e^3)/12 (e*(b-2*e)^3)/12 (h*e^3)/12];
        Ac = [h*e e*(b-2*e) h*e];
        yc = [e/2 b/2 b-(e/2)];
        I = inercia(Ic,Ac,yc);
        c = b/2;
        if( sum(isnan([b h e])) > 0 )
            msgbox('Un campo de la sección está vacío o posee datos invalidos','Advertencia','warn');
            return;
        end
    case 6  %Tipo I.
        b = str2double(handles.v1_sec.String)/1000;
        h = str2double(handles.v2_sec.String)/1000;
        ea = str2double(handles.v3_sec.String)/1000;
        ep = str2double(handles.v4_sec.String)/1000;
        T = str2double(handles.v5_sec.String)/1000;
        Ic = [(b*ep^3)/12 (ea*(h-2*ep)^3)/12 (b*ep^3)/12];
        Ac = [b*ep ea*(h-2*ep) b*ep];
        yc = [ep/2 h/2 h-(ep/2)];
        I = inercia(Ic,Ac,yc);
        c = h/2;
        if( sum(isnan([b h ea ep T])) > 0 )
            msgbox('Un campo de la sección está vacío o posee datos invalidos','Advertencia','warn');
            return;
        end
    case 7  %Tipo L. 
        b = str2double(handles.v1_sec.String)/1000;
        h = str2double(handles.v2_sec.String)/1000;
        e = str2double(handles.v3_sec.String)/1000;
        Ic = [(e*h^3)/12 ((b-e)*e^3)/12];
        Ac = [h*e (b-e)*e];
        yc = [h/2 e/2];
        I = inercia(Ic,Ac,yc);
        c = h/2;
        if( sum(isnan([b h e])) > 0 )
            msgbox('Un campo de la sección está vacío o posee datos invalidos','Advertencia','warn');
            return;
        end
    case 8  %Tipo T.
        b = str2double(handles.v1_sec.String)/1000;
        h = str2double(handles.v2_sec.String)/1000;
        e = str2double(handles.v3_sec.String)/1000;
        Ic = [(e*(h-e)^3)/12 (b*e^3)/12];
        Ac = [(h-e)*e b*e];
        yc = [e+((h-e)/2) e/2];
        I = inercia(Ic,Ac,yc);
        c = h/2;
        if( sum(isnan([b h e])) > 0 )
            msgbox('Un campo de la sección está vacío o posee datos invalidos','Advertencia','warn');
            return;
        end
end
mat = handles.listMat.Value;
Rc = handles.rcksi(mat)*6.8944;
Ru = handles.ruksi(mat)*6.8944;
g = 9.81;           %constante de aceleraci?n de la gravedad [m/s^2].
L1p = 1.642;        %Largo del panel (dimension mas larga) [m].
L2p = 0.994;        %Ancho del panel (dimension mas corta) [m].
mp = str2double(get(handles.mp,'String'));      %Masa de cada panel [Kg].
vmax = str2double(get(handles.vmax,'String'));  %Velocidad maxima del viento [m/s].
Ce = str2double(get(handles.Ce,'String'));      %Coeficiente eolico.
Fp = mp*g;  
pAire = 1.22;           %Densidad teórica del aire.
Aefp = L1p*L2p;         %Área efectiva de contacto del panel.
Aefv = (3*L1p)*(4*L2p); %Área efectiva de contacto del viento.
Fv = Ce*(0.5*pAire*(vmax^2)*Aefv);  %Fviento = Ce[(1/2)*p*v^2]*A.

bandFp = handles.cargaPanel.Value;
bandFv = handles.cargaViento.Value;
band2 = handles.bandFuerzas.Value;
dist = []; ff = [];
if(band2 == 1)
    dist = str2num(handles.distancias.String);
    ff = str2num(handles.fuerzas.String);
end

syms x a;
La = str2double(handles.distLa.String);	%Distancia desde el borde izquierdo de la viga al primer apoyo.
Lb = str2double(handles.distLb.String); %Distancia desde el borde derecho de la viga al segundo apoyo.
L1 = str2double(handles.Lviga.String);  %Longitud de la viga.
L2 = 4*L2p;
ms = str2double(handles.ms.String)*L1;  %Masa de la viga (soporte) [Kg].
if(isnan(ms))
    msgbox('Campo de masa por metro está vacío o posee datos inválidos','Advertencia','warn');
    return;
end

% Cálculo de reacciones y momentos en la viga considerando el peso
% y el material de la viga.
Fs = ms*g;                      %Peso de la viga.
wv = bandFv*(0.2*Fv/Aefv)*L2;
wp = bandFp*(Fp/Aefp)*L2p;
resa = abs(dot(dist-La,ff));
resb = abs(dot(dist-(L1-Lb),ff));
Ra1 = ((wv*L1)*(0.5*L1-Lb))/(L1-La-Lb);
Rb1 = ((wv*L1)*(0.5*L1-La))/(L1-La-Lb);
Ra2 = ((wp*L1 + Fs)*(0.5*L1-Lb) + resb)/(L1-La-Lb);
Rb2 = ((wp*L1 + Fs)*(0.5*L1-La) + resa)/(L1-La-Lb);
ResF = @(x) sum(double((ff.*(dist - x)).*(x>=dist)));
M1 = @(x) double(Ra1*(x - La).*(x>=La) + Rb1*(x - L1 + Lb).*(x>=(L1-Lb))...
    - wv*(0.5*x.^2));
M2 = @(x) double(Ra2*(x - La).*(x>=La) + Rb2*(x - L1 + Lb).*(x>=(L1-Lb))...
    - Fs*(x - 0.5*L1).*(x>=(0.5*L1)) - wp*(0.5*x.^2) + ResF(x));

%Despliega información en la interfaz.
handles.datosFp2.String = sprintf(strcat('Carga por paneles\n','Fp = ',...
    num2str(wp*L1),'[N]\n','wp = ',num2str(wp),'[N/m]'));
handles.datosFv2.String = sprintf(strcat('Carga por viento\n','Fv = ',...
    num2str(wv*L1),'[N]\n','wv = ',num2str(wv),'[N/m]'));
handles.datosFs2.String = sprintf(strcat('Peso de la viga\n','Wviga = ',...
    num2str(Fs),'[N]\n','masa = ',num2str(ms),'[kg]'));
Rra = sqrt(Ra1^2 + Ra2^2);  Rrb = sqrt(Rb1^2 + Rb2^2);  f = Ra2*sind(a) + Ra1*sind(90-a);  
aCrit = double(solve( diff(f,a)==0,a ));   aCrit = min(abs(aCrit));
if( bandFv*bandFp == 1 )
    handles.reacciones.String = sprintf(strcat('Reacciones en los apoyos\n',...
        'Rra = ',num2str(Rra),'[N]\n','Rrb = ',num2str(Rrb),'[N]\n','Ángulo crítico: ',num2str(aCrit,'%.2f'),'°'));
else
    handles.reacciones.String = sprintf(strcat('Reacciones en los apoyos\n',...
        'Rra = ',num2str(Rra),'[N]\n','Rrb = ',num2str(Rrb),'[N]\n'));
end

%Genera la gráfica del momento flexionante.
xdat = 0:0.001:L1;
datM = zeros(1,length(xdat));
for i=1:length(xdat)
    datM(i) = sqrt(M1(xdat(i))^2 + M2(xdat(i))^2);
end
[Mfmax,xmax] = max(datM);  Mfmax = Mfmax(1);
Sf = Mfmax*c/I;
axes(handles.grafMfA2);
plot(xdat,datM,'b','linewidth',2); grid on; hold on;
plot(xdat(xmax),Mfmax,'*r','linewidth',5);
title('Gráfica del momento flexionante producido'); xlim([-0.3 L1+0.3]); hold off;
legend('Momento flector',strcat('Momento máximo producido=',num2str(Mfmax),'[Nm]'));

datos = [num2cell(round(Rc)) num2cell(round(Ru)) num2cell(Sf) num2cell(round(Rc*1000000/Sf,2))];
set(handles.tablaA2,'Data',datos);
set(handles.tablaA2,'visible','on');
set(handles.textArea,'String',strcat('Area de la sección: ',num2str(sum(Ac)*1000000),' [mm^2]'));



% --- Executes on button press in volvA2.
function volvA2_Callback(hObject, eventdata, handles)
set(handles.panelA2,'visible','off');
set(handles.panelAnalisis,'visible','on');
set(handles.salir,'visible','on');


function v1_sec_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function v1_sec_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in listMat.
function listMat_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function listMat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in seccion.
function seccion_Callback(hObject, eventdata, handles)
set(handles.t1_sec,'visible','off');
set(handles.t2_sec,'visible','off');
set(handles.t3_sec,'visible','off');
set(handles.t4_sec,'visible','off');
set(handles.t5_sec,'visible','off');
set(handles.v1_sec,'visible','off');
set(handles.v2_sec,'visible','off');
set(handles.v3_sec,'visible','off');
set(handles.v4_sec,'visible','off');
set(handles.v5_sec,'visible','off');
sec = hObject.Value;
switch(sec)
    case 1
        handles.t1_sec.Visible = 'on';  handles.t1_sec.String = 'Largo';
        handles.t2_sec.Visible = 'on';  handles.t2_sec.String = 'Alto';
        handles.v1_sec.Visible = 'on';  handles.v2_sec.Visible = 'on';
    case 2
        handles.t1_sec.Visible = 'on';  handles.t1_sec.String = 'Largo';
        handles.t2_sec.Visible = 'on';  handles.t2_sec.String = 'Alto';
        handles.t3_sec.Visible = 'on';  handles.t3_sec.String = 'Espesor';
        handles.v1_sec.Visible = 'on';  handles.v2_sec.Visible = 'on';  handles.v3_sec.Visible = 'on';
    case 3
        handles.t1_sec.Visible = 'on';  handles.t1_sec.String = 'Diámetro';
        handles.v1_sec.Visible = 'on';
    case 4
        handles.t1_sec.Visible = 'on';  handles.t1_sec.String = 'Diámetro externo';
        handles.t2_sec.Visible = 'on';  handles.t2_sec.String = 'Espesor';
        handles.v1_sec.Visible = 'on';  handles.v2_sec.Visible = 'on';
        
    case 5
        handles.t1_sec.Visible = 'on';  handles.t1_sec.String = 'Largo';
        handles.t2_sec.Visible = 'on';  handles.t2_sec.String = 'Alto';
        handles.t3_sec.Visible = 'on';  handles.t3_sec.String = 'Espesor';
        handles.v1_sec.Visible = 'on';  handles.v2_sec.Visible = 'on';  handles.v3_sec.Visible = 'on';
    case 6
        handles.t1_sec.Visible = 'on';  handles.t1_sec.String = 'Patin';
        handles.t2_sec.Visible = 'on';  handles.t2_sec.String = 'Peralte';
        handles.t3_sec.Visible = 'on';  handles.t3_sec.String = 'Esp Alma';
        handles.t4_sec.Visible = 'on';  handles.t4_sec.String = 'Esp Patin';
        handles.t5_sec.Visible = 'on';  handles.t5_sec.String = 'Altura T';
        handles.v1_sec.Visible = 'on';  handles.v2_sec.Visible = 'on';  handles.v3_sec.Visible = 'on';  handles.v4_sec.Visible = 'on';  handles.v5_sec.Visible = 'on';
    case 7
        handles.t1_sec.Visible = 'on';  handles.t1_sec.String = 'Largo';
        handles.t2_sec.Visible = 'on';  handles.t2_sec.String = 'Alto';
        handles.t3_sec.Visible = 'on';  handles.t3_sec.String = 'Espesor';
        handles.v1_sec.Visible = 'on';  handles.v2_sec.Visible = 'on';  handles.v3_sec.Visible = 'on';
    case 8
        handles.t1_sec.Visible = 'on';  handles.t1_sec.String = 'Largo';
        handles.t2_sec.Visible = 'on';  handles.t2_sec.String = 'Alto';
        handles.t3_sec.Visible = 'on';  handles.t3_sec.String = 'Espesor';
        handles.v1_sec.Visible = 'on';  handles.v2_sec.Visible = 'on';  handles.v3_sec.Visible = 'on';
end


% --- Executes during object creation, after setting all properties.
function seccion_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function v2_sec_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function v2_sec_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function v3_sec_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function v3_sec_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function v4_sec_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function v4_sec_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function v5_sec_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function v5_sec_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ms_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function ms_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in cargaPanel.
function cargaPanel_Callback(hObject, eventdata, handles)

function distLa_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function distLa_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function distLb_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function distLb_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Lviga_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function Lviga_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
