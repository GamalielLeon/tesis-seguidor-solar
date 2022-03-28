/* USER CODE BEGIN Header */
/**
  * Seguidor Solar Industrial de Uso Ligero, Generador de 
  * Energia Eléctrica mediante un Sistema Fotovoltaico
  *
  * Este software implementa el protocolo de comunicaciones 
  * EtherNet/IP para indicarle a la NUCLEO-STM32F767ZI el modo 
  * de operación por el cual obtendrá los ángulos de posición 
  * del robot, así­ como el control de su movimiento.
  *
  **************************************************************
  */
/* USER CODE END Header */

/* Librerías --------------------------------------------------*/
#include "main.h"
#include "lwip.h"
#include "stdio.h"
#include "stdlib.h"
#include <math.h>
#include "lwip/apps/httpd.h"// incluyendo httpd.h [- HTTPd #1 -]
#include <String.h>	
#include <Stdbool.h>	// utilizar variables booleanas para SSI

#ifndef PI
#define PI 3.14159265358979
#endif

#define PI2 6.28318530717959    // 2*PI
#define PIM 1.57079632679490    // PI/2

/* Variables privadas y globales -----------------------------*/
UART_HandleTypeDef huart3;
PCD_HandleTypeDef hpcd_USB_OTG_FS;

bool onON = false;
bool autoON = false;
bool manuON = false;
bool bcON = false;
bool delayON=false;
bool alfaON=false;
bool gamaON=false;
bool fhsON=false;
char *delays;
char *alfas;
char *gamas;
char *FHS;
unsigned short int D=0;

/* USER CODE BEGIN PV */
// simplemente declarando la función para el compilador
const char* LedCGIhandler(int iIndex,int iNumParams,char *pcParam[],
                            char *pcValue[]);
// en nuestro archivo SHTML <form method="get" action="/leds.cgi">
const tCGI LedCGI = { "/leds.cgi", LedCGIhandler };
// [= CGI #4 =]
tCGI theCGItable[1];
// simplemente declarando la función del controlador SSI
u16_t mySSIHandler(int iIndex, char *pcInsert, int iInsertLen);
#define numSSItags 8
char const *theSSItags[numSSItags]={"tag1","tag2","tag3","tag4",
                                    "tag5","tag6","tag7","tag8"};
/* USER CODE END PV */

/* Prototipos de funciones privadas. ---------------------------*/
void SystemClock_Config(void);
static void MX_GPIO_Init(void);
static void MX_USART3_UART_Init(void);
static void MX_USB_OTG_FS_PCD_Init(void);

bool tiempo(float);			    //Determina si es de día o de noche
float presion();				//Sensor de presión
float temperatura();			//Sensor de temperatura
unsigned short int humedad();	//Mide el umbral de humedad
void GPS(float punto[2]);		//Sensor GPS
void fecha(unsigned short int datos[6]);//Registra la fecha y hora
//Tiempo de espera para realizar el siguiente movimiento
void espera(unsigned short int,int);	
//Modo de operación Automático
void Automatico(float,unsigned short int *,float angle[2]);
//Función que calcula los ángulos de posicionamiento solar				
void ASS(float,unsigned short int *,float,float *,float,float,
            float ang[2]);

/* Código de usuario privado -------------------------------------*/
/* USER CODE BEGIN 0 */

// la función real para manejar CGI [= CGI #5 =]
const char* LedCGIhandler(int iIndex, int iNumParams, char *pcParam[],
                            char *pcValue[]) {
  unsigned short int estado=0;
  int delay=atoi(delays);
  float hora=0.0; //hora entera más parte decimal (minutos y segundos)
  float alfa=0.0;	//Ángulo de Elevación
  float gama=0.0;	//Ángulo Azimutal
  float angulos[2]={0.0,0.0};

  if (iIndex == 0) {
	HAL_GPIO_WritePin(LD1_GPIO_Port, LD1_Pin, GPIO_PIN_RESET);
	HAL_GPIO_WritePin(LD2_GPIO_Port, LD2_Pin, GPIO_PIN_RESET);
	HAL_GPIO_WritePin(LD3_GPIO_Port, LD3_Pin, GPIO_PIN_RESET);
	onON = false;
	autoON = false;
	manuON = false;
	bcON = true;
	delayON=false;
	alfaON=true;
	gamaON=true;
	fhsON=false;
  }
  if (strcmp(pcValue[0], "0") == 0) {
	HAL_GPIO_WritePin(LD1_GPIO_Port, LD1_Pin, GPIO_PIN_SET);
	onON = true;
	if (strcmp(pcValue[1], "1") == 0) {
	  HAL_GPIO_WritePin(LD2_GPIO_Port, LD2_Pin, GPIO_PIN_SET);
	  autoON = true;
	  manuON = false;
	  bcON = false;
	  delayON=true;
	  alfaON=true;
	  gamaON=true;
	  fhsON=true;
	  D=1;
	}
	else if (strcmp(pcValue[1], "2") == 0) {
	  HAL_GPIO_WritePin(LD3_GPIO_Port, LD3_Pin, GPIO_PIN_SET);
	  autoON = false;
	  manuON = true;
	  bcON = false;
	  delayON=false;
	  alfaON=true;
	  gamaON=true;
	  fhsON=false;
	  D=2;
	}
	else if (strcmp(pcValue[1], "3") == 0) {
	  HAL_GPIO_WritePin(LD3_GPIO_Port, LD3_Pin, GPIO_PIN_SET);
	  HAL_GPIO_WritePin(LD2_GPIO_Port, LD2_Pin, GPIO_PIN_SET);
	  autoON = false;
	  manuON = false;
	  bcON = true;
	  delayON=false;
	  alfaON=true;
	  gamaON=true;
	  fhsON=true;
	  D=0;
	}

	char lecturaFH[30];
	uint16_t leidoFH;
	char aux[16];
	char *p;
	p=pcValue[5];
	for(int j=0; j<13; j++)
	  aux[j]=*(p+j);
	aux[13]=':';
	aux[14]=*(p+16);
	aux[15]=*(p+17);
	strcpy(pcValue[5],"");
	strncat(pcValue[5],aux,16);
	FHS = pcValue[5];

	unsigned short int dia[6],H;
	fecha(dia);     //Se obtiene un arreglo con los datos leídos
    //Se calcula la hora como dato flotante
	hora=(float)dia[3]+(float)dia[4]/60+(float)dia[5]/3600;	
	H=humedad();
	estado=3*H+D;                       //Cálculo del estado
    //Si el estado no es ni 1 ni 2, se pasará al modo de Bajo Consumo
	if(estado==1 || estado==2){		                
      //Si el estado es 1 y es de día, se pasará al modo Automático        
	  if(estado==1 && tiempo(hora)){	                    

		char lectura[30];
		uint16_t leido;
		delays=pcValue[2];

		float angle[2];
        //Se generan los ángulos calculados por el modo Automático
		Automatico(hora,dia,angle);
		alfa=*angle;
		gama=*(angle+1);
		//Los ángulos calculados se guardan en el arreglo global 
        //de angulos[]
		angulos[0]=alfa;
		angulos[1]=gama;
		sprintf(alfas, "%.3f", alfa);
		sprintf(gamas, "%.3f", gama);
		pcValue[4] = gamas;
        //Se aplica el control de posición con los ángulos calculados
		//control(angulos);	
        //Se hace esperar al robot para la siguiente iteración de la rutina
		//espera(D,delay);	
	  }
	  else{
		HAL_Delay(100);	    //Tiempo de espera para el uC

		char lecturaA[30];
		uint16_t leidoA;
		alfas=pcValue[3];
		char lecturaZ[30];
		uint16_t leidoZ;
		gamas=pcValue[4];

        //Se toma el ángulo de Elevación desde la Interfaz
		alfa=atof(alfas);
        //Se toma el ángulo Azimutal desde la Interfaz
		gama=atof(gamas);
		angulos[0]=alfa;
		angulos[1]=gama;
        //Se aplica el control de posición con los ángulos calculados
		//control(angulos);	
	  }
	}
	else{
	  HAL_Delay(100);				//Tiempo de espera para el uC
	  alfa=0.0;		    //Forza al ángulo de Elevación igual a 0°
      //Conserva el último ángulo Azimutal calculado
	  sprintf(alfas, "%.3f", alfa);	
	  gamas = pcValue[4];
	  angulos[0]=alfa;
	  angulos[1]=atof(gamas);
      //Se aplica el control de posición con los ángulos calculados
	  //control(angulos);
	}
  }
  else{
	autoON = false;
	manuON = false;
	bcON = false;
	delayON=false;
	alfaON=false;
	gamaON=false;
	fhsON=false;
  }
  return "/index.shtml";// la extensión .shtml para que SSI funcione
}

// función para inicializar CGI [= CGI #6 =]
void myCGIinit(void) {
	//agregar control LED CGI a la mesa
	theCGItable[0] = LedCGI;
	//dar la tabla al servidor HTTP
	http_set_cgi_handlers(theCGItable, 1);
} // END [= CGI #6 =]

//Función que evalúa si es de dí­a o de noche en función de la 
//hora en formato decimal
bool tiempo(float h){
  //Si la hora se encuentra entre las horas del amanecer y anochecer,
  //devuelve "verdadero"
  if(h>=7.0 && h<=20.0)	
	return true;
  else				//Si no, será de noche y devolverá "falso"
	return false;
}

//Función que mide la presión atmosférica en atmósferas [atm] de 0 a 1
float presion(){
	float p=0.761;
	return p;
}

//Función que mide la temperatura del ambiente en grados Celsius [°C]
float temperatura(){
	float t=25.0;
	return t;
}

//Función que devuelve un número entero corto sin signo si la 
//humedad es alta
unsigned short int humedad(){
	float lectura=50.0;		//Lectura del sensor de humedad
	if(lectura>=70.0)		//Umbral de humedad en la CDMX
		return 1;
	else
		return 0;
}

//Función que utiliza un GPS para obtener las coordenadas del 
//robot en el sistema estándar WSG84
void GPS(float punto[2]){
  //Crea el arreglo en el que se guardarán las coordenadas
  float longitud=-1.730060623,latitud=0.3405426746;	
  //Guarda las lecturas del GPS en el arreglo
  punto[0]=longitud; punto[1]=latitud; 
  //Devuelve el arreglo con las coordenadas obtenidas 
  return;								
}

//Función que obtiene los datos de la fecha y hora en un arreglo
void fecha(unsigned short int datos[6]){
  char anio[4]={0},mes[2]={0},day[2]={0},hour[2]={0},min[2]={0};
  //Toma los datos de la etiqueta "datetime-local" de la Interfaz
  anio[0]=*FHS;
  anio[1]=*(FHS+1);
  anio[2]=*(FHS+2);
  anio[3]=*(FHS+3);
  mes[0]=*(FHS+5);
  mes[1]=*(FHS+6);
  day[0]=*(FHS+8);
  day[1]=*(FHS+9);
  hour[0]=*(FHS+11);
  hour[1]=*(FHS+12);
  min[0]=*(FHS+14);
  min[1]=*(FHS+15);
  datos[0]=atoi(day);	//Día actual
  datos[1]=atoi(mes);	//Mes actual
  datos[2]=atoi(anio);	//Año actual
  datos[3]=atoi(hour);	//Hora de la lectura
  datos[4]=atoi(min);	//Minutos de la lectura
  datos[5]=0;			//Segundos de la lectura
  return;	//Devuelve el arreglo con los datos requeridos
}

//Función que hace esperar al seguidor para su siguiente acción
void espera(unsigned short int b,int c){
  int fin=c;//*60;
  //El contador se detendrá hasta que llegue al tiempo requerido 
  //desde la Interfaz
  for(int i=1;i<=fin;i++){	
    //Si se apaga el robot o si se cambia de modo de operación, 
    //la espera se interrumpe
	if(onON==false || b!=1)	
	  i=fin;
	else
      //Es posible cambiar el modo de operación durante la espera
	  b=decision();	
	HAL_Delay(1000);
  }
}

//Modo de operación Automático. Devuelve los ángulos calculados 
//en un arreglo flotante.
void Automatico(float h,unsigned short int *datos,float angle[2]){
  HAL_Delay(100);//Tiempo de espera para el uC
  float UT=5.0;	 //Tiempo Universal Coordinado (UTC) para la CDMX
  //Diferencia entre el Tiempo Terrestre y el Tiempo Universal
  float Dt=96.4+0.567*((float)datos[2]-2061.0);	
  float GTM=h+UT;	//Hora desde el Meridiano de Greenwwich
  float punto[2];	//Arreglo que devuelve longitud y latitud
  GPS(punto);
  //Cálculo de los ángulos por medio del Algoritmo de Seguimiento 
  //Solar en un arreglo, en función de los parámetros requeridos
  //Se guardan los ángulos calculados en el arreglo contenedor
  float ang[2];	
  ASS(GTM,datos,Dt,punto,presion(),temperatura(),ang);
  angle[0]=*ang;
  angle[1]=*(ang+1);
  return;			//Devuelve el arreglo con los ángulos calculados
}

//Algoritmo de Seguimiento Solar Grena 2012 #5
void ASS(float h,unsigned short int *dia,float dt,float *gps,
            float p,float T,float ang[2]){
  // input data:
  double UT;
  int Day;
  int Month;
  int Year;
  double Dt;
  double Longitude;
  double Latitude;
  double Pressure;
  double Temperature;
  //output data
  double RightAscension;
  double Declination;
  double HourAngle;
  double Zenith;
  double Azimuth;
  //auxiliary
  double t, te, wte, s1, c1, s2, c2, s3, c3, sp, cp, sd, cd, sH, cH;
  double se0, ep, De, lambda, epsi, sl, cl, se, ce, L, nu, Dlam;
  int yt, mt;

  float conv=180/PI;	//Conversor de radianes a grados
  Day = dia[0];		    //Se carga el día del movimiento
  Month = dia[1];		//Se carga el mes del movimiento
  Year = dia[2];		//Se carga el año del movimiento
  Dt = dt;			    //Se carga la diferencia de horarios
  Longitude = gps[0];	//Se carga la longitud geográfica del robot
  Latitude = gps[1];	//Se carga la latitud geográfica del robot
  Pressure = p;		    //Se carga la presión medida
  Temperature = T;	    //Se carga la temperatura medida
  UT = h;				//Se carga la hora en el sistema GTM

  if (Month <= 2) {
    mt = Month + 12;
    yt = Year - 1;
  } else {
	mt = Month;
    yt = Year;
  }

  t=(double)((int)(365.25*(double)(yt-2000))+(int)(30.6001*(double)(mt+1))
    -(int)(0.01*(double)yt)+Day)+0.0416667*UT-21958.0;
  te = t + 1.1574e-5*Dt;

  wte = 0.0172019715*te;

  s1 = sin(wte);
  c1 = cos(wte);
  s2 = 2.0*s1*c1;
  c2 = (c1+s1)*(c1-s1);
  s3 = s2*c1 + c2*s1;
  c3 = c2*c1 - s2*s1;

  L=1.7527901+1.7202792159e-2*te+3.33024e-2*s1-2.0582e-3*c1+
    3.512e-4*s2-4.07e-5*c2+5.2e-6*s3-9e-7*c3
    -8.23e-5*s1*sin(2.92e-5*te)+1.27e-5*sin(1.49e-3*te-2.337)    
    +1.21e-5*sin(4.31e-3*te+3.065)+2.33e-5*sin(1.076e-2*te-1.533)
    +3.49e-5*sin(1.575e-2*te-2.358)+2.67e-5*sin(2.152e-2*te+0.074)
    +1.28e-5*sin(3.152e-2*te+1.547)+3.14e-5*sin(2.1277e-1*te-0.488);

  nu = 9.282e-4*te - 0.8;
  Dlam = 8.34e-5*sin(nu);
  lambda = L + PI + Dlam;

  epsi = 4.089567e-1 - 6.19e-9*te + 4.46e-5*cos(nu);

  sl = sin(lambda);
  cl = cos(lambda);
  se = sin(epsi);
  ce = sqrt(1-se*se);

  RightAscension = atan2(sl*ce, cl);
  if (RightAscension < 0.0)
    RightAscension += PI2;

  Declination = asin(sl*se);

  HourAngle=1.7528311+6.300388099*t+Longitude-RightAscension+0.92*Dlam;
  HourAngle = fmod(HourAngle + PI, PI2) - PI;
  if (HourAngle < -PI) HourAngle += PI2;

  sp = sin(Latitude);
  cp = sqrt((1-sp*sp));
  sd = sin(Declination);
  cd = sqrt(1-sd*sd);
  sH = sin(HourAngle);
  cH = cos(HourAngle);
  se0 = sp*sd + cp*cd*cH;
  ep = asin(se0) - 4.26e-5*sqrt(1.0-se0*se0);
  Azimuth = atan2(sH, cH*sp - sd*cp/cd);

  if (ep > 0.0)
    De=(0.08422*Pressure)/((273.0+Temperature)*tan(ep+0.003138/(ep+0.08919)));
  else
    De = 0.0;

  Zenith = PIM - ep - De;

  //Guarda el ángulo de Elevación calculado
  ang[0]=(float)(90.0-Zenith*conv);	
  //Guarda el ángulo Azimutal calculado
  ang[1]=(float)(Azimuth*conv);		
  return;   //Devuelve el arreglo con los ángulos calculados
}

// la función real para SSI [* SSI #4 *]
u16_t mySSIHandler(int iIndex, char *pcInsert, int iInsertLen){
  if (iIndex == 0) {
	if (onON == false) {
	  char myStr1[]=
        "<input value=\"0\" name=\"on\" type=\"checkbox\">";
	  strcpy(pcInsert, myStr1);
	  return strlen(myStr1);
	}
	else if (onON == true) {
	  char myStr1[]=
        "<input value=\"0\" name=\"on\" type=\"checkbox\" checked>";
	  strcpy(pcInsert, myStr1);
	  return strlen(myStr1);
	}
  }
  else if (iIndex == 1){
	if (autoON == false) {
	  char myStr2[]=
        "<input value=\"1\" name=\"led\" type=\"radio\">";
	  strcpy(pcInsert, myStr2);
	  return strlen(myStr2);
	}
	else if (autoON == true) {
	  char myStr2[]=
        "<input value=\"1\" name=\"led\" type=\"radio\" checked>";
	  strcpy(pcInsert, myStr2);
	  return strlen(myStr2);
	  D=1;
	}
  }
  else if (iIndex == 2){
	if (manuON == false) {
	  char myStr3[]=
        "<input value=\"2\" name=\"led\" type=\"radio\">";
	  strcpy(pcInsert, myStr3);
	  return strlen(myStr3);
	}
	else if (manuON == true) {
	  char myStr3[]=
        "<input value=\"2\" name=\"led\" type=\"radio\" checked>";
	  strcpy(pcInsert, myStr3);
	  return strlen(myStr3);
	  D=2;
	}
  }
  else if (iIndex == 3){
	if (bcON == false) {
	  char myStr4[]=
        "<input value=\"3\" name=\"led\" type=\"radio\">";
	  strcpy(pcInsert, myStr4);
	  return strlen(myStr4);
	}
	else if (bcON == true) {
	  char myStr4[]=
        "<input value=\"3\" name=\"led\" type=\"radio\" checked>";
	  strcpy(pcInsert, myStr4);
	  return strlen(myStr4);
	  D=0;
	}
  }
  else if (iIndex == 4){
	if (delayON == false) {
	  char myStr5[]="<input name=\"delay\" type=\"number\">";
	  strcpy(pcInsert, myStr5);
	  return strlen(myStr5);
	}
	else if (delayON == true) {
	  char myStr5[]=
        "<input name=\"delay\" type=\"number\" value=\"";
	  strcat(myStr5,delays);
	  strcat(myStr5,"\">");
	  strcpy(pcInsert, myStr5);
	  return strlen(myStr5);
	}
  }
  else if (iIndex == 5){
	if (alfaON == false) {
	  char myStr6[]=
        "<input name=\"alfa\" type=\"number\" min=\"0\" max=\"90\">";
	  strcpy(pcInsert, myStr6);
	  return strlen(myStr6);
	}
	else if (alfaON == true) {
	  char myStr6[]="<input name=\"alfa\" type=\"number\" value=\"";
	  strcat(myStr6,alfas);
	  strcat(myStr6,"\">");
	  strcpy(pcInsert, myStr6);
	  return strlen(myStr6);
	}
  }
  else if (iIndex == 6){
	if (gamaON == false) {
	  char myStr7[]=
        "<input name=\"gama\" type=\"number\" min=\"-180\" max=\"180\">";
	  strcpy(pcInsert, myStr7);
	  return strlen(myStr7);
	}
	else if (gamaON == true) {
	  char myStr7[]="<input name=\"gama\" type=\"number\" value=\"";
	  strcat(myStr7,gamas);
	  strcat(myStr7,"\">");
	  strcpy(pcInsert, myStr7);
	  return strlen(myStr7);
	}
  }
  else if (iIndex == 7){
	if (fhsON == false) {
	  char myStr8[]=
        "<input type=\"datetime-local\" name=\"movido\" value=\"2020-07-14T00:00\">";
	  strcpy(pcInsert, myStr8);
	  return strlen(myStr8);
	}
	else if (fhsON == true) {
	  char myStr8[]=
        "<input type=\"datetime-local\" name=\"movido\" value=\"";
	  strcat(myStr8,FHS);
	  strcat(myStr8,"\">");
	  strcpy(pcInsert, myStr8);
	  return strlen(myStr8);
	}
  }
  return 0;
}

// función para inicializar SSI [* SSI #5 *]
void mySSIinit(void) {
  http_set_ssi_handler(mySSIHandler, 
    (char const**) theSSItags,numSSItags);
}

/* USER CODE END 0 */

int main(void){
  /* Configuración de MCU------------------------------------*/
  /* Reinicio de todos los periféricos, inicializa la interfaz Flash y el Systick. */
  HAL_Init();               
  SystemClock_Config();     /* Configura el reloj del sistema */
  /* Inicializa todos los periféricos configurados */
  MX_GPIO_Init();
  MX_USART3_UART_Init();
  MX_USB_OTG_FS_PCD_Init();
  MX_LWIP_Init();

  /* USER CODE BEGIN 2 */
  httpd_init();	// inicializando el HTTPd [-HTTPd #2-]
  myCGIinit();	// inicializando el CGI  [= CGI #7 =]
  mySSIinit();	// inicializando el SSI [* SSI #6 *]
  /* USER CODE END 2 */

  /* Bucle infinito */
  while (1){
	// comenzando el proceso LWIP [-HTTPd #3-]
	MX_LWIP_Process();
  }
}

void SystemClock_Config(void){
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};
  RCC_PeriphCLKInitTypeDef PeriphClkInitStruct = {0};
  /** Configure LSE Drive Capability */
  HAL_PWR_EnableBkUpAccess();
  /** Configure the main internal regulator output voltage */
  __HAL_RCC_PWR_CLK_ENABLE();
  __HAL_PWR_VOLTAGESCALING_CONFIG(PWR_REGULATOR_VOLTAGE_SCALE1);
  /** Initializes the CPU, AHB and APB busses clocks */
  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSE;
  RCC_OscInitStruct.HSEState = RCC_HSE_BYPASS;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_ON;
  RCC_OscInitStruct.PLL.PLLSource = RCC_PLLSOURCE_HSE;
  RCC_OscInitStruct.PLL.PLLM = 4;
  RCC_OscInitStruct.PLL.PLLN = 216;
  RCC_OscInitStruct.PLL.PLLP = RCC_PLLP_DIV2;
  RCC_OscInitStruct.PLL.PLLQ = 9;
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK){
    Error_Handler();
  }
  /** Activate the Over-Drive mode */
  if (HAL_PWREx_EnableOverDrive() != HAL_OK){
    Error_Handler();
  }
  /** Initializes the CPU, AHB and APB busses clocks */
  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK|RCC_CLOCKTYPE_SYSCLK|RCC_CLOCKTYPE_PCLK1|RCC_CLOCKTYPE_PCLK2;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_PLLCLK;
  RCC_ClkInitStruct.AHBCLKDivider = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV4;
  RCC_ClkInitStruct.APB2CLKDivider = RCC_HCLK_DIV2;
  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_7) != HAL_OK){
    Error_Handler();
  }
  PeriphClkInitStruct.PeriphClockSelection = RCC_PERIPHCLK_USART3|RCC_PERIPHCLK_CLK48;
  PeriphClkInitStruct.Usart3ClockSelection = RCC_USART3CLKSOURCE_PCLK1;
  PeriphClkInitStruct.Clk48ClockSelection = RCC_CLK48SOURCE_PLL;
  if (HAL_RCCEx_PeriphCLKConfig(&PeriphClkInitStruct) != HAL_OK){
    Error_Handler();
  }
}

static void MX_USART3_UART_Init(void){
  huart3.Instance = USART3;
  huart3.Init.BaudRate = 9600;
  huart3.Init.WordLength = UART_WORDLENGTH_8B;
  huart3.Init.StopBits = UART_STOPBITS_1;
  huart3.Init.Parity = UART_PARITY_NONE;
  huart3.Init.Mode = UART_MODE_TX_RX;
  huart3.Init.HwFlowCtl = UART_HWCONTROL_NONE;
  huart3.Init.OverSampling = UART_OVERSAMPLING_16;
  huart3.Init.OneBitSampling = UART_ONE_BIT_SAMPLE_DISABLE;
  huart3.AdvancedInit.AdvFeatureInit = UART_ADVFEATURE_NO_INIT;
  if (HAL_UART_Init(&huart3) != HAL_OK){
    Error_Handler();
  }
}

static void MX_USB_OTG_FS_PCD_Init(void){
  hpcd_USB_OTG_FS.Instance = USB_OTG_FS;
  hpcd_USB_OTG_FS.Init.dev_endpoints = 6;
  hpcd_USB_OTG_FS.Init.speed = PCD_SPEED_FULL;
  hpcd_USB_OTG_FS.Init.dma_enable = DISABLE;
  hpcd_USB_OTG_FS.Init.phy_itface = PCD_PHY_EMBEDDED;
  hpcd_USB_OTG_FS.Init.Sof_enable = ENABLE;
  hpcd_USB_OTG_FS.Init.low_power_enable = DISABLE;
  hpcd_USB_OTG_FS.Init.lpm_enable = DISABLE;
  hpcd_USB_OTG_FS.Init.vbus_sensing_enable = ENABLE;
  hpcd_USB_OTG_FS.Init.use_dedicated_ep1 = DISABLE;
  if (HAL_PCD_Init(&hpcd_USB_OTG_FS) != HAL_OK){
    Error_Handler();
  }
}

static void MX_GPIO_Init(void){
  GPIO_InitTypeDef GPIO_InitStruct = {0};
  /* GPIO Ports Clock Enable */
  __HAL_RCC_GPIOC_CLK_ENABLE();
  __HAL_RCC_GPIOH_CLK_ENABLE();
  __HAL_RCC_GPIOA_CLK_ENABLE();
  __HAL_RCC_GPIOB_CLK_ENABLE();
  __HAL_RCC_GPIOD_CLK_ENABLE();
  __HAL_RCC_GPIOG_CLK_ENABLE();
  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(GPIOB, LD1_Pin|LD3_Pin|LD2_Pin, GPIO_PIN_RESET);
  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(USB_PowerSwitchOn_GPIO_Port, USB_PowerSwitchOn_Pin, GPIO_PIN_RESET);
  /*Configure GPIO pin : USER_Btn_Pin */
  GPIO_InitStruct.Pin = USER_Btn_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_IT_RISING;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  HAL_GPIO_Init(USER_Btn_GPIO_Port, &GPIO_InitStruct);
  /*Configure GPIO pins : LD1_Pin LD3_Pin LD2_Pin */
  GPIO_InitStruct.Pin = LD1_Pin|LD3_Pin|LD2_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_PULLUP;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(GPIOB, &GPIO_InitStruct);
  /*Configure GPIO pin : USB_PowerSwitchOn_Pin */
  GPIO_InitStruct.Pin = USB_PowerSwitchOn_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(USB_PowerSwitchOn_GPIO_Port, &GPIO_InitStruct);
  /*Configure GPIO pin : USB_OverCurrent_Pin */
  GPIO_InitStruct.Pin = USB_OverCurrent_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_INPUT;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  HAL_GPIO_Init(USB_OverCurrent_GPIO_Port, &GPIO_InitStruct);
}

void Error_Handler(void){
}

#ifdef  USE_FULL_ASSERT
void assert_failed(uint8_t *file, uint32_t line){
}
#endif

/************************ (C) COPYRIGHT STMicroelectronics *****END OF FILE****/
