//programma finalizzato a trovare il transiente del sistema una volta assegnati i valori di periodo e ampiezze per le tre coordinate

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define COMPATIBILITA 0.01 //questa define stabilisce con quanta precisione voglio che coincidano i valori delle ampiezze e dei periodi con i rispettivi valori stazionari. Poiché in precedenza mi sono stimato T e i massimi con una precisione dell' 1%, mantengo questa precisione per stabilire quando terminerà il transiente

//////////////////////////////

struct vettore_tri
{
  double x;
  double y;
  double z;
};

struct spazio_fasi
{
  struct vettore_tri r;
  struct vettore_tri v;
};

///////////////////////////

struct vettore_tri inizializzazione(double x_0, double y_0, double z_0);
struct vettore_tri derivata (struct vettore_tri r, double a, double b, double c);
struct vettore_tri Runge_Kutta_2(struct spazio_fasi xEv, double dt, double a, double b, double c);
struct spazio_fasi uguaglianza_spaziofasi(struct spazio_fasi old);
double ricerca_zeri(double v_z_old, double v_z_new, double dt, int i);


/*********************************************************/

int main (int argc, char *argv[])
{
  struct spazio_fasi xEv_old, xEv_new;
  int n_passi, i = 0, n_zeri_v_z = 0, n_zeri_v_x = 0, n_zeri_v_y = 0, j_z = 0, k_z = 0, j_x = 0, k_x = 0, j_y = 0, k_y = 0;
  double T_max, dt, x_0, y_0, z_0, a, b, c, T, A_x, A_y, A_z, transiente, periodo_gen_z, periodo_gen_x, periodo_gen_y;
  double TEMPI_ZERI_v_z[3], TEMPI_ZERI_v_x[3], TEMPI_ZERI_v_y[3];
  FILE *pFdati;

  //controllo che sia stato inserito il numero  corretto di valori da linea di comando
  if(argc!=9)
    {
      fprintf(stderr, "#Si inseriscano i valori delle seguenti grandezze: %s, a, b, c, x0=x(0), y0=y(0), z0=z(0), delta t,  Tmax di integrazione\n", argv[0]);
      exit(1);
    }

  
  //inizializzo i valori delle condizioni iniziali e il numero di passi
  pFdati = fopen("ampiezze_periodo.dat", "r");
  fscanf(pFdati, "%lf %lf %lf %lf %lf", &dt, &A_x, &A_y, &A_z, &T);
  fclose(pFdati);
  a = atof(argv[1]);
  b = atof(argv[2]);
  c = atof(argv[3]);
  x_0 = atof(argv[4]);
  y_0 = atof(argv[5]);
  z_0 = atof(argv[6]);
  dt = atof(argv[7]);
  T_max = atof(argv[8]);
  n_passi = (int)(T_max/dt + 0.5);
  xEv_old.r = inizializzazione(x_0, y_0, z_0);
  xEv_old.v = derivata(xEv_old.r, a, b, c);
  xEv_new = uguaglianza_spaziofasi(xEv_old);

  //Stampo le condizioni iniziali
  printf("#Condizioni iniziali: a = %f, b = %f, c = %f, dt = %f, t_max = %f, x_0 = %f, y_0 = %f, z_0 = %f\n", a, b, c, dt, T_max, x_0, y_0, z_0);
  
  //implemento l'integrazione e scrivo i valori delle coordinate e delle loro velocità ad ogni istante
  for(i=1; i<=n_passi; i++) 
    {
      xEv_new.r = Runge_Kutta_2(xEv_old, dt, a, b, c);
      xEv_new.v = derivata(xEv_new.r, a, b, c);

      /********************* CALCOLO DEL TRANSIENTE  ***********************/ 
      if(xEv_old.v.z * xEv_new.v.z < 0) //ogni volta che questa condizione è soddisfatta si trova uno zero della velocità dunque un massimo o minimo della z
	{
	  n_zeri_v_z ++; //si incrementa un contatore
	  TEMPI_ZERI_v_z[n_zeri_v_z - 1] = ricerca_zeri(xEv_old.v.z, xEv_new.v.z, dt, i); //si assegna all'elemento n-1 il valore del tempo in cui si trova l' n-esimo zero della velocità       
	}
      if(n_zeri_v_z == 3) //non appena si trovano 3 zeri si calcola la distanza temporale tra il primo e il terzo e la si salva in periodo_gen_z
	{
	  periodo_gen_z = TEMPI_ZERI_v_z[2] - TEMPI_ZERI_v_z[0];
	  if(fabs((periodo_gen_z - T)/T) < COMPATIBILITA) //se il periodo appena trovato è compatibile con il periodo stazionario T nella soglia impostata tramite #define, si imposta un contatore = 1
	    j_z = 1;
	  n_zeri_v_z = 0; //si riinizializza il contatore degli zeri così da poter continuare il processo
	}
      if(fabs((xEv_old.r.z - A_z)/A_z) < COMPATIBILITA) //se l'ampiezza di z è compatibile con l'ampiezza stazionaria di z nella soglia impostata tramite #define, si imposta un nuovo contatore = 1
	{
	  k_z = 1;
	}
      /*********** si ripete lo stesso ragionamento per le x e per le y *************/  
      if(xEv_old.v.x * xEv_new.v.x < 0)
	{
	  n_zeri_v_x ++;
	  TEMPI_ZERI_v_x[n_zeri_v_x - 1] = ricerca_zeri(xEv_old.v.x, xEv_new.v.x, dt, i);
	}
      if(n_zeri_v_x == 3)
	{
	  periodo_gen_x = TEMPI_ZERI_v_x[2] - TEMPI_ZERI_v_x[0];
	  if(fabs((periodo_gen_x - T)/T) < COMPATIBILITA)
	    j_x = 1;
	  n_zeri_v_x = 0;
	}
      if(fabs((xEv_old.r.x - A_x)/A_x) < COMPATIBILITA)
	{
	  k_x = 1;
	}
      if(xEv_old.v.y * xEv_new.v.y < 0)
	{
	  n_zeri_v_y ++;
	  TEMPI_ZERI_v_y[n_zeri_v_y - 1] = ricerca_zeri(xEv_old.v.y, xEv_new.v.y, dt, i);
	}
      if(n_zeri_v_y == 3)
	{
	  periodo_gen_y = TEMPI_ZERI_v_y[2] - TEMPI_ZERI_v_y[0];
	  if(fabs((periodo_gen_y - T)/T) < COMPATIBILITA)
	    j_y = 1;
	  n_zeri_v_y = 0;
	}
      if(fabs((xEv_old.r.y - A_y)/A_y) < COMPATIBILITA)
	{
	  k_y = 1;
	}      
      if(j_z == 1 && k_z == 1 && j_x == 1 && k_x == 1 && j_y == 1 && k_y == 1) //se i periodi e le ampiezze di tutte le varibili si sono stabilizzati entro l'1% si è usciti dal transiente.
	{
	  transiente = i*dt; //si salva il valore di tempo in cui termina il transiente
	  break; //per non fare calcoli inutili, trovato il transiente si termina volontariamente il ciclo di integrazione
	}
      
      /*************** FINE CALCOLO DEL TRANSIENTE  *********************/
     
      xEv_old = uguaglianza_spaziofasi(xEv_new);
    }

  printf("Il transiente vale %.3f secondi.\n", transiente);
  
  return EXIT_SUCCESS;
}

/*********************************************************/

//si implementa una funzione che mi inizializza il vettore (x, y, z) con le condizioni iniziali 
struct vettore_tri inizializzazione(double x_0, double y_0, double z_0) 
{
  struct vettore_tri r;
  
  r.x = x_0;
  r.y = y_0;
  r.z = z_0;
  
  return r;
}

//si implementa una funzione che calcola la derivata prima del vettore (x, y, z)
struct vettore_tri derivata (struct vettore_tri r, double a, double b, double c) 
{
  struct vettore_tri v;
  
  v.x = -r.y -r.z;
  v.y = r.x + a*r.y;
  v.z = b + (r.x - c)*r.z;
  
  return v;
}

//si implementa il metodo di integrazione di Runge-Kutta al secondo ordine 
struct vettore_tri Runge_Kutta_2(struct spazio_fasi xEv, double dt, double a, double b, double c)
{
  struct vettore_tri nuovo_r, delta_r1, delta_r2, r_deltar1mezzi, f2;

  delta_r1.x = xEv.v.x*dt;
  delta_r1.y = xEv.v.y*dt;
  delta_r1.z = xEv.v.z*dt;  
  r_deltar1mezzi.x = xEv.r.x + 0.5*delta_r1.x;
  r_deltar1mezzi.y = xEv.r.y + 0.5*delta_r1.y;
  r_deltar1mezzi.z = xEv.r.z + 0.5*delta_r1.z;
  f2 = derivata(r_deltar1mezzi, a, b, c);
  delta_r2.x = f2.x*dt;
  delta_r2.y = f2.y*dt;
  delta_r2.z = f2.z*dt; 
  nuovo_r.x = xEv.r.x + delta_r2.x;
  nuovo_r.y = xEv.r.y + delta_r2.y;
  nuovo_r.z = xEv.r.z + delta_r2.z;

  return nuovo_r;
}

//si implementa una funzione per equaglare due struct
struct spazio_fasi uguaglianza_spaziofasi(struct spazio_fasi old)
{
  struct spazio_fasi new;
  
  new.r.x = old.r.x;
  new.r.y = old.r.y;
  new.r.z = old.r.z;
  new.v.x = old.v.x;
  new.v.y = old.v.y;
  new.v.z = old.v.z;

  return new;
}

//si implementa una funzione che fa un'interpolazione lineare per cercare il tempo esatto in cui la velocità si annulla
double ricerca_zeri(double v_z_old, double v_z_new, double dt, int i) 
{
  double y1, y2, t1, t2, ybar = 0, tbar, a, b;

  y1 = v_z_old;
  y2 = v_z_new;
  t1 = (i-1)*dt;
  t2 = i*dt;
  a = y1-(y1-y2)/(t1-t2)*t1;
  b = (y1-y2)/(t1-t2);
  tbar = (ybar-a)/b;

  return tbar;
}
////////////////////////////
