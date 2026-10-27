//Programma finalizzato a calcolare massimi e periodo del sistema

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

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
double interpolazione(double v_z_old, double v_z_new, double dt, int i);


/*********************************************************/

int main (int argc, char *argv[])
{
  struct spazio_fasi xEv_old, xEv_new;
  int n_passi, i = 0, n_zeri_v_z = 0;
  double T_max, dt, x_0, y_0, z_0, a, b, c, T_z, A_z, A_y, A_x;
  FILE *fpr, *pFris;
  double *TEMPI_ZERI_v_z, *periodi_v_z;
  
  //controllo che sia stato inserito il numero  corretto di valori da linea di comando
  if(argc!=9)
    {
      fprintf(stderr, "#Si inseriscano i valori delle seguenti grandezze: %s, a, b, c, x0=x(0), y0=y(0), z0=z(0), delta t,  Tmax di integrazione\n", argv[0]);
      exit(1);
    }
  
  //inizializzo, ad un valore iniziale generico, la memoria degli array che conterranno i tempi degli zeri delle velocità lungo z
  TEMPI_ZERI_v_z = (double *) malloc(n_zeri_v_z * sizeof(double));
  if (TEMPI_ZERI_v_z == NULL)
    {
      printf("malloc a TEMPI_ZERI_v_z FALLITA.\n");
      exit(EXIT_FAILURE);
    }
  
  //inizializzo i valori delle condizioni iniziali e il numero di passi
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
  
  //scrivo su file i valori al tempo zero e i valori dei parametri inseriti
  fpr=fopen("../0_ToPython/integratore.dat", "w");
  if(fpr==NULL)
    {
      printf("#Errore apertura file scrittura dati\n");
      exit(EXIT_FAILURE);
    }
  
  printf("#Condizioni iniziali: a = %f, b = %f, c = %f, dt = %f, t_max = %f, x_0 = %f, y_0 = %f, z_0 = %f\n", a, b, c, dt, T_max, x_0, y_0, z_0);
  fprintf(fpr, "# t_n         x_n           y_n         z_n          v_x_n         v_y_n        v_z_n\n");
  fprintf(fpr, "%.10f %.10f %.10f %.10f %.10f %.10f %.10f\n", i*dt, xEv_old.r.x, xEv_old.r.y, xEv_old.r.z, xEv_old.v.x, xEv_old.v.y, xEv_old.v.z);
  
  //implemento l'integrazione e scrivo i valori delle coordinate e delle loro velocità ad ogni istante
  for(i=1; i<=n_passi; i++) 
    {
      xEv_new.r = Runge_Kutta_2(xEv_old, dt, a, b, c); //chiamo la funzione che implementa RK2
      xEv_new.v = derivata(xEv_new.r, a, b, c); //chiamo la funzione che implementa il sistema assegnato
      fprintf(fpr, "%.10f %.10f %.10f %.10f %.10f %.10f %.10f\n", i*dt, xEv_new.r.x, xEv_new.r.y, xEv_new.r.z, xEv_new.v.x, xEv_new.v.y, xEv_new.v.z); //stampo su file il tempo, i valori delle coordinate e delle loro velocità ad ogni passo di integrazione
      
      /********************* PARTE MASSIMI E MINIMI DI X Y E Z ***********************/ 
      if(xEv_old.v.z * xEv_new.v.z < 0) //ogni volta che questa condizione è soddisfatta si trova uno zero della velocità dunque un massimo o minimo della z
	{
	  n_zeri_v_z ++; //ogni volta che la velocità cambia di segno si incrementa un contatore del numero di zeri
	  TEMPI_ZERI_v_z = (double *) realloc (TEMPI_ZERI_v_z, n_zeri_v_z * sizeof(double)); //si alloca nuova memoria a TEMPI_ZERI pari al numero di zeri trovato fino a questo punto
	  if (TEMPI_ZERI_v_z == NULL) //si controlla il corretto funzionamento della realloc
	    {
	      printf("realloc a tempi_zeri_fin fallita.\n");
	      exit(EXIT_FAILURE);
	    }
	  TEMPI_ZERI_v_z[n_zeri_v_z - 1] = interpolazione(xEv_old.v.z, xEv_new.v.z, dt, i); //si assegna all'elemento n-1 il valore del tempo in cui si trova l' n-esimo zero della velocità
	  
	  //ogni volta che si trova uno zero della velocità si salva il valore della z in A_z così da avere una buona stima del massimo nello stato stazionario senza ricorrere ad un'interpolazione
	  if(xEv_old.r.z < xEv_new.r.z)
	    A_z = xEv_new.r.z;
	  else
	    A_z = xEv_old.r.z;
	}
      //allo stesso modo si salvano i valori dei massimi di x e y ogni volta che si trovano gli zeri delle velocià corrispondenti
      if(xEv_old.v.x * xEv_new.v.x < 0)
	{
	  if(xEv_old.r.x < xEv_new.r.x)
	    A_x = xEv_new.r.x;
	  else
	    A_x = xEv_old.r.x;
	}
      if(xEv_old.v.y * xEv_new.v.y < 0)
	{
	  if(xEv_old.r.y < xEv_new.r.y)
	    A_y = xEv_new.r.y;
	  else
	    A_y = xEv_old.r.y;
	}
      /*************** FINE PARTE MASSIMI E MINIMI X Y E Z *********************/
      xEv_old = uguaglianza_spaziofasi(xEv_new);
    }
  fclose(fpr);
  printf("Il valore massimo della coordinata x è A_x = %.15f\n", A_x);
  printf("Il valore massimo della coordinata y è A_y = %.15f\n", A_y);
  printf("Il valore massimo della coordinata z è A_z = %.15f\n", A_z);

  /***************** PARTE PERIODO ***********************/
  
  //inizializzo tramite una malloc un array dove salvare i valori di periodo
  periodi_v_z = (double *) malloc((n_zeri_v_z - 2) * sizeof(double));
  if (periodi_v_z == NULL)//controllo il corretto funzionamento della malloc
    {
      printf("malloc a periodi_v_z FALLITA.\n");
      exit(EXIT_FAILURE);
    }  
  //calcolo la distanza tra uno zero e quello due volte dopo e, per ogni coppia, salvo i valori di tempo in un'array
  for(i=0; i<n_zeri_v_z-2; i++)
    {
      periodi_v_z [i] = TEMPI_ZERI_v_z[i+2] - TEMPI_ZERI_v_z[i];
    }
  T_z = fabs(periodi_v_z [n_zeri_v_z-2] - periodi_v_z [n_zeri_v_z-3]); //assegno a T_z l'ultimo  valore di periodo
  printf("Il periodo stimato vale T = %.15f s.\n", T_z);
  /************* FINE PARTE PERIODI DI Z *********************/

  //scrivo tutti i riaultati su file
  pFris=fopen("../0_ToPython/ampiezze_periodo.dat", "w");  //apro il file in scrittura
  if(pFris==NULL) //verifico la corretta apertura del file
    {
      printf("#Errore apertura file scrittura risultati\n");
      exit(EXIT_FAILURE);
    }
  fprintf(pFris, "%f %.15f %.15f %.15f %.15f\n", dt, A_x, A_y, A_z, T_z);
  fprintf(pFris, "# dt     A_x               A_y               A_z               T\n");
  fclose(pFris);
  
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

//si implementa una funzione per eguagliare due struct
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
double interpolazione(double v_z_old, double v_z_new, double dt, int i) 
{
  double y1, y2, t1, t2, ybar = 0,  tbar, a, b;

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
