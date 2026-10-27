//programma finalizzato a calcolare il vettore (x, y, z) ad un tempo T_max scelto per diversi valori di dt anche questi inseriti da terminale

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//////////////////////////////

int lettura_int (char *print);
double lettura_double (char *print);

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

/*********************************************************/

int main (int argc, char *argv[])
{
  struct spazio_fasi xEv;
  int  i = 0, j = 0, n_dt;
  double T_max, x_0, y_0, z_0, a, b, c;
  double *dt;
  int *n_passi;
  FILE *fpr;

  //verifico che il numero di parametri inseriti da linea di comando sia corretto
  if(argc!=8)
    {
      fprintf(stderr, "Si inseriscano i valori delle seguenti grandezze: %s, a, b, c, x0=x(0), y0=y(0), z0=z(0), Tmax di integrazione\n", argv[0]);
      exit(1);
    }

  //inizializzo le condizioni iniziali
  a = atof(argv[1]);
  b = atof(argv[2]);
  c = atof(argv[3]);
  x_0 = atof(argv[4]);
  y_0 = atof(argv[5]);
  z_0 = atof(argv[6]);
  T_max = atof(argv[7]);  
  
  n_dt = lettura_int("Inserire il numero di dt per i quali si vuole calcolare x(max):"); //si richiede all'utente di inserire il numero di dt
  dt = (double *) malloc(n_dt * sizeof(double)); //una volta conosciuto il numero di dt si crea un array che conterrà i loro valori
  if (dt == NULL) //si verifica il corretto funzionamento della malloc
    {
      printf("malloc a dt FALLITA.\n");
      exit(EXIT_FAILURE);
    }
  printf("# Inserire di seguito i valori dei dt:\n"); //si richiedono all'utente i valori di dt da utilizzare e li si salva nell'array dt
  for(i=0; i<n_dt; i++)
    dt[i] = lettura_double("");
  n_passi = (int *) malloc(n_dt * sizeof(int)); //si crea un array che conterrà il numero di passi corrispondente ad ogni dt
  if (dt == NULL) //si verifica il corretto funzionamento della malloc
    {
      printf("malloc a n_passi FALLITA.\n");
      exit(EXIT_FAILURE);
    }
  for(i=0; i<n_dt; i++) //per ogni dt si calcola il numero di passi corrispondente e lo si salva nell'array
    n_passi[i] = (int)(T_max/dt[i] + 0.5);

  fpr = fopen("errore_metodo.dat", "w"); //si apre il file dove verranno scritti i risultati
  if(fpr==NULL) //si controlla che il file sia stato aperto correttamente
    {
      printf("#Errore apertura file scrittura risultati\n");
      exit(EXIT_FAILURE);
    }
  //si stampano le condizioni iniziali sia su terminale che su file
  printf("# Condizioni iniziali: a = %f, b = %f, c = %f, T_max = %f, x_0 = %f, y_0 = %f, z_0 = %f\n", a, b, c, T_max, x_0, y_0, z_0);
  fprintf(fpr, "# Condizioni iniziali: a = %f, b = %f, c = %f, T_max = %f, x_0 = %f, y_0 = %f, z_0 = %f\n", a, b, c, T_max, x_0, y_0, z_0);
  printf("#T_max          x_n          y_n          z_n          dt\n");
  fprintf(fpr, "#T_max               x_n               y_n               z_n               dt\n");
  for(j=0; j<n_dt; j++) //si integra il sistema per ogni dt
    {
      xEv.r = inizializzazione(x_0, y_0, z_0);
      xEv.v = derivata(xEv.r, a, b, c);
      for(i=1; i<=n_passi[j]; i++)
	{
	  xEv.r = Runge_Kutta_2(xEv, dt[j], a, b, c);
	  xEv.v = derivata(xEv.r, a, b, c);
	}
      //si stampano, sia su terminale che su file, i valori di: T_max, (x, y, z)(T_max) e dt.
      printf("%.10f %.10f %.10f %.10f %.10f\n", T_max, xEv.r.x, xEv.r.y, xEv.r.z, dt[j]);
      fprintf(fpr, "%.15f %.15f %.15f %.15f %.15f\n", T_max, xEv.r.x, xEv.r.y, xEv.r.z, dt[j]);
    }
  return EXIT_SUCCESS;
}
  

/*********************************************************/

//funzione che mi inizializza il vettore (x, y, z) con le condizioni iniziali
struct vettore_tri inizializzazione(double x_0, double y_0, double z_0)
{
  struct vettore_tri r;
  
  r.x = x_0;
  r.y = y_0;
  r.z = z_0;
  
  return r;
}

//funzione che calcola la derivata prima del vettore (x, y, z)
struct vettore_tri derivata (struct vettore_tri r, double a, double b, double c)
{
  struct vettore_tri v;
  
  v.x = -r.y -r.z;
  v.y = r.x + a*r.y;
  v.z = b + (r.x - c)*r.z;
  
  return v;
}

//funzione che implementa il metodo di integrazione di Runge-Kutta al secondo ordine
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

////////////////////////////

//funzione per semplificare la lettura di un numero intero da terminale
int lettura_int (char *print)
{
  int input;

  printf("# %s\n", print);
  fflush(stdout);
  scanf("%d", &input);

  return input;
}

//funzione per semplificare la lettura di un numero double da terminale
double lettura_double (char *print)
{
  double input;

  printf("%s", print);
  fflush(stdout);
  scanf("%lf", &input);

  return input;
}
