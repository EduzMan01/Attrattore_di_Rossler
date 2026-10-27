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

/*********************************************************/

int main (int argc, char *argv[])
{
  struct spazio_fasi xEv;
  int n_passi, i = 0;
  double T_max, dt, x_0, y_0, z_0, a, b, c;

  //controllo che il numero di argomenti inseriti da linea di comando sia corretto, in caso contrario viene comunicata la sintassi da seguire e si termina il programma
  if(argc!=9)
    {
      fprintf(stderr, "Si inseriscano i valori delle seguenti grandezze: %s, a, b, c, x0=x(0), y0=y(0), z0=z(0), delta t,  T_max di integrazione\n", argv[0]);
      exit(EXIT_FAILURE);
    }
  //si inizializzano i valori passati da riga di comando ed alcune grandezze iniziali quali il numero di passi
  a = atof(argv[1]);
  b = atof(argv[2]);
  c = atof(argv[3]);
  x_0 = atof(argv[4]);
  y_0 = atof(argv[5]);
  z_0 = atof(argv[6]);
  dt = atof(argv[7]);
  T_max = atof(argv[8]);

  n_passi = (int)(T_max/dt + 0.5);
  xEv.r = inizializzazione(x_0, y_0, z_0);
  xEv.v = derivata(xEv.r, a, b, c);

  //si stampano le condizioni iniziali
  printf("#Condizioni iniziali: a = %f, b = %f, c = %f, dt = %f, T_max = %f, x_0 = %f, y_0 = %f, z_0 = %f\n", a, b, c, dt, T_max, x_0, y_0, z_0);
  printf("# t_n   x_n   y_n   z_n\n");
  printf("%.10f %.10f %.10f %.10f\n", i*dt,  xEv.r.x, xEv.r.y, xEv.r.z);
  for(i=1; i<=n_passi; i++)
    {
      xEv.r = Runge_Kutta_2(xEv, dt, a, b, c); //si richiama la funzione RK2
      xEv.v = derivata(xEv.r, a, b, c); //si calcola la velocità ad ogni istante
      printf("%.10f %.10f %.10f %.10f\n", i*dt, xEv.r.x, xEv.r.y, xEv.r.z); //si stampano i valori di tempo e posizione per ogni istante
    }
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
