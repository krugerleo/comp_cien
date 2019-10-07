#ifndef __METHOD_H__
#define __METHOD_H__

#include <stdlib.h>
#include <sys/time.h>
#define PI 3.14159265358979323846
#define ENE 4*PI*PI



typedef struct {
  double U;
  unsigned int nx; 
  unsigned int ny;
} sys_t;


double timestamp(void);
double calcNorma ();
double timestamp(void);
sys_t* inic(unsigned int nx, unsigned int ny);
void prnSistLinear (sys_t *SL);
double calFunc(int x, int y);
double coefum(double hx, double hy);
double coefdois(double hx, double hy);
double coeftres(double hx, double hy);
double coefqua(double hx, double hy);
double divi(double hx, double hy);

#endif
