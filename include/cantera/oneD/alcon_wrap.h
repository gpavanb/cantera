#ifndef alcon_wrap_h
#define alcon_wrap_h

extern "C" {

void c_f(double*,double*,double*);
void c_alcon1(void (Cantera::Sim1D::*c_f)(double*,double*,double*),int*,double*,double*,
              double*,double*,double*,double*,int*,double*,int*,int*,int*);

}

#endif
