#ifndef gc_wrap_h
#define gc_wrap_h

///////////////////////
// GroupContribution //
//////////////////////
extern "C" {
void init_gc(const char* fuel);
int numSpecies();
// Dependent properties
void rho_l(double*,double*);
void rho_l_PR(double*,double*,double*);
void Hv(double*,double*);
void c_l(double*,double*);
void mu_g(double*,double*,double*);
void mu_l(double*,double*,double*);
void D(double*,double*,double*);
void PSat(double*,double*);
void surf(double*,double*);

// Independent properties
void MW(double*);

void clean_gc();
}

void initialize_gc(std::string fuel) {
  init_gc(fuel.data());
}

#endif /* gc_wrap_h */
