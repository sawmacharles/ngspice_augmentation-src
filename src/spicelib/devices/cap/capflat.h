/* card #532 (round 15): flattened capacitor instance arrays for the tran load */
#ifndef CAPFLAT_H_
#define CAPFLAT_H_
typedef struct { int n; int *pos, *neg, *q; double *m, *cap; } CAPflat;
int CAPbuildBase(GENmodel *, const double *, unsigned int, double *, CAPflat *);
int CAPloadFlat(CKTcircuit *, const CAPflat *);
void CAPflatFree(CAPflat *);
#endif
