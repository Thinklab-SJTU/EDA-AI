#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "geosteiner/geosteiner.h"

double eval(double* points, int* res, int degree);
int call_gst_rsmt(int nterms, double* terms, double* length, int* nsps, double* sps, int* nedges, int* edges) ;