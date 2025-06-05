#include "eval.h"

double eval(double* points, int* res, int degree) {
    double length = 0;
    double xs[degree], ys[degree];
    double xl[degree], xr[degree];
    double yl[degree], yh[degree];
    int i, j;
    for (i = 0, j = 0; i < degree;) {
        xs[i] = points[j++];
        ys[i] = points[j++];
        i++;
    }
    size_t cpy_len = sizeof(double) * degree;
    memcpy(xl, xs, cpy_len);
    memcpy(xr, xs, cpy_len);
    memcpy(yl, ys, cpy_len);
    memcpy(yh, ys, cpy_len);
    int x_idx, y_idx;
    double e_x, e_y, val;
    int limit = degree + degree - 2;
    for (i = 0; i < limit; ) {
        x_idx = res[i++];
        y_idx = res[i++];
        e_x = xs[x_idx];
        e_y = ys[y_idx];
        val = yl[x_idx]; yl[x_idx] = (val < e_y ? val : e_y);
        val = xl[y_idx]; xl[y_idx] = (val < e_x ? val : e_x);
        val = yh[x_idx]; yh[x_idx] = (val > e_y ? val : e_y);
        val = xr[y_idx]; xr[y_idx] = (val > e_x ? val : e_x);
    }
    for (i = 0; i < degree; i++) {
        length += yh[i] - yl[i] + xr[i] - xl[i];
    }
    return length;
}

int call_gst_rsmt(int nterms, double* terms, double* length, int* nsps, double* sps, int* nedges, int* edges) {
    return gst_rsmt(nterms, terms, length, nsps, sps, nedges, edges, NULL, NULL);
}
