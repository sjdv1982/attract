#include "lib_corr.h"

extern "C" double get_product_moment(double *x, double *y, unsigned int elements) {
  double sumx = 0, sumy = 0, sumxy = 0;
  for (unsigned int n = 0; n < elements; n++) {
    double &xx = x[n];
    double &yy = y[n];
    sumx += xx; sumy += yy;
    sumxy += xx * yy;
  }
  return sumxy - sumx * sumy/ elements;
}

extern "C" double get_r(double product_moment_xx, double product_moment_yy, double product_moment_xy) {
  double r = product_moment_xy/sqrt(product_moment_xx * product_moment_yy);
  return r;
}

extern "C" double get_corr(double *x, double *y, unsigned int elements) {
  double Sxx = get_product_moment(x,x,elements);
  double Syy = get_product_moment(y,y,elements);
  double Sxy = get_product_moment(x,y,elements);
  double r = get_r(Sxx,Syy,Sxy);
  return r;
}

extern "C" int *get_mask(double *x, unsigned int elements, double threshold) {
  int *mask = new int[elements];
  for (unsigned int n = 0; n < elements; n++) {
    mask[n] = (x[n] > threshold);
  }
  return mask;
}

extern "C" unsigned int apply_mask(double *source, unsigned int elements, double *target, int *mask) {
  unsigned int target_elements = 0;
  for (unsigned int n = 0; n < elements; n++) {
    if (mask[n]) {
      target[target_elements] = source[n];
      target_elements++;      
    }
  }
  return target_elements;
}
