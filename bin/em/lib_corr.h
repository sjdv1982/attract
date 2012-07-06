#include "math.h"

#ifdef __cplusplus
extern "C" {
#endif

double get_product_moment(double *x, double *y, unsigned int elements);

double get_r(double product_moment_xx, double product_moment_yy, double product_moment_xy);

double get_corr(double *x, double *y, unsigned int elements);

int *get_mask(double *x, unsigned int elements, double threshold);

unsigned int apply_mask(double *source, unsigned int elements, double *target, int *mask);
#ifdef __cplusplus
} /* closing brace for extern "C" */
#endif
