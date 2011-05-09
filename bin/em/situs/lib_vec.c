/*********************************************************************
*                           L I B _ V E C                            *
**********************************************************************
* Library is part of the Situs package URL: situs.biomachina.org     *
* (c) Pablo Chacon and Willy Wriggers, 2001-2003                     *
**********************************************************************
*                                                                    *
* Creating, resetting, copying arrays and other data structures.     * 
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

#include "situs.h"
#include "lib_vec.h"
#include "lib_err.h"

/*====================================================================*/
void zero_vect(double *vect, unsigned long len) {
	unsigned long i;
	for(i=0;i<len;++i) vect[i] = 0.0;
}

/*====================================================================*/
void do_vect(double **vect, unsigned long len) {
	char *program = "lib_vec";
	*vect = (double *) malloc(len*sizeof(double));
	if (*vect == NULL) {
		error_memory_allocation(18010, program);
	}
	zero_vect(*vect,len);
}

/*====================================================================*/
void zero_mat(double **mat,unsigned long len_i,unsigned long len_j) {
	unsigned long i;
	for(i=0;i<len_i;++i) zero_vect(mat[i],len_j);
}

/*====================================================================*/
void do_mat(double ***pmat,unsigned long len_i,unsigned long len_j) {
	unsigned long i;
	char *program = "lib_vec";
	
	*pmat = (double **) malloc(len_i*sizeof(double *));
	if (*pmat == NULL) {
		error_memory_allocation(18020, program);
	}
	for(i=0;i<len_i;i++) do_vect(&((*pmat)[i]),len_j);
}

/*====================================================================*/
void cp_vect(double **vect1,double **vect2,unsigned long len) {
	memcpy(*vect1,*vect2, len*sizeof(double));
}

/*====================================================================*/
/* destroys memory allocated to vect2 after copying */
void cp_vect_destroy(double **vect1,double **vect2,unsigned long len) {
	free(*vect1); 
	do_vect(vect1,len); 
	cp_vect(vect1,vect2,len); 
	free(*vect2);
}

/*====================================================================*/
void add_scaled_vect(double *to_vect, double *from_vect, double scalar, unsigned long len) { 
	unsigned long i;
	for(i=0;i<len;++i) to_vect[i] += scalar * from_vect[i];
}



