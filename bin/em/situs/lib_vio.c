/*********************************************************************
*                            L I B _ V I O                           *
**********************************************************************
* Program is part of the Situs package (c) Willy Wriggers, 1998-2009 *
* URL: situs.biomachina.org                                          *
**********************************************************************
*                                                                    * 
* Auxiliary program to read and write EM maps in Situs format        *
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

/* 7 March 2010: Adapted by Sjoerd de Vries (SJdV), adaptations marked*/

#include "situs.h"
#include "lib_vio.h" 
#include "lib_err.h"

void read_vol(char *vol_file, double *width, double *origx, double *origy, double *origz, unsigned *extx, unsigned *exty, unsigned *extz, double **phi) {
	unsigned long nvox, count;
	double dorigx, dorigy, dorigz, dwidth;
	double phitest, dtemp;
	char *program = "lib_vio";
	FILE *fin;
	
	fin = fopen(vol_file, "r");
	if( fin == NULL ) {
		error_open_filename(13010, program, vol_file);
	}
	
	/* read header and print information */
	fscanf(fin, "%le %le %le %le %d %d %d", &dwidth, &dorigx, &dorigy, &dorigz, extx, exty, extz);
	*width = dwidth; *origx = dorigx; *origy = dorigy; *origz = dorigz;
/* Adapted by SJdV: old:
	printf ("lib_vio> File %s - Header information: \n", vol_file);
	printf ("lib_vio> Columns, rows, and sections: x=1-%d, y=1-%d, z=1-%d\n",*extx,*exty,*extz);
	printf ("lib_vio> 3D coordinates of first voxel: (%f,%f,%f)\n",*origx,*origy,*origz);
	printf ("lib_vio> Voxel size in Angstrom: %f \n", *width);
new: */
	fprintf (stderr,"lib_vio> File %s - Header information: \n", vol_file);
	fprintf (stderr,"lib_vio> Columns, rows, and sections: x=1-%d, y=1-%d, z=1-%d\n",*extx,*exty,*extz);
	fprintf (stderr,"lib_vio> 3D coordinates of first voxel: (%f,%f,%f)\n",*origx,*origy,*origz);
	fprintf (stderr,"lib_vio> Voxel size in Angstrom: %f \n", *width);

	nvox = *extx * *exty * *extz;
	
	/* allocate memory and read data */
/* Adapted by SJdV: old:	
	printf ("lib_vio> Reading density data... \n");
new: */	
        fprintf (stderr, "lib_vio> Reading density data... \n");
	*phi = (double *) malloc(nvox * sizeof(double));
	if (*phi == NULL) {
		error_memory_allocation(13020, program);
	}
	
	for (count=0;count<nvox;count++) {
		if (fscanf(fin,"%le", &dtemp) != 1) {
			error_unreadable_file_long(13030, program, vol_file);
		} else *(*phi+count) = dtemp;
	}
	if (fscanf(fin,"%le", &phitest) != EOF) {
		error_unreadable_file_long(13040, program, vol_file);
	}
	fclose(fin);
/* Adapted by SJdV: old:	
	printf ("lib_vio> Volumetric data read from file %s\n", vol_file); 
new: */	
        fprintf (stderr, "lib_vio> Volumetric data read from file %s\n", vol_file); 
	return;
}

void write_vol(char *vol_file, double width, double origx, double origy, double origz, unsigned extx, unsigned exty, unsigned extz, double *phi) {
	unsigned long nvox, count;
	char *program = "lib_vio";
	FILE *fout;
	 
	nvox = extx * exty * extz;
	fout = fopen(vol_file, "w");
	if( fout == NULL ) {
		error_open_filename(13210, program, vol_file);
	}
	
	printf ("lib_vio> Writing density data... \n");
	fprintf(fout, "%f %f %f %f %d %d %d\n", width, origx, origy, origz, extx, exty, extz);
	fprintf(fout, "\n");
	
	for(count=0;count<nvox;count++) {
		if ((count+1)%10 == 0) fprintf (fout," %10.6f \n",*(phi+count));
		else fprintf (fout," %10.6f ",*(phi+count));
	}
	fclose(fout);
	printf ("lib_vio> Volumetric data written to file %s \n", vol_file);
	
	/* header information */
	printf ("lib_vio> File %s - Header information: \n", vol_file);
	printf ("lib_vio> Columns, rows, and sections: x=1-%d, y=1-%d, z=1-%d\n",extx,exty,extz);
	printf ("lib_vio> 3D coordinates of first voxel: (%f,%f,%f)\n",origx,origy,origz);
	printf ("lib_vio> Voxel size in Angstrom: %f \n", width);

	return;
}
