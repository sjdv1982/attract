/*********************************************************************
*                         M A P 2 M A P                              *
**********************************************************************
* Program is part of the Situs package (c) Willy Wriggers, 1998-2009 *
* URL: situs.biomachina.org                                          *
**********************************************************************
*                                                                    *
* Map file format conversion.                                        * 
*                                                                    *
**********************************************************************
* See legal statement for terms of distribution                      *
*********************************************************************/

#include "situs.h"
#include "lib_vio.h"
#include "lib_vwk.h"
#include "lib_vec.h"
#include "lib_std.h"
#include "lib_err.h"

static void read_ascii (char *, unsigned long, double **);
static void read_xplor (char *, int *, int *, double *, double *, double *, double *, double *, double *, int *, int *, int *, unsigned *, unsigned *, unsigned *, double **);
static void xplor_skip_to_number (FILE **, char **);  
static void read_mrc2000 (char *, int *, int *, int*, unsigned *, unsigned *, unsigned *, int *, int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double **);
static void read_spider (char *, unsigned *, unsigned *, unsigned *, double **);
static void dump_binary_and_exit (char *, char *, int);
static int read_float (float *, FILE *, int);
static int read_float_empty (FILE *);
static int read_char_float (float *, FILE *);
static int read_char (char *, FILE *);
static int read_int (int *, FILE *, int);
static unsigned long count_floats (FILE **);
static int test_registration(float, float, float, float);
static int test_mrc2000 (char *, int);
static int test_spider (char *, int);
static unsigned long permuted_index (int, unsigned long, unsigned, unsigned, unsigned);
static void permute_map (int, unsigned, unsigned, unsigned, unsigned *, unsigned *, unsigned*, int, int, int, int *, int *, int *, double *, double **);
static void permute_dimensions (int, unsigned, unsigned, unsigned, unsigned *, unsigned *, unsigned *, int, int, int, int *, int *, int *);
static void permute_print_info (int, unsigned, unsigned, unsigned, unsigned, unsigned, unsigned, int, int, int, int, int, int); 
static int set_origin_get_mode (int, int, int, double, double, double, double, double, double, double *, double *, double *);
static void assert_cubic_map (int, int, double, double, double, double, double, double, unsigned, unsigned, unsigned, int, int, int, double, double, double, double *, double *, double *, double *, unsigned *, unsigned *, unsigned *, double **);
static void interpolate_skewed_map_to_cubic (double **, unsigned *, unsigned *, unsigned *, double *, double *, double *, double *, double *, unsigned, unsigned, unsigned, int, int, int, double, double, double, double, double, double, double, double, double);
static void write_xplor (char *, double, double, double, double, unsigned, unsigned, unsigned, double *);
static void write_mrc2000 (int, char *, double, double, double, double, unsigned, unsigned, unsigned, double *);
static void write_spider (char *, double, double, double, double, unsigned, unsigned, unsigned, double *);


int main (int argc, char *argv[]) {
  double porigx, porigy, porigz;  
  double *phi;
  double *pphi;
  unsigned pextx, pexty, pextz;
  unsigned long nvox;
  int menumode, ordermode = 7, swapmode, cubic = 1, orom = 1;
  double pwidth, widthx, widthy, widthz;
  double alpha, beta, gamma;
  int nxstart=0, nystart=0, nzstart=0;
  int ncstart=0, nrstart=0, nsstart=0;
  unsigned nc, nr, ns;
  unsigned nx, ny, nz;
  int currext;
  double xorigin, yorigin, zorigin;
  char ac='X', ar='Y', as='Z';
  
  if (argc != 3) {
    printf ("map2map> Usage: map2map inputfile outputfile\n");
    exit (1);
  }
  
  printf ("map2map> \n");
  printf ("map2map> Situs file format conversion for cubic lattice maps typically used in EM.\n");
  printf ("map2map> Choose one of the following 10 options.\n");
  printf ("map2map> \n");
  printf ("map2map> OTHER -> Situs: \n");
  printf ("map2map> \n");
  printf ("map2map>      1: ASCII (editable text) file, sequential list of map densities** \n");
  printf ("map2map>      2: MRC2000 or CCP4 binary (auto)* \n");
  printf ("map2map>      3: MRC2000 or CCP4 binary (manual)** \n");
  printf ("map2map>      4: SPIDER binary** \n");
  printf ("map2map>      5: X-PLOR map (ASCII editable text)* \n");
  printf ("map2map>      6: Generic 4-byte field binary (unknown map or header parameters) \n");
  printf ("map2map> \n");
  printf ("map2map> Situs -> OTHER: \n");
  printf ("map2map> \n");
  printf ("map2map>      7: Situs -> MRC2000 / CCP4 binary (auto)*\n");
  printf ("map2map>      8: Situs -> MRC2000 / CCP4 binary (manual)**\n"); 
  printf ("map2map>      9: Situs -> SPIDER binary* \n");
  printf ("map2map>     10: Situs -> X-PLOR (ASCII editable text)* \n"); 
  printf ("map2map> \n");
  printf ("map2map>      *: automatic fill of header fields \n");
  printf ("map2map>     **: manual assignment of header fields \n");
  printf ("map2map> \n");
  printf ("map2map> ");
  
  menumode = readln_int();
  
  /************** process input map ******************/ 
  
  switch (menumode) {
  case 1: /* free format ASCII */    
    /* get order and map size parameters */
    printf ("map2map> \n");
    printf ("map2map> Data order and axis permutation in file %s.\n", argv[1]);
    printf ("map2map> Assign columns (C, fastest), rows (R), and sections (S, slowest) to X, Y, Z:\n");
    printf ("map2map> \n");
    printf ("map2map>         C  R  S = \n");
    printf ("map2map> ------------------\n");
    printf ("map2map>      1: X  Y  Z (no permutation)\n");
    printf ("map2map>      2: X  Z  Y\n");
    printf ("map2map>      3: Y  X  Z\n");
    printf ("map2map>      4: Y  Z  X\n");
    printf ("map2map>      5: Z  X  Y\n");
    printf ("map2map>      6: Z  Y  X\n");
    printf ("map2map> Enter selection number [1-6]: ");
    ordermode = readln_int();
    
    /* assign axis character */
    switch (ordermode) {
    case 1: 
      ac = 'X'; ar = 'Y'; as = 'Z';
      break;
    case 2:
      ac = 'X'; ar = 'Z'; as = 'Y';
      break;
    case 3:
      ac = 'Y'; ar = 'X'; as = 'Z';
      break;
    case 4:
      ac = 'Y'; ar = 'Z'; as = 'X';
      break; 
    case 5:
      ac = 'Z'; ar = 'X'; as = 'Y';
      break;
    case 6:
      ac = 'Z'; ar = 'Y'; as = 'X';
      break;
    default:
      error_option(70209, "map2map");
    }

    printf ("map2map> Enter number of columns (%c fields): ",ac);
    currext = readln_int();
    if (currext < 1) {
      error_number_columns(70010, "map2map");
    } else nc = currext;
    printf ("map2map> Enter number of rows (%c fields): ",ar);
    currext = readln_int();
    if (currext < 1) {
      error_number_rows(70020, "map2map");
    } else nr = currext;
    printf ("map2map> Enter number of sections (%c fields): ",as);
    currext = readln_int();
    if (currext < 1) {
      error_number_sections(70030, "map2map"); 
    } else ns = currext;

    /* read map and permute */
    nvox = nc * nr * ns;  
    read_ascii(argv[1],nvox, &phi);
    permute_map (ordermode, nc, nr, ns, &pextx, &pexty, &pextz, ncstart, nrstart, nsstart, &nxstart, &nystart, &nzstart, phi, &pphi);
  
    /* get lattice spacing and origin */
    printf ("map2map> Enter desired cubic grid (lattice) spacing in Angstrom: ");
    pwidth = readln_double();
    if (pwidth <= 0) {
      error_number_spacing(70040, "map2map");
    } 
    printf("map2map> \n");
    printf("map2map> Enter X-origin (coord of first voxel) in Angstrom (0 is recommended, a different value shifts the resulting map accordingly): ");
    porigx = readln_double();
    printf("map2map> \n");
    printf("map2map> Enter Y-origin (coord of first voxel) in Angstrom (0 is recommended, a different value shifts the resulting map accordingly): ");
    porigx = readln_double();
    printf("map2map> \n");
    printf("map2map> Enter Z-origin (coord of first voxel) in Angstrom (0 is recommended, a different value shifts the resulting map accordingly): ");
    porigx = readln_double();
    printf("map2map> \n");
    break;
    
  case 2: /* MRC2000 or CCP4 format with automatically filled Situs header information */

    read_mrc2000 (argv[1], &orom, &cubic, &ordermode, &nc, &nr, &ns, &ncstart, &nrstart, &nsstart, &widthx, &widthy, &widthz, &xorigin, &yorigin, &zorigin, &alpha, &beta, &gamma, &phi);
    permute_map (ordermode, nc, nr, ns, &nx, &ny, &nz, ncstart, nrstart, nsstart, &nxstart, &nystart, &nzstart, phi, &pphi);
    permute_print_info(ordermode,nc,nr,ns,nx,ny,nz,ncstart,nrstart,nsstart,nxstart,nystart,nzstart);
    assert_cubic_map (orom, cubic, alpha, beta, gamma, widthx, widthy, widthz, nx, ny, nz, nxstart, nystart, nzstart, xorigin, yorigin, zorigin, &pwidth, &porigx, &porigy, &porigz, &pextx, &pexty, &pextz, &pphi);
    break;
    
  case 3: /* MRC2000 or CCP4 format with manual override of Situs header information */
    
    read_mrc2000 (argv[1], &orom, &cubic, &ordermode, &nc, &nr, &ns, &ncstart, &nrstart, &nsstart, &widthx, &widthy, &widthz, &xorigin, &yorigin, &zorigin, &alpha, &beta, &gamma, &phi);
    nvox = nc * nr * ns;
    permute_map (ordermode, nc, nr, ns, &nx, &ny, &nz, ncstart, nrstart, nsstart, &nxstart, &nystart, &nzstart, phi, &pphi);
    permute_print_info(ordermode,nc,nr,ns,nx,ny,nz,ncstart,nrstart,nsstart,nxstart,nystart,nzstart);
    assert_cubic_map (orom, cubic, alpha, beta, gamma, widthx, widthy, widthz, nx, ny, nz, nxstart, nystart, nzstart, xorigin, yorigin, zorigin, &pwidth, &porigx, &porigy, &porigz, &pextx, &pexty, &pextz, &pphi);
 
    /* offer to override Situs header parameters */
    printf("map2map> \n");
    printf("map2map> The currently assigned voxel spacing is %f Angstrom\n",pwidth);
    printf("map2map> To keep the field, enter this value again (different value rescales the map): ");
    pwidth = readln_double();
    if (pwidth <= 0) {
      error_number_spacing(70140, "map2map");
    } 
    printf("map2map> \n");
    printf("map2map> The currently assigned X-origin (coord of first voxel) is %f Angstrom\n",porigx);
    printf("map2map> To keep the field, enter this value again (different value shifts the map): ");
    porigx = readln_double();
    printf("map2map> \n");
    printf("map2map> The currently assigned Y-origin (coord of first voxel) is %f Angstrom\n",porigy);
    printf("map2map> To keep the field, enter this value again (different value shifts the map): ");
    porigy = readln_double();
    printf("map2map> \n");
    printf("map2map> The currently assigned Z-origin (coord of first voxel) is %f Angstrom\n",porigz);
    printf("map2map> To keep the field, enter this value again (different value shifts the map): ");
    porigz = readln_double();
    printf("map2map> \n");
    printf("map2map> The currently assigned NX (# X fields) value is %d \n",pextx);
    printf("map2map> Enter the same or a new value: ");
    pextx = readln_int();
    printf("map2map> The currently assigned NY (# Y fields) value is %d \n",pexty);
    printf("map2map> Enter the same or a new value: ");
    pexty = readln_int();
    printf("map2map> The currently assigned NZ (# Z fields) value is %d \n",pextz);
    printf("map2map> Enter the same or a new value: ");
    pextz = readln_int();
    if (pextx*pexty*pextz != nvox) {
      fprintf (stderr,"map2map> Sorry, NX * NY * NZ must match the total number of voxels, %ld,  in the map.\n",nvox);
      fprintf (stderr,"map2map> Please try again, bye bye.\n");
      exit (1);
    }
    break;
    
  case 4: /* SPIDER format with manual assignment of Situs header information */
    read_spider (argv[1], &pextx, &pexty, &pextz, &pphi);
    printf("map2map> \n");
    printf("map2map> SPIDER maps don't typically contain scale and origin information relative to a PDB coordinate system.\n");
    printf("map2map> Enter desired grid (lattice) spacing in Angstrom: ");
    pwidth = readln_double();
    if (pwidth <= 0) {
      error_number_spacing(70141, "map2map");
    } 
    printf("map2map> \n");
    printf("map2map> Enter X-origin (coord of first voxel) in Angstrom (0 is recommended, a different value shifts the resulting map accordingly): ");
    porigx = readln_double();
    printf("map2map> \n");
    printf("map2map> Enter Y-origin (coord of first voxel) in Angstrom (0 is recommended, a different value shifts the resulting map accordingly): ");
    porigx = readln_double();
    printf("map2map> \n");
    printf("map2map> Enter Z-origin (coord of first voxel) in Angstrom (0 is recommended, a different value shifts the resulting map accordingly): ");
    porigx = readln_double();
    printf("map2map> \n");
    break;
    
  case 5: /* X-PLOR ASCII format */
    read_xplor (argv[1], &orom, &cubic, &widthx, &widthy, &widthz, &alpha, &beta, &gamma, &nxstart, &nystart, &nzstart, &nx, &ny, &nz, &pphi);
    assert_cubic_map (orom, cubic, alpha, beta, gamma, widthx, widthy, widthz, nx, ny, nz, nxstart, nystart, nzstart, 0.0, 0.0, 0.0, &pwidth, &porigx, &porigy, &porigz, &pextx, &pexty, &pextz, &pphi);
    break;
    
  case 6: /* generic binary */
    printf ("map2map> \n");
    printf ("map2map> Do you want to swap the byte order (endianism)?\n");
    printf ("map2map> \n");
    printf ("map2map>      1: No (order is correct for current machine architecture)\n");
    printf ("map2map>      2: Yes (foreign machine binary)\n");
    printf ("map2map> ");
    swapmode = readln_int()-1;
    dump_binary_and_exit(argv[1],argv[2],swapmode);
    break; 
    
  case 7: case 8: case 9: case 10: /* Situs format */
    ordermode = 1;
    read_vol(argv[1], &pwidth, &porigx, &porigy, &porigz, &pextx, &pexty, &pextz, &pphi);
    break;
  }
  
  /************** write output map ******************/ 
  
  switch (menumode) {
  case 1: case 2: case 3: case 4: case 5: /* Situs format output */
    write_vol(argv[2], pwidth, porigx, porigy, porigz, pextx, pexty, pextz, pphi);
    printf("map2map> \n");
    printf("map2map> All done.\n");
    break; 
    
  case 6: /* generic binary already dumped, this option should not be reached */
    fprintf(stderr, "map2map> Error: Unable to process option [e.c. 70205]\n");
    exit(70205);
    break;
    
  case 7: /* hybrid MRC2000 / CCP4 format, auto */
    write_mrc2000 (1, argv[2], pwidth, porigx, porigy, porigz, pextx, pexty, pextz, pphi);
    break;

  case 8: /* hybrid MRC2000 / CCP4 format, manual */    
    write_mrc2000 (0, argv[2], pwidth, porigx, porigy, porigz, pextx, pexty, pextz, pphi);
    break;

  case 9: /* SPIDER format */
    write_spider (argv[2], pwidth, porigx, porigy, porigz, pextx, pexty, pextz, pphi);
    break;
    
  case 10: /* X-PLOR format */
    write_xplor (argv[2], pwidth, porigx, porigy, porigz, pextx, pexty, pextz, pphi);
    break;
  }
  return 0;
}
	


/* reads ASCII stream */ 
static void read_ascii (char *vol_file, unsigned long nvox, double **fphi) {
  unsigned long count;
  FILE *fin;
  float currfloat;
  
  fin = fopen(vol_file, "r");
  if( fin == NULL ) {
    error_open_filename(70310, "map2map", vol_file);
  }
  printf ("map2map> Reading ASCII data... \n");
  do_vect(fphi,nvox);
  for (count=0;count<nvox;count++) { 
    if (fscanf(fin,"%e", &currfloat) != 1) {
      error_unreadable_file_short(70330 , "map2map", vol_file);
    } else {
      *(*fphi+count) = currfloat;
    }
  }
  if (fscanf(fin,"%e", &currfloat) != EOF) {
    error_unreadable_file_long(70340 , "map2map", vol_file);
  }
  printf ("map2map> Volumetric data read from file %s\n", vol_file);
  fclose (fin);
}


/* reads X-PLOR ASCII file */ 
static void read_xplor (char *vol_file, int *orom, int *cubic,
			double *widthx, double *widthy, double *widthz,
			double *alpha, double *beta, double *gamma, 
			int *nxstart, int *nystart, int *nzstart,
			unsigned *extx, unsigned *exty, unsigned *extz, double **fphi) {
  
  unsigned long nvox;
  FILE *fin;
  int idummy, done;
  char *nextline;
  float a_tmp, b_tmp, g_tmp, currfloat;
  long mx, my, mz, mxend, myend, mzend;
  long mxstart, mystart, mzstart; 
  int testa, testb, testg;
  float xlen, ylen, zlen;
  int indx, indy, indz;
  
  nextline = (char *) malloc(FLENGTH * sizeof(char));
  if (nextline == NULL) {
    error_memory_allocation(70410, "map2map");
  }
  
  fin = fopen(vol_file, "r");
  if( fin == NULL ) {
    error_open_filename(70420, "map2map", vol_file);
  }
  
  /* ignore header length line */
  xplor_skip_to_number (&fin, &nextline);
  
  /* read index line */
  xplor_skip_to_number (&fin, &nextline);
  if (sscanf (nextline, "%8ld%8ld%8ld%8ld%8ld%8ld%8ld%8ld%8ld",&mx,&mxstart,&mxend,&my,&mystart,&myend,&mz,&mzstart,&mzend) != 9) {
    error_xplor_file_indexing(70430, "map2map");
  } 
  *extx = mxend-mxstart+1;*exty = myend-mystart+1;*extz = mzend-mzstart+1; 
  nvox = *extx * *exty * *extz;
  
  printf("map2map> X-PLOR map indexing (counting from 0): \n");
  printf("map2map>       NA = %8ld  (# of X intervals in unit cell) \n",mx);
  printf("map2map>     AMIN = %8ld  (start index X) \n",mxstart);
  printf("map2map>     AMAX = %8ld  (end index X) \n",mxend);
  printf("map2map>       NB = %8ld  (# of Y intervals in unit cell) \n",my);
  printf("map2map>     BMIN = %8ld  (start index Y) \n",mystart);
  printf("map2map>     BMAX = %8ld  (end index Y) \n",myend);
  printf("map2map>       NC = %8ld  (# of Z intervals in unit cell) \n",mz);
  printf("map2map>     CMIN = %8ld  (start index Z) \n",mzstart);
  printf("map2map>     CMAX = %8ld  (end index Z) \n",mzend);
  
  /* read unit cell info and determine grid width and origin */
  xplor_skip_to_number (&fin, &nextline);
  if (sscanf (nextline, "%12f%12f%12f%12f%12f%12f",&xlen,&ylen,&zlen,&a_tmp,&b_tmp,&g_tmp) != 6) {
    error_xplor_file_unit_cell(70440, "map2map");
  } else {
    *alpha = a_tmp; *beta = b_tmp; *gamma = g_tmp;
  }
  
  printf("map2map> X-PLOR unit cell info: \n");
  printf("map2map>        A = %8.3f  (unit cell dimension) \n",xlen);
  printf("map2map>        B = %8.3f  (unit cell dimension) \n",ylen);
  printf("map2map>        C = %8.3f  (unit cell dimension) \n",zlen);
  printf("map2map>    ALPHA = %8.3f  (unit cell angle) \n",*alpha);
  printf("map2map>     BETA = %8.3f  (unit cell angle) \n",*beta);
  printf("map2map>    GAMMA = %8.3f  (unit cell angle) \n",*gamma);
  printf("map2map> \n");
  
  /* assign voxel spacing parameters */
  *widthx = xlen / (double) mx;
  *widthy = ylen / (double) my;
  *widthz = zlen / (double) mz;
  *nxstart = mxstart; *nystart = mystart; *nzstart = mzstart; 

  /* test for orthogonal and cubic lattice */
  testa = floor(100* *alpha+0.5); testb = floor(100* *beta+0.5); testg = floor(100* *gamma+0.5);
  if (testa != 9000 || testb != 9000 || testg != 9000) *orom = 0;  
  if (*orom == 0 || floor((*widthx-*widthy)*1000+0.5)!=0 || floor((*widthy-*widthz)*1000+0.5)!=0 || floor((*widthx-*widthz)*1000+0.5)!=0) *cubic = 0;
    
  /* read ZYX info */
  for (done=0;done==0;) {/* read next line and check if it contains ZYX */ 
    if(fgets(nextline,FLENGTH,fin)==NULL) {
      error_EOF_ZYX_mode(70450, "map2map");
    }
    if (*nextline=='Z' && *(nextline+1)=='Y' && *(nextline+2)=='X') { done = 1; break; }
    if (*nextline=='z' && *(nextline+1)=='y' && *(nextline+2)=='x') { done = 1; break; }
  }
  
  /* read sections */
  do_vect(fphi,nvox);
  for(indz=0;indz<*extz;indz++) {
    /* skip section header and read section number */
    xplor_skip_to_number (&fin, &nextline);
    if (sscanf (nextline, "%8d",&idummy) != 1) {
      error_xplor_file_map_section(70470,"map2map");
    }
    if (idummy!=indz) {
      error_xplor_file_map_section_number(70480,"map2map");
    } 
    
    /* read section data */
    for(indy=0;indy<*exty;indy++) for(indx=0;indx<*extx;indx++) {
      if (fscanf(fin,"%12f",&currfloat) != 1) {
	error_unreadable_file_short(70490, "map2map", vol_file);
      } else {
	*(*fphi+gidz_general(indz,indy,indx,*exty,*extx))=currfloat;
      }
    }
  }

  /* read end of data marker */
  xplor_skip_to_number (&fin, &nextline);
  if (sscanf (nextline, "%8d",&idummy) != 1 || idummy != -9999) {
    error_xplor_maker("map2map");
  }
  
  fclose(fin);
  free (nextline);
}


/* reads lines from current position in opened X-PLOR file */
/* stops if current line contains numbers before any '!' */
static void xplor_skip_to_number (FILE **fin, char **nextline) {
  int i, done, foundletter;
  
  for (done=0;done==0;) {/* read next line */ 
    if(fgets(*nextline,FLENGTH,*fin)==NULL) {
      error_EOF(70510, "map2map");
    }
    for (i=0;i<FLENGTH;++i) if ((*(*nextline+i)=='!') || (*(*nextline+i)=='\n')) { *(*nextline+i)='\0'; break; } 
    foundletter = 0;
    for (i=0;*(*nextline+i)!='\0';++i) {/* check if it contains ABC FGHIJKLMNOPQRSTUVWXYZ (or lower case) */
      if (*(*nextline+i) > 64 && *(*nextline+i) < 68) { foundletter = 1; break; }
      if (*(*nextline+i) > 69 && *(*nextline+i) < 91) { foundletter = 1; break; }
      if (*(*nextline+i) > 96 && *(*nextline+i) < 100) { foundletter = 1; break; }
      if (*(*nextline+i) > 101 && *(*nextline+i) < 123) { foundletter = 1; break; }
    }
    if (foundletter == 0) for (i=0;*(*nextline+i)!='\0';++i) 
      if (*(*nextline+i) >= '0' && *(*nextline+i) <= '9') { done = 1; break; }
  } 
}


/* reads MRC2000 or CCP4 binary file and swaps bytes automatically */ 
static void read_mrc2000 (char *vol_file, int *orom, int *cubic, int *ordermode,
			  unsigned *nc, unsigned *nr, unsigned *ns, 
			  int *ncstart, int *nrstart, int *nsstart,
			  double *widthx, double *widthy, double *widthz, 
			  double *xorigin, double *yorigin, double *zorigin,
			  double *alpha, double *beta, double *gamma,
			  double **fphi) {
  
  unsigned long count, nvox;
  FILE *fin;
  int nc_tmp, nr_tmp, ns_tmp, mx, my, mz;
  int mode; 
  float a_tmp, b_tmp, g_tmp;
  float x_tmp, y_tmp, z_tmp;
  int testa, testb, testg;
  float xlen, ylen, zlen;
  int i, swap, header_ok = 1;
  int mapc, mapr, maps;
  float dmin, dmax, dmean, dummy, currfloat;
  int n_range_viol0, n_range_viol1;
  char mapchar1, mapchar2, mapchar3, mapchar4;
  
  n_range_viol0 = test_mrc2000(vol_file,0);
  n_range_viol1 = test_mrc2000(vol_file,1);
  
  if (n_range_viol0 < n_range_viol1) { /* guess endianism */
    swap = 0;
    if (n_range_viol0 > 0) {
      printf("map2map> Warning: %i header field range violations detected \n", n_range_viol0); 
    }
  } else {
    swap = 1;
    if (n_range_viol1 > 0) {
      printf("map2map> Warning: %i header field range violations detected \n", n_range_viol1); 
    }
  }
  
  /* read header */
  fin = fopen(vol_file, "rb");
  if( fin == NULL ) {
    error_open_filename(70620, "map2map", vol_file);
  }
  printf("map2map> Reading header information from MRC2000 or CCP4 file %s \n", vol_file); 
  header_ok *= read_int(&nc_tmp,fin,swap);
  header_ok *= read_int(&nr_tmp,fin,swap);
  header_ok *= read_int(&ns_tmp,fin,swap);
  *nc = nc_tmp; *nr = nr_tmp; *ns = ns_tmp;
  header_ok *= read_int(&mode,fin,swap);
  header_ok *= read_int(ncstart,fin,swap);
  header_ok *= read_int(nrstart,fin,swap);
  header_ok *= read_int(nsstart,fin,swap);
  header_ok *= read_int(&mx,fin,swap);
  header_ok *= read_int(&my,fin,swap);
  header_ok *= read_int(&mz,fin,swap);
  header_ok *= read_float(&xlen,fin,swap);
  header_ok *= read_float(&ylen,fin,swap);
  header_ok *= read_float(&zlen,fin,swap);
  header_ok *= read_float(&a_tmp,fin,swap);
  header_ok *= read_float(&b_tmp,fin,swap);
  header_ok *= read_float(&g_tmp,fin,swap);
  *alpha = a_tmp; *beta = b_tmp; *gamma = g_tmp;
  header_ok *= read_int(&mapc,fin,swap);
  header_ok *= read_int(&mapr,fin,swap);
  header_ok *= read_int(&maps,fin,swap);
  header_ok *= read_float(&dmin,fin,swap);
  header_ok *= read_float(&dmax,fin,swap);
  header_ok *= read_float(&dmean,fin,swap);
  for (i=23; i<50; ++i) header_ok *= read_float(&dummy,fin,swap);
  header_ok *= read_float(&x_tmp,fin,swap);
  header_ok *= read_float(&y_tmp,fin,swap);
  header_ok *= read_float(&z_tmp,fin,swap);
  *xorigin = x_tmp; *yorigin = y_tmp; *zorigin = z_tmp;
  header_ok *= read_char(&mapchar1,fin);
  header_ok *= read_char(&mapchar2,fin);
  header_ok *= read_char(&mapchar3,fin);
  header_ok *= read_char(&mapchar4,fin);
  if (header_ok == 0) {
    error_file_header(70650, "map2map", vol_file);
  } 
  
  /* print some info */
  printf("map2map>       NC = %8d  (# columns)\n",*nc); 
  printf("map2map>       NR = %8d  (# rows)\n",*nr); 
  printf("map2map>       NS = %8d  (# sections)\n",*ns); 
  printf("map2map>     MODE = %8d  (data type: 0: 1-byte char, 2: 4-byte float)\n",mode); 
  printf("map2map>  NCSTART = %8d  (index of first column, counting from 0)\n",*ncstart); 
  printf("map2map>  NRSTART = %8d  (index of first row, counting from 0)\n",*nrstart);
  printf("map2map>  NSSTART = %8d  (index of first section, counting from 0)\n",*nsstart);
  printf("map2map>       MX = %8d  (# of X intervals in unit cell)\n",mx); 
  printf("map2map>       MY = %8d  (# of Y intervals in unit cell)\n",my);
  printf("map2map>       MZ = %8d  (# of Z intervals in unit cell)\n",mz);
  printf("map2map> X length = %8.3f  (unit cell dimension)\n",xlen);
  printf("map2map> Y length = %8.3f  (unit cell dimension)\n",ylen);
  printf("map2map> Z length = %8.3f  (unit cell dimension)\n",zlen);
  printf("map2map>    Alpha = %8.3f  (unit cell angle)\n",*alpha);
  printf("map2map>     Beta = %8.3f  (unit cell angle)\n",*beta);
  printf("map2map>    Gamma = %8.3f  (unit cell angle)\n",*gamma);
  printf("map2map>     MAPC = %8d  (columns axis: 1=X,2=Y,3=Z)\n",mapc); 
  printf("map2map>     MAPR = %8d  (rows axis: 1=X,2=Y,3=Z)\n",mapr);
  printf("map2map>     MAPS = %8d  (sections axis: 1=X,2=Y,3=Z)\n",maps);
  printf("map2map>     DMIN = %8.3f  (minimum density value)\n",dmin);
  printf("map2map>     DMAX = %8.3f  (maximum density value)\n",dmax);
  printf("map2map>    DMEAN = %8.3f  (mean density value)\n",dmean);
  printf("map2map>  XORIGIN = %8.3f  (X origin - MRC2000 only)\n",*xorigin);
  printf("map2map>  YORIGIN = %8.3f  (Y origin - MRC2000 only)\n",*yorigin);
  printf("map2map>  ZORIGIN = %8.3f  (Z origin - MRC2000 only)\n",*zorigin);
  printf("map2map>      MAP = %c%c%c%c      (map string)\n",mapchar1,mapchar2,mapchar3,mapchar4);
  
  /* extract data based on file type mode, currently supports 1-byte char and 4-byte float */
  /* note: fphi will contain data in the order of the input file, no axis permutation yet */
  nvox = *nc * *nr * *ns;
  switch (mode) {
  case 0: /* char - converted to float */
    rewind (fin);
    do_vect(fphi,nvox);
    for (count=0; count<256; ++count) if (read_float_empty(fin)==0) {
      error_file_convert(70642, "map2map", vol_file);
    } 
    for (count=0; count<nvox; ++count) {
      if (read_char_float(&currfloat,fin)==0) {
	error_file_convert(70651, "map2map", vol_file);
      } else {
	*(*fphi+count)=currfloat;
      }
    } 
    fclose (fin);
    break;
  case 2: /* float */
    rewind (fin);
    do_vect(fphi,nvox);
    for (count=0; count<256; ++count) if (read_float_empty(fin)==0) {
      error_file_convert(70642, "map2map", vol_file);
    } 
    for (count=0; count<nvox; ++count) {
      if (read_float(&currfloat,fin,swap)==0) {
	error_file_convert(70652, "map2map", vol_file);
      } else {
	*(*fphi+count)=currfloat;
      }
    } 
    fclose (fin);
    break;
  default:
    error_file_float_mode(70653,"map2map", vol_file);
  }

  /* assign voxel spacing parameters */
  *widthx = xlen / (double) mx;
  *widthy = ylen / (double) my;
  *widthz = zlen / (double) mz;
  
  /* test for orthogonal and cubic lattice */
  testa = floor(100* *alpha+0.5); testb = floor(100* *beta+0.5); testg = floor(100* *gamma+0.5);
  if (testa != 9000 || testb != 9000 || testg != 9000) *orom = 0;
  if (*orom == 0 || floor((*widthx-*widthy)*1000+0.5)!=0 || floor((*widthy-*widthz)*1000+0.5)!=0 || floor((*widthx-*widthz)*1000+0.5)!=0) *cubic = 0;
  
  /* set axis permutation mode */
  *ordermode = 7;
  if (mapc==1 && mapr==2 && maps==3) {
    *ordermode = 1;
  }
  if (mapc==1 && mapr==3 && maps==2) {
    *ordermode = 2;
  }
  if (mapc==2 && mapr==1 && maps==3) {
    *ordermode = 3;
  }
  if (mapc==2 && mapr==3 && maps==1) {
    *ordermode = 4;    
  }
  if (mapc==3 && mapr==1 && maps==2) {
    *ordermode = 5;
  }
  if (mapc==3 && mapr==2 && maps==1) {
    *ordermode = 6;
  }
  if (*ordermode == 7) {
    error_axis_assignment(70680, "map2map");
  }
}


/* reads SPIDER binary file and swaps bytes automatically */ 
static void read_spider (char *vol_file, unsigned *extx, unsigned *exty, unsigned *extz, double **fphi) {
  unsigned long nvox;
  FILE *fin;
  int i, swap, header_ok = 1, headlen;
  unsigned long count; 
  float dummy, currfloat;
  float nslice, nrow, iform, imami, fmax, fmin, av, sig, nsam, headrec;
  float iangle, phi, theta, gamma, xoff, yoff, zoff, scale, labbyt, lenbyt;
  float istack, inuse, maxim, kangle, phi1, theta1, psi1, phi2, theta2, psi2;
  int n_range_viol0, n_range_viol1;
  
  n_range_viol0 = test_spider(vol_file,0);
  n_range_viol1 = test_spider(vol_file,1);
  
  if (n_range_viol0 < n_range_viol1) { /* guess endianism */
    swap = 0;
    if (n_range_viol0 > 0) {
      printf("map2map> Warning: %i header field range violations detected \n", n_range_viol0); 
    }
  } else {
    swap = 1;
    if (n_range_viol1 > 0) {
      printf("map2map> Warning: %i header field range violations detected \n", n_range_viol1); 
    }
  }
  
  /* read header */
  fin = fopen(vol_file, "rb");
  if( fin == NULL ) {
    error_open_filename(70820, "map2map", vol_file);
  } 	
  printf("map2map> Reading header information from SPIDER file %s \n", vol_file); 
  header_ok *= read_float(&nslice,fin,swap);
  header_ok *= read_float(&nrow,fin,swap);
  header_ok *= read_float(&dummy,fin,swap);
  header_ok *= read_float(&dummy,fin,swap);
  header_ok *= read_float(&iform,fin,swap);
  header_ok *= read_float(&imami,fin,swap);
  header_ok *= read_float(&fmax,fin,swap);
  header_ok *= read_float(&fmin,fin,swap);
  header_ok *= read_float(&av,fin,swap);
  header_ok *= read_float(&sig,fin,swap);
  header_ok *= read_float(&dummy,fin,swap);
  header_ok *= read_float(&nsam,fin,swap);
  header_ok *= read_float(&headrec,fin,swap);
  header_ok *= read_float(&iangle,fin,swap);
  header_ok *= read_float(&phi,fin,swap);
  header_ok *= read_float(&theta,fin,swap);
  header_ok *= read_float(&gamma,fin,swap);
  header_ok *= read_float(&xoff,fin,swap);
  header_ok *= read_float(&yoff,fin,swap);
  header_ok *= read_float(&zoff,fin,swap);
  header_ok *= read_float(&scale,fin,swap);
  header_ok *= read_float(&labbyt,fin,swap);
  header_ok *= read_float(&lenbyt,fin,swap);
  header_ok *= read_float(&istack,fin,swap);
  header_ok *= read_float(&inuse,fin,swap);
  header_ok *= read_float(&maxim,fin,swap);
  for (i=0; i<4; ++i) header_ok *= read_float(&dummy,fin,swap);
  header_ok *= read_float(&kangle,fin,swap);
  header_ok *= read_float(&phi1,fin,swap);
  header_ok *= read_float(&theta1,fin,swap);
  header_ok *= read_float(&psi1,fin,swap);
  header_ok *= read_float(&phi2,fin,swap);
  header_ok *= read_float(&theta2,fin,swap);
  header_ok *= read_float(&psi2,fin,swap);
  if (header_ok == 0) {
    error_file_header(70850, "map2map", vol_file);
  } 
  
  /* print some info */
  printf("map2map>   NSLICE = %8.f  (# sections)\n",nslice); 
  printf("map2map>     NROW = %8.f  (# rows)\n",nrow); 
  printf("map2map>    IFORM = %8.f  (ignored: file type specifier)\n",iform); 
  printf("map2map>    IMAMI = %8.f  (flag: 1 = the maximum and minimum are computed)\n",imami); 
  printf("map2map>     FMAX = %8.3f  (maximum density value)\n",fmax); 
  printf("map2map>     FMIN = %8.3f  (minimum density value)\n",fmin);
  printf("map2map>       AV = %8.3f  (average density value)\n",av);
  printf("map2map>      SIG = %8.3f  (standard deviation of density distribution)\n",sig); 
  printf("map2map>     NSAM = %8.f  (# columns)\n",nsam);
  printf("map2map>  HEADREC = %8.f  (number of records in file header)\n",headrec);
  printf("map2map>   IANGLE = %8.f  (flag: 1 = tilt angles filled)\n",iangle);
  printf("map2map>      PHI = %8.3f  (ignored: tilt angle)\n",phi);
  printf("map2map>    THETA = %8.3f  (ignored: tilt angle)\n",theta);
  printf("map2map>    GAMMA = %8.3f  (ignored: tilt angle)\n",gamma);
  printf("map2map>     XOFF = %8.3f  (ignored: X offset)\n",xoff);
  printf("map2map>     YOFF = %8.3f  (ignored: Y offset)\n",yoff);
  printf("map2map>     ZOFF = %8.3f  (ignored: Z offset)\n",zoff); 
  printf("map2map>    SCALE = %8.3f  (ignored: scale factor)\n",scale);
  printf("map2map>   LABBYT = %8.f  (total number of bytes in header)\n",labbyt);
  printf("map2map>   LENBYT = %8.f  (record length in bytes)\n",lenbyt);
  printf("map2map>   ISTACK = %8.f  (flag; file contains a stack of images)\n",istack);
  printf("map2map>    INUSE = %8.f  (flag; this image in stack is used)\n",inuse);
  printf("map2map>    MAXIM = %8.f  (ignored: maximum image used in stack)\n",maxim);
  printf("map2map>   KANGLE = %8.f  (flag; additional angles set)\n",kangle);
  printf("map2map>     PHI1 = %8.3f  (ignored: additional rotation)\n",phi1);
  printf("map2map>   THETA1 = %8.3f  (ignored: additional rotation)\n",theta1);
  printf("map2map>     PSI1 = %8.3f  (ignored: additional rotation)\n",psi1);
  printf("map2map>     PHI2 = %8.3f  (ignored: additional rotation)\n",phi2);
  printf("map2map>   THETA2 = %8.3f  (ignored: additional rotation)\n",theta2);
  printf("map2map>     PSI2 = %8.3f  (ignored: additional rotation)\n",psi2);
  printf("map2map> \n");   
  
  *extx = nsam; *exty = nrow; *extz = nslice; nvox = *extx * *exty * *extz;
  headlen = *extx * ceil(256/(*extx*1.0));
  do_vect(fphi,nvox);
  rewind(fin);
  for (count=0; count<headlen; ++count) if (read_float_empty(fin)==0) {
    error_file_convert(70841, "map2map", vol_file);
  } 
  for (count=0; count<nvox; ++count) if (read_float(&currfloat,fin,swap)==0) {
    error_file_convert(70842, "map2map", vol_file);
  } else {
    *(*fphi+count) = currfloat;
  }
  fclose(fin);

  /* notes and error checks */
  if (((int) floor(100*(sizeof(float)*headlen-labbyt)+0.5)) != 0) {
    error_spider_header(70860, "map2map");
  }
}


/* ASCII dumps binary file, flips bytes if swap==1 */ 
static void dump_binary_and_exit (char *in_file, char *out_file, int swap) {
  FILE *fin, *fout;
  float *phi;
  unsigned long nfloat;
  unsigned long count;
  
  fin = fopen(in_file, "rb");
  if( fin == NULL ) {
    error_open_filename(70910, "map2map", in_file);
  }
  
  nfloat = count_floats(&fin);
  
  phi = (float *) malloc(nfloat * 4); 
  if (phi == NULL) {
    error_memory_allocation(70920, "map2map");
  }
  
  for(count=0;count<nfloat;count++) read_float(phi+count,fin,swap);
  
  printf("map2map> Binary data read as 'float' type from file %s\n", in_file);
  if (swap) printf("map2map> The byte order (endianism) has been swapped.\n");
  fclose (fin);
  
  fout = fopen(out_file, "w");
  if( fout == NULL ) {
    error_open_filename(70940, "map2map", out_file);
  }
  
  printf ("map2map> Writing ASCII data... \n"); 
  for(count=0;count<nfloat;count++) {
    if ((count+1)%10 == 0) fprintf (fout," %10.6f \n",*(phi+count));
    else fprintf (fout," %10.6f ",*(phi+count)); 
  }
  fclose(fout);
  
  printf("map2map> ASCII dump of binary file %s written to file %s. \n", in_file, out_file);
  printf("map2map> Open / check ASCII file with text editor and extract map densities. \n");
  printf("map2map> Then convert to Situs format with option 1: ASCII.\n");
  printf("map2map> \n");
  printf("map2map> All done.\n");
  exit(1);
}

/* reads float and swaps bytes if swap==1 */ 
static int read_float (float *currfloat, FILE *fin, int swap) {
  unsigned char *cptr, tmp;
  
  if (fread(currfloat,4,1,fin)!=1) return 0;
  if (swap == 1) {
    cptr = (unsigned char *)currfloat;
    tmp = cptr[0];
    cptr[0]=cptr[3];
    cptr[3]=tmp;
    tmp = cptr[1];
    cptr[1]=cptr[2];
    cptr[2]=tmp;
  }
  return 1;
}

/* reads header float and does nothing */ 
static int read_float_empty (FILE *fin) {
  float currfloat;
  
  if (fread(&currfloat,4,1,fin)!=1) return 0;
  return 1;
}

/* reads char and assigns to float */ 
static int read_char_float (float *currfloat, FILE *fin) {
  char currchar;
  
  if (fread(&currchar,1,1,fin)!=1) return 0;
  *currfloat=(float)currchar;
  return 1;
}

/* reads char and assigns to char */ 
static int read_char (char *currchar, FILE *fin) {
	
  if (fread(currchar,1,1,fin)!=1) return 0;
  return 1;
}

/* reads int and swaps bytes if swap==1 */ 
static int read_int (int *currlong, FILE *fin, int swap) {
  unsigned char *cptr, tmp;
  
  if (fread(currlong,4,1,fin)!=1) return 0;
  if (swap == 1) {
    cptr = (unsigned char *)currlong;
    tmp = cptr[0];
    cptr[0]=cptr[3];
    cptr[3]=tmp;
    tmp = cptr[1];
    cptr[1]=cptr[2];
    cptr[2]=tmp;
  }
  return 1;
}

/* counts number of 4-byte floats of open file */ 
static unsigned long count_floats (FILE **fin) {
  unsigned long finl = 0;
  
  rewind(*fin);
  for( ; ; ) {
    if( fgetc(*fin) == EOF ) break;
    ++finl;
  }
  rewind(*fin);
  return finl / 4;
}
	
/* checks for registration of grid with origin of coordinate system */ 
static int test_registration(float origx, float origy, float origz, float width) {
  float xreg, xreg1, yreg, yreg1, zreg, zreg1;

  xreg = fabs(fmod(origx+0.00001*width,width));
  xreg1 = fabs(fmod(origx-0.00001*width,width));
  yreg = fabs(fmod(origy+0.00001*width,width));
  yreg1 = fabs(fmod(origy-0.00001*width,width));
  zreg = fabs(fmod(origz+0.00001*width,width));
  zreg1 = fabs(fmod(origz-0.00001*width,width));
  if (xreg1<xreg) xreg = xreg1;
  if (yreg1<yreg) yreg = yreg1;
  if (zreg1<zreg) zreg = zreg1;
  if (xreg+yreg+zreg>0.0001) return 0;
  else return 1;
}
	


/* tests MRC2000 / CCP4 header and returns (non-exhaustive) number of range violations */ 
static int test_mrc2000 (char *vol_file, int swap) {
  FILE *fin;
  int nc, nr, ns, mx, my, mz;
  int mode, ncstart, nrstart, nsstart;
  float xlen, ylen, zlen;
  int i, header_ok = 1, n_range_viols = 0;
  int mapc, mapr, maps;
  float alpha, beta, gamma;
  float dmin, dmax, dmean, dummy, xorigin, yorigin, zorigin;
  
  fin = fopen(vol_file, "rb");
  if( fin == NULL ) {
    error_open_filename(71010, "map2map", vol_file);
  }
  
  /* read header info */
  header_ok *= read_int(&nc,fin,swap);
  header_ok *= read_int(&nr,fin,swap);
  header_ok *= read_int(&ns,fin,swap);
  header_ok *= read_int(&mode,fin,swap);
  header_ok *= read_int(&ncstart,fin,swap);
  header_ok *= read_int(&nrstart,fin,swap);
  header_ok *= read_int(&nsstart,fin,swap);
  header_ok *= read_int(&mx,fin,swap);
  header_ok *= read_int(&my,fin,swap);
  header_ok *= read_int(&mz,fin,swap);
  header_ok *= read_float(&xlen,fin,swap);
  header_ok *= read_float(&ylen,fin,swap);
  header_ok *= read_float(&zlen,fin,swap);
  header_ok *= read_float(&alpha,fin,swap);
  header_ok *= read_float(&beta,fin,swap);
  header_ok *= read_float(&gamma,fin,swap);
  header_ok *= read_int(&mapc,fin,swap);
  header_ok *= read_int(&mapr,fin,swap);
  header_ok *= read_int(&maps,fin,swap);
  header_ok *= read_float(&dmin,fin,swap);
  header_ok *= read_float(&dmax,fin,swap);
  header_ok *= read_float(&dmean,fin,swap);
  for (i=23; i<50; ++i) header_ok *= read_float(&dummy,fin,swap);
  header_ok *= read_float(&xorigin,fin,swap);
  header_ok *= read_float(&yorigin,fin,swap);
  header_ok *= read_float(&zorigin,fin,swap);
  fclose (fin);
  if (header_ok == 0) {
    error_file_header(71020, "map2map", vol_file);
  } 
  
  n_range_viols += (nc>5000); n_range_viols += (nc<0);
  n_range_viols += (nr>5000); n_range_viols += (nr<0);
  n_range_viols += (ns>5000); n_range_viols += (ns<0);
  n_range_viols += (ncstart>5000); n_range_viols += (ncstart<-5000);
  n_range_viols += (nrstart>5000); n_range_viols += (nrstart<-5000);
  n_range_viols += (nsstart>5000); n_range_viols += (nsstart<-5000);
  n_range_viols += (mx>5000); n_range_viols += (mx<0);
  n_range_viols += (my>5000); n_range_viols += (my<0);
  n_range_viols += (mz>5000); n_range_viols += (mz<0);
  n_range_viols += (alpha>360.0f); n_range_viols += (alpha<-360.0f);
  n_range_viols += (beta>360.0f); n_range_viols += (beta<-360.0f);
  n_range_viols += (gamma>360.0f); n_range_viols += (gamma<-360.0f);
  
  return n_range_viols;
}


/* tests SPIDER header and returns (non-exhaustive) number of range violations */
static int test_spider (char *vol_file, int swap) {
  FILE *fin;
  int i, header_ok = 1, n_range_viols = 0, headlen;
  float dummy;
  float nslice, nrow, iform, imami, fmax, fmin, av, sig, nsam, headrec;
  float iangle, phi, theta, gamma, xoff, yoff, zoff, scale, labbyt, lenbyt;
  float istack, inuse, maxim, kangle, phi1, theta1, psi1, phi2, theta2, psi2;
  
  fin = fopen(vol_file, "rb");
  if( fin == NULL ) {
    error_open_filename(71210, "map2map", vol_file);
  }
  
  /* read header info */
  header_ok *= read_float(&nslice,fin,swap);
  header_ok *= read_float(&nrow,fin,swap);
  header_ok *= read_float(&dummy,fin,swap);
  header_ok *= read_float(&dummy,fin,swap);
  header_ok *= read_float(&iform,fin,swap);
  header_ok *= read_float(&imami,fin,swap);
  header_ok *= read_float(&fmax,fin,swap);
  header_ok *= read_float(&fmin,fin,swap);
  header_ok *= read_float(&av,fin,swap);
  header_ok *= read_float(&sig,fin,swap);
  header_ok *= read_float(&dummy,fin,swap);
  header_ok *= read_float(&nsam,fin,swap);
  header_ok *= read_float(&headrec,fin,swap);
  header_ok *= read_float(&iangle,fin,swap);
  header_ok *= read_float(&phi,fin,swap);
  header_ok *= read_float(&theta,fin,swap);
  header_ok *= read_float(&gamma,fin,swap);
  header_ok *= read_float(&xoff,fin,swap);
  header_ok *= read_float(&yoff,fin,swap);
  header_ok *= read_float(&zoff,fin,swap);
  header_ok *= read_float(&scale,fin,swap);
  header_ok *= read_float(&labbyt,fin,swap);
  header_ok *= read_float(&lenbyt,fin,swap);
  header_ok *= read_float(&istack,fin,swap);
  header_ok *= read_float(&inuse,fin,swap);
  header_ok *= read_float(&maxim,fin,swap);
  for (i=0; i<4; ++i) header_ok *= read_float(&dummy,fin,swap);
  header_ok *= read_float(&kangle,fin,swap);
  header_ok *= read_float(&phi1,fin,swap);
  header_ok *= read_float(&theta1,fin,swap);
  header_ok *= read_float(&psi1,fin,swap);
  header_ok *= read_float(&phi2,fin,swap);
  header_ok *= read_float(&theta2,fin,swap);
  header_ok *= read_float(&psi2,fin,swap);
  fclose (fin);
  if (header_ok == 0) {
    error_file_header(71220, "map2map", vol_file);
  } 
  headlen = nsam * ceil(256/(nsam*1.0)); 
  n_range_viols += (((int) floor(100*(4*headlen-labbyt)+0.5)) != 0);
  n_range_viols += (headrec>1e10); n_range_viols += (headrec<=1e-10);
  n_range_viols += (labbyt>1e10); n_range_viols += (labbyt<=1e-10);
  n_range_viols += (lenbyt>1e10); n_range_viols += (lenbyt<=1e-10);
  n_range_viols += (nslice>1e10); n_range_viols += (nslice<1e-10);
  n_range_viols += (nrow>1e10); n_range_viols += (nrow<1e-10);
  n_range_viols += (nsam>1e10); n_range_viols += (nsam<1e-10);
  
  return n_range_viols;
}


/* returns the permuted index */
/* nc, nr, ns, are the unperturbed dimensions (# colums, # rows, # sections) */
/* computes first from 'count' the indices ic, ir, is, assumed to correspond to original, unpermuted input map */
/* returns new index m = ix + iy * nx + iz * nx * ny, where the x,y,z mapping is done (implicitly) as in permute_dimensions */
/* new(m) = old(count) will create the permuted (x,y,z) map */ 
static unsigned long permuted_index (int ordermode, unsigned long count, unsigned nc, unsigned nr, unsigned ns) {
  unsigned ic, ir, is; 
  unsigned long ncr, q;
  
  ncr = nc*nr;
  is = count / ncr;
  q = count - is*ncr;
  ir = q / nc;
  ic = q - ir*nc;
  
  switch (ordermode) {
  case 1:
    return ic+ir*nc+is*nc*nr;
  case 2:
    return ic+is*nc+ir*nc*ns;
  case 3:
    return ir+ic*nr+is*nr*nc;
  case 4:
    return is+ic*ns+ir*ns*nc;
  case 5:
    return ir+is*nr+ic*nr*ns;
  case 6:
    return is+ir*ns+ic*ns*nr;
  default:
    error_option(70210, "map2map");
    return 0;
  }
}

/* create permuted map and reorder dimensions */
static void permute_map (int ordermode, unsigned nc, unsigned nr, unsigned ns, unsigned *nx, unsigned *ny, unsigned *nz, 
			 int ncstart, int nrstart, int nsstart, int *nxstart, int *nystart, int *nzstart, double *phi, double **pphi) {
  
  unsigned long nvox, count;

  /* create permuted map pphi */
  nvox = nc * nr * ns;
  do_vect(pphi,nvox);
  for(count=0;count<nvox;count++) *(*pphi+permuted_index(ordermode,count,nc,nr,ns)) = *(phi+count);
  free(phi);
  
  /* permute the map dimensions */
  permute_dimensions(ordermode,nc,nr,ns,nx,ny,nz,ncstart,nrstart,nsstart,nxstart,nystart,nzstart);
}


/* create permuted dimensions */
static void permute_dimensions (int ordermode, unsigned nc, unsigned nr, unsigned ns, unsigned *nx, unsigned *ny, unsigned *nz, int ncstart, int nrstart, int nsstart, int *nxstart, int *nystart, int *nzstart) {
  
  switch (ordermode) {
  case 1: 
    *nx = nc; *ny = nr; *nz = ns;
    *nxstart = ncstart; *nystart = nrstart; *nzstart = nsstart;
    break;
  case 2:
    *nx = nc; *ny = ns; *nz = nr;
    *nxstart = ncstart; *nystart = nsstart; *nzstart = nrstart;
    break;
  case 3:
    *nx = nr; *ny = nc; *nz = ns;
    *nxstart = nrstart; *nystart = ncstart; *nzstart = nsstart;
    break;
  case 4:
    *nx = ns; *ny = nc; *nz = nr;
    *nxstart = nsstart; *nystart = ncstart; *nzstart = nrstart;
    break;
  case 5:
    *nx = nr; *ny = ns; *nz = nc;
    *nxstart = nrstart; *nystart = nsstart; *nzstart = ncstart;
    break;
  case 6:
    *nx = ns; *ny = nr; *nz = nc;
    *nxstart = nsstart; *nystart = nrstart; *nzstart = ncstart;
    break;
  default:
    error_option(70212, "map2map");
  }  
}


/* print information on permuted dimensions */
static void permute_print_info (int ordermode, unsigned nc, unsigned nr, unsigned ns, unsigned nx, unsigned ny, unsigned nz, int ncstart, int nrstart, int nsstart, int nxstart, int nystart, int nzstart) {
  
  if (ordermode > 1) { 
    printf("map2map> \n");
    printf("map2map> Map parameters BEFORE selected axis permutation: \n");
    printf("map2map>       NC = %8d  (# columns)\n",nc); 
    printf("map2map>       NR = %8d  (# rows)\n",nr); 
    printf("map2map>       NS = %8d  (# sections)\n",ns); 
    printf("map2map>  NCSTART = %8d  (index of first column, counting from 0)\n",ncstart); 
    printf("map2map>  NRSTART = %8d  (index of first row, counting from 0)\n",nrstart);
    printf("map2map>  NSSTART = %8d  (index of first section, counting from 0)\n",nsstart);
    printf("map2map> \n");
    printf("map2map> Map parameters AFTER selected axis permutation: \n");
    printf("map2map>       NX = %8d  (# X fields)\n",nx); 
    printf("map2map>       NY = %8d  (# Y fields)\n",ny); 
    printf("map2map>       NZ = %8d  (# Z fields)\n",nz); 
    printf("map2map>  NXSTART = %8d  (X index of first voxel, counting from 0)\n",nxstart); 
    printf("map2map>  NYSTART = %8d  (Y index of first voxel, counting from 0)\n",nystart);
    printf("map2map>  NZSTART = %8d  (Z index of first voxel, counting from 0)\n",nzstart);   
    printf("map2map> \n");
  } else {
    printf("map2map> \n");
    printf("map2map> No axis permutation present, assigned map dimension parameters: \n");
    printf("map2map>            NC,NX = %8d  (# columns = # X fields)\n",nc); 
    printf("map2map>            NR,NY = %8d  (# rows = # Y fields)\n",nr); 
    printf("map2map>            NS,NZ = %8d  (# sections # Z fields)\n",ns); 
    printf("map2map> NCSTART, NXSTART = %8d  (X index of first voxel, counting from 0)\n",ncstart); 
    printf("map2map> NRSTART, NYSTART = %8d  (Y index of first voxel, counting from 0)\n",nrstart);
    printf("map2map> NSSTART, NSSTART = %8d  (Z index of first voxel, counting from 0)\n",nsstart);
    printf("map2map> \n");
  }
}


/* sets origin based on crystallographic or MRC2000 conventions, returns 0 for crystallographic, 1 for MRC2000 style */
static int set_origin_get_mode (int nxstart, int nystart, int nzstart, double widthx, double widthy, double widthz, double xorigin, double yorigin, double zorigin, double *origx, double *origy, double *origz) {
 
  if (fabs(xorigin)<0.0001 && fabs(yorigin)<0.0001 && fabs(zorigin)<0.0001) { /* seem to have crystallographic origin */
    printf("map2map> Using crystallographic style origin defined by unit cell start indices.\n");
    *origx = nxstart * widthx; *origy = nystart * widthy; *origz = nzstart * widthz;
    return 0;
  } else { /* seem to have MRC2000 */
    printf("map2map> Using MRC2000 style origin defined by [X,Y,Z]ORIGIN fields.\n");
    *origx = xorigin; *origy = yorigin; *origz = zorigin;
    return 1;
  }
}

/* if necessary, map non-cubic maps to cubic lattice by interpolation */
static void assert_cubic_map (int orom, int cubic, double alpha, double beta, double gamma, double widthx, double widthy, double widthz, unsigned nx, unsigned ny, unsigned nz, 
			      int nxstart, int nystart, int nzstart, double xorigin, double yorigin, double zorigin, double *pwidth, double *porigx, double *porigy, 
			      double *porigz, unsigned *pextx, unsigned *pexty, unsigned *pextz, double **pphi) {

  unsigned long nvox;
  double *pp2;
  double origx, origy, origz;
  
  /* distinguish between cubic, orthonormal, and skewed */
  if (cubic) {
    set_origin_get_mode(nxstart,nystart,nzstart,widthx,widthy,widthz,xorigin,yorigin,zorigin,porigx,porigy,porigz);
    printf("map2map> Cubic lattice present. \n");
    *pwidth = widthx; *pextx = nx; *pexty = ny; *pextz = nz; 
  } else {
    if (orom) {
      set_origin_get_mode(nxstart,nystart,nzstart,widthx,widthy,widthz,xorigin,yorigin,zorigin,&origx,&origy,&origz);
      printf("map2map> Orthogonal lattice with unequal spacings in x, y, z detected.\n");
      printf("map2map> Interpolating map to cubic lattice using smallest detected spacing rounded to nearest 0.1 Angstrom.\n");
      /* for interpolation, find minimum voxel spacing, rounded to nearest 0.1 Angstrom for good measure */
      *pwidth=widthx; if(widthy<*pwidth) *pwidth=widthy; if(widthz<*pwidth) *pwidth=widthz;
      *pwidth = floor(*pwidth*10.0+0.5)/10.0;
      /* the new map will have the maximum box dimensions that do not exceed the old map */
      interpolate_map (&pp2, pextx, pexty, pextz, porigx, porigy, porigz, 
		       *pwidth, *pwidth, *pwidth, *pphi, nx, ny, nz, origx, 
		       origy, origz, widthx, widthy, widthz);
      nvox = *pextx * *pexty * *pextz;
      cp_vect_destroy(pphi,&pp2,nvox);
    } else {
      printf("map2map> Skewed, non-orthogonal lattice detected.\n");
      printf("map2map> Interpolating skewed map to cubic lattice using smallest detected spacing rounded to nearest 0.1 Angstrom.\n");
      /* here, the new map will have the minimum box dimensions that fully enclose the old map */
      interpolate_skewed_map_to_cubic (&pp2, pextx, pexty, pextz, porigx, porigy, porigz, 
				       pwidth, *pphi, nx, ny, nz, nxstart, nystart, nzstart, xorigin, yorigin, zorigin, 
				       widthx, widthy, widthz, alpha, beta, gamma);
      nvox = *pextx * *pexty * *pextz;
      cp_vect_destroy(pphi,&pp2,nvox);
    }
  }   

  /* check for registration of lattice with origin of coordinate system */
  if (test_registration(*porigx,*porigy,*porigz,*pwidth)==0) 
    fprintf (stderr,"map2map> Input grid not in register with origin of coordinate system.\n");
  printf("map2map> \n");
}


/* interpolate skewed map to cubic lattice within rectangular bounding box */
/* pphiout is allocated and new output map parameters are returned */
static void interpolate_skewed_map_to_cubic (double **pphiout, unsigned *pextx, unsigned *pexty, unsigned *pextz, 
					     double *porigx, double *porigy, double *porigz, double *pwidth,
					     double *phiin, 
					     unsigned nx, unsigned ny, unsigned nz, 
					     int nxstart, int nystart, int nzstart, 
					     double xorigin, double yorigin, double zorigin,
					     double widthx, double widthy, double widthz,
					     double alpha, double beta, double gamma) {
  
  unsigned long pnvox;
  double ax, ay, az, bx, by, bz, cx, cy, cz, aix, aiy, aiz, bix, biy, biz, cix, ciy, ciz;
  double t1x, t1y, t2x, t2y, t3x, t3y, t4x, t4y, cdet, scz;
  double ux[8], uy[8], uz[8], uxmin, uxmax, uymin, uymax, uzmin, uzmax;
  double xpos, ypos, zpos, gx, gy, gz, a, b, c;
  int x0, y0, z0, x1, y1, z1;
  double endox, endoy, endoz;
  int i, indx, indy, indz;
  double origx, origy, origz;
 
  /* forward transform skewed -> rectangular (ax,ay,az,bx,by,bz,cx,cy,cz) */
  /* compute unit cell vectors a b c by intersection of projections in a,b plane; Bronstein & Semendjajew 2.6.6.1 */
  ax = 1.0; ay = 0.0; az = 0.0;
  bx = cos(PI*gamma/180.0); by = sin(PI*gamma/180.0); bz = 0.0;
  t1x = cos(PI*(gamma+alpha)/180.0); t1y = sin(PI*(gamma+alpha)/180.0);
  t2x = cos(PI*(gamma-alpha)/180.0); t2y = sin(PI*(gamma-alpha)/180.0);
  t3x = cos(PI*beta/180.0); t3y = sin(PI*beta/180.0);
  t4x = cos(PI*beta/180.0); t4y = -1.0 * sin(PI*beta/180.0);  
  cdet = (t4y-t3y)*(t2x-t1x)-(t2y-t1y)*(t4x-t3x);
  if(fabs(cdet)<1E-15) {
    error_divide_zero(71330, "map2map");
  }
  cx = ((t4x-t3x)*(t1y*t2x-t2y*t1x)-(t2x-t1x)*(t3y*t4x-t4y*t3x))/cdet;
  cy = ((t4y-t3y)*(t1y*t2x-t2y*t1x)-(t2y-t1y)*(t3y*t4x-t4y*t3x))/cdet;
  scz = 1.0 - (cx*cx+cy*cy);
  if(scz < 0.0) {
    error_sqrt_negative(71340, "map2map");
  }
  cz = sqrt(scz);
  
  /* inverse transform rectangular -> skewed (aix,aiy,aiz,bix,biy,biz,cix,ciy,ciz) */
  aix = 1.0; aiy = 0.0; aiz = 0.0;
  bix = -bx / by; biy = 1.0 / by; biz = 0.0;
  cix = (bx * cy - cx * by) / (by * cz); ciy = -cy / (by * cz); ciz = 1.0 / cz;

  /* assign origin and map it to skewed coordinates if it is not yet skewed */
  if (set_origin_get_mode(nxstart,nystart,nzstart,widthx,widthy,widthz,xorigin,yorigin,zorigin,&origx,&origy,&origz)) { 
    gx = aix * origx + bix * origy + cix * origz;
    gy = aiy * origx + biy * origy + ciy * origz;
    gz = aiz * origx + biz * origy + ciz * origz;
    origx = gx; origy = gy; origz = gz;
  }

  /* compute actual x y z extent of the skewed map: */
  endox = origx + (nx-1) * widthx;
  endoy = origy + (ny-1) * widthy;
  endoz = origz + (nz-1) * widthz;
  ux[0] = ax * origx + bx * origy + cx * origz;
  uy[0] = ay * origx + by * origy + cy * origz;
  uz[0] = az * origx + bz * origy + cz * origz;
  ux[1] = ax * endox + bx * origy + cx * origz;
  uy[1] = ay * endox + by * origy + cy * origz;
  uz[1] = az * endox + bz * origy + cz * origz;
  ux[2] = ax * origx + bx * endoy + cx * origz;
  uy[2] = ay * origx + by * endoy + cy * origz;
  uz[2] = az * origx + bz * endoy + cz * origz;
  ux[3] = ax * origx + bx * origy + cx * endoz;
  uy[3] = ay * origx + by * origy + cy * endoz;
  uz[3] = az * origx + bz * origy + cz * endoz;
  ux[4] = ax * endox + bx * endoy + cx * origz;
  uy[4] = ay * endox + by * endoy + cy * origz;
  uz[4] = az * endox + bz * endoy + cz * origz;
  ux[5] = ax * origx + bx * endoy + cx * endoz;
  uy[5] = ay * origx + by * endoy + cy * endoz;
  uz[5] = az * origx + bz * endoy + cz * endoz;
  ux[6] = ax * endox + bx * origy + cx * endoz;
  uy[6] = ay * endox + by * origy + cy * endoz;
  uz[6] = az * endox + bz * origy + cz * endoz;
  ux[7] = ax * endox + bx * endoy + cx * endoz;
  uy[7] = ay * endox + by * endoy + cy * endoz;
  uz[7] = az * endox + bz * endoy + cz * endoz;
  uxmin = 1E20; for(i=0;i<8;i++) if (ux[i] < uxmin) uxmin = ux[i];
  uymin = 1E20; for(i=0;i<8;i++) if (uy[i] < uymin) uymin = uy[i];
  uzmin = 1E20; for(i=0;i<8;i++) if (uz[i] < uzmin) uzmin = uz[i];
  uxmax = -1E20; for(i=0;i<8;i++) if (ux[i] > uxmax) uxmax = ux[i];
  uymax = -1E20; for(i=0;i<8;i++) if (uy[i] > uymax) uymax = uy[i];
  uzmax = -1E20; for(i=0;i<8;i++) if (uz[i] > uzmax) uzmax = uz[i];

  /* for interpolation, find minimum voxel spacing, rounded to nearest 0.1 Angstrom for good measure */
  *pwidth = widthx; /* ax = 1 */
  if ((widthy*by) < *pwidth) *pwidth = widthy*by; 
  if ((widthz*cz) < *pwidth) *pwidth = widthz*cz;
  *pwidth = floor(*pwidth*10.0+0.5)/10.0;

  /* compute output origin */
  /* we start pphiout at or below lower bound of phiin, and assert the new origin is in register with origin of the orthogonal coordinate system */ 
  *porigx = *pwidth * floor(uxmin / *pwidth);
  *porigy = *pwidth * floor(uymin / *pwidth);
  *porigz = *pwidth * floor(uzmin / *pwidth);

  /* compute output map dimensions */
  /* we end pphiout at or above upper bound of phiin */
  *pextx = (int) (ceil(uxmax/ *pwidth) - floor(uxmin/ *pwidth) + 1.5); 
  *pexty = (int) (ceil(uymax/ *pwidth) - floor(uymin/ *pwidth) + 1.5);
  *pextz = (int) (ceil(uzmax/ *pwidth) - floor(uzmin/ *pwidth) + 1.5);
  pnvox = *pextx * *pexty * *pextz;
   
  /* create output map, and loop through its cubic lattice to perform interpolation */
  do_vect(pphiout,pnvox);
  for (indz=0;indz<*pextz;indz++)
    for (indy=0;indy<*pexty;indy++)
      for (indx=0;indx<*pextx;indx++) {
	/* determine position in orthogonal coordinates */
	xpos = *porigx + indx * *pwidth;
	ypos = *porigy + indy * *pwidth; 
	zpos = *porigz + indz * *pwidth; 
	
	/* compute position of probe cube within skewed map in voxel units */
	gx = (aix * xpos + bix * ypos + cix * zpos - origx) / widthx;
	gy = (aiy * xpos + biy * ypos + ciy * zpos - origy) / widthy;
	gz = (aiz * xpos + biz * ypos + ciz * zpos - origz) / widthz;
	x0 = floor (gx); 
	y0 = floor (gy); 
	z0 = floor (gz); 
	x1 = x0+1;
	y1 = y0+1; 
	z1 = z0+1;
	
	/* if probe cube is fully within skewed map, do interpolate */
	if (x0>=0 && x1<nx && y0>=0 && y1<ny && z0>=0 && z1<nz) {	
	  a = gx-x0;
	  b = gy-y0;
	  c = gz-z0;
	  *(*pphiout+gidz_general(indz,indy,indx,*pexty,*pextx)) =
	    a * b * c * *(phiin+gidz_general(z1,y1,x1,ny,nx)) +
	    (1-a) * b * c * *(phiin+gidz_general(z1,y1,x0,ny,nx)) +
	    a * (1-b) * c * *(phiin+gidz_general(z1,y0,x1,ny,nx)) +
	    a * b * (1-c) * *(phiin+gidz_general(z0,y1,x1,ny,nx)) +
	    a * (1-b) * (1-c) * *(phiin+gidz_general(z0,y0,x1,ny,nx)) +
	    (1-a) * b * (1-c) * *(phiin+gidz_general(z0,y1,x0,ny,nx)) +
	    (1-a) * (1-b) * c * *(phiin+gidz_general(z1,y0,x0,ny,nx)) +
	    (1-a) * (1-b) * (1-c) * *(phiin+gidz_general(z0,y0,x0,ny,nx));
	}
      }	 
  printf("map2map> Conversion to cubic lattice completed. \n");
}



/* write X-PLOR map to vol_file */
static void write_xplor (char *vol_file, double pwidth, double porigx, double porigy,
			 double porigz, unsigned pextx, unsigned pexty, unsigned pextz, double *pphi) {
  FILE *fout; 
  long mxstart, mystart, mzstart, mxend, myend, mzend;
  unsigned indx, indy, indz;
  unsigned long count;
  unsigned extx2, exty2, extz2;
  double origx2, origy2, origz2;
  double *pphi2;
  
  fout = fopen(vol_file, "w");
  if(fout == NULL) {
    error_open_filename(71410, "map2map", vol_file);
  }
  
  /* X-PLOR does not support free origin definitions, must work within indexing constraints */
  /* if necessary, bring into register with coordinate system origin (for integer start indices) */
  if (test_registration(porigx,porigy,porigz,pwidth)==0) {
    printf("map2map> Input grid not in register with origin of coordinate system.\n");
    printf("map2map> Data will be interpolated to fit crystallographic X-PLOR format.\n");
  }
  interpolate_map (&pphi2, &extx2, &exty2, &extz2, &origx2, &origy2, &origz2,
		   pwidth, pwidth, pwidth, pphi, pextx, pexty, pextz, porigx,
		   porigy, porigz, pwidth, pwidth, pwidth);
  
  /* compute indices */
  mxstart = floor((origx2/pwidth)+0.5);
  mystart = floor((origy2/pwidth)+0.5);
  mzstart = floor((origz2/pwidth)+0.5);
  mxend = mxstart + extx2 - 1;
  myend = mystart + exty2 - 1;
  mzend = mzstart + extz2 - 1;
  
  /* write map */
  printf("map2map>\n");
  printf("map2map> Writing X-PLOR formatted (ASCII) volumetric map \n");
  fprintf(fout, " \n");
  fprintf(fout, " 2 !NTITLE \n");
  fprintf(fout, "REMARKS FILENAME=\"%s\" \n",vol_file);
  fprintf(fout, "REMARKS created by the Situs map2map utility\n");
  fprintf(fout, "%8d%8ld%8ld%8d%8ld%8ld%8d%8ld%8ld\n",extx2,mxstart,mxend,exty2,mystart,myend,extz2,mzstart,mzend);
  fprintf(fout, "%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\n",extx2*pwidth,exty2*pwidth,extz2*pwidth,90.0,90.0,90.0);
  fprintf(fout, "ZYX\n");
  for(indz=0;indz<extz2;indz++) {
    fprintf(fout, "%8d\n",indz);
    count = 0;
    for(indy=0;indy<exty2;indy++) for(indx=0;indx<extx2;indx++) {
      if ((count+1)%6 == 0) fprintf (fout,"%12.5E \n",*(pphi2+gidz_general(indz,indy,indx,exty2,extx2)));
      else fprintf (fout,"%12.5E",*(pphi2+gidz_general(indz,indy,indx,exty2,extx2))); 
      ++count;
    }
    if ((count)%6 != 0) fprintf (fout," \n");
  }
  fprintf(fout, "%8d\n",-9999);
  fclose(fout);
  
  /* print some info */
  printf("map2map> X-PLOR map written to file %s \n",vol_file);
  printf("map2map> X-PLOR map indexing (counting from 0): \n");
  printf("map2map>       NA = %8d  (# of X intervals in unit cell) \n",extx2);
  printf("map2map>     AMIN = %8ld  (start index X) \n",mxstart);
  printf("map2map>     AMAX = %8ld  (end index X) \n",mxend);
  printf("map2map>       NB = %8d  (# of Y intervals in unit cell) \n",exty2);
  printf("map2map>     BMIN = %8ld  (start index Y) \n",mystart);
  printf("map2map>     BMAX = %8ld  (end index Y) \n",myend);
  printf("map2map>       NC = %8d  (# of Z intervals in unit cell) \n",extz2);
  printf("map2map>     CMIN = %8ld  (start index Z) \n",mzstart);
  printf("map2map>     CMAX = %8ld  (end index Z) \n",mzend);
  printf("map2map> X-PLOR unit cell (based on the extent of Situs map): \n");
  printf("map2map>        A = %8.3f  (unit cell dimension) \n",extx2*pwidth);
  printf("map2map>        B = %8.3f  (unit cell dimension) \n",exty2*pwidth);
  printf("map2map>        C = %8.3f  (unit cell dimension) \n",extz2*pwidth);
  printf("map2map>    ALPHA = %8.3f  (unit cell angle) \n",90.0);
  printf("map2map>     BETA = %8.3f  (unit cell angle) \n",90.0);
  printf("map2map>    GAMMA = %8.3f  (unit cell angle) \n",90.0);
  printf("map2map> \n");  
  printf("map2map> All done.\n");
}



/* write MRC2000 / CCP4 map to vol_file in auto or manual mode */
static void write_mrc2000 (int automode, char *vol_file, double pwidth, double porigx, double porigy,
			   double porigz, unsigned pextx, unsigned pexty, unsigned pextz, 
			   double *pphi) {
  FILE *fout; 
  long nxstart, nystart, nzstart;
  float xorigin, yorigin, zorigin;
  unsigned long count, pnvox; 
  long nx, ny, nz, mx, my, mz;
  long mode;
  float xlen, ylen, zlen, alpha, beta, gamma;
  long mapc, mapr, maps;
  float dummy, fispg;
  double dmax, dmin, dav;
  float fmax, fmin, fav;
  char dummychar;
  char mapstring[4];
  int i, wmaperr;
	
  fout = fopen(vol_file, "wb");
  if( fout == NULL ) {
    error_open_filename(71420, "map2map", vol_file);
  }
  
  /* set initial values of map parameters */
  pnvox = pextx * pexty * pextz;
  nx = pextx; ny = pexty; nz = pextz; 
  fispg = 1.0f;
  dummy = 0.0f;
  dummychar = '\0';
  mapstring[0]  = 'M';
  mapstring[1]  = 'A';
  mapstring[2]  = 'P';
  mapstring[3]  = ' ';
  calc_map_info_short (pphi, pnvox, &dmax, &dmin, &dav);
  fmax = dmax; fmin = dmin; fav = dav;
  xorigin = porigx;
  yorigin = porigy;
  zorigin = porigz;
  mode = 2;
  mx = pextx; my = pexty; mz = pextz;
  xlen = pwidth * mx; ylen = pwidth * my; zlen = pwidth * mz;
  alpha = 90.0f; beta = 90.0f; gamma = 90.0f;
  mapc = 1; mapr = 2; maps = 3;

  /* if grid is not in register with origin we set the CCP4 start indices to zero */
  if (test_registration(porigx,porigy,porigz,pwidth)==0) {
    nxstart = 0; nystart = 0; nzstart = 0;
  } else {
    nxstart = floor((porigx/pwidth)+0.5);
    nystart = floor((porigy/pwidth)+0.5);
    nzstart = floor((porigz/pwidth)+0.5);
  }
  
  /* optional manual override of variable map parameters */
  if (automode==0) {
    printf("map2map> \n");
    printf("map2map> Manual override of most MRC2000 / CCP4 header fields for debugging or testing.\n");
    printf("map2map> With the exception of the MODE field, the written data is not affected by manually editing header values.\n");
    printf("map2map> \n");
    printf("map2map> The currently assigned MODE field is %ld \n",mode);
    printf("map2map> Pick one of the following two data modes: \n");
    printf("map2map>      0: density stored as signed 1-byte chars \n");
    printf("map2map>      2: density stored as 4-byte floats \n");
    printf("map2map> \n");
    printf("map2map> Enter either 0 or 2: ");
    mode = readln_int();
    if (mode != 0 && mode != 2) {
      mode = 2;
      printf("map2map> Did not recognize value, assuming MODE field of %ld \n",mode);
    }
    printf("map2map> The currently assigned CCP4-style NC field (# columns) is %ld \n",nx);
    printf("map2map> Enter the same or a new value: ");
    nx = readln_int();
    printf("map2map> The currently assigned CCP4-style NR field (# rows) is %ld \n",ny);
    printf("map2map> Enter the same or a new value: ");
    ny = readln_int();
    printf("map2map> The currently assigned CCP4-style NS field (# sections) is %ld \n",nz);
    printf("map2map> Enter the same or a new value: ");
    nz = readln_int();
    if (nx * ny * nz != pnvox) {
      printf("map2map> \n");
      fprintf (stderr,"map2map> Warning: NC * NR * NS does not match the total number of voxels in the map.\n");
      printf("map2map> \n");
    }
    printf("map2map> The currently assigned CCP4-style NCSTART field (column start index) is %ld \n",nxstart);
    if (nxstart == 0 && test_registration(porigx,porigy,porigz,pwidth)==0) {
      printf("map2map> Set to zero by default because input grid not in register with origin of coordinate system.\n");
    }
    printf("map2map> Enter the same or a new value: ");
    nxstart = readln_int();
    printf("map2map> The currently assigned CCP4-style NRSTART field (row start index) is %ld \n",nystart);
    if (nxstart == 0 && test_registration(porigx,porigy,porigz,pwidth)==0) {
      printf("map2map> Set to zero by default because input grid not in register with origin of coordinate system.\n");
    }
    printf("map2map> Enter the same or a new value: ");
    nystart = readln_int();
    printf("map2map> The currently assigned CCP4-style NSSTART field (section start index) is %ld \n",nzstart);
    if (nxstart == 0 && test_registration(porigx,porigy,porigz,pwidth)==0) {
      printf("map2map> Set to zero by default because input grid not in register with origin of coordinate system.\n");
    }
    printf("map2map> Enter the same or a new value: ");
    nzstart = readln_int();
    printf("map2map> The currently assigned MX field (# of X intervals in the unit cell) is %ld \n",mx);
    printf("map2map> Enter the same or a new value: ");
    mx = readln_int();
    printf("map2map> The currently assigned MY field (# of Y intervals in the unit cell) is %ld \n",my);
    printf("map2map> Enter the same or a new value: ");
    my = readln_int();
    printf("map2map> The currently assigned MZ field (# of Z intervals in the unit cell) is %ld \n",mz);
    printf("map2map> Enter the same or a new value: ");
    mz = readln_int();
    printf("map2map> The currently assigned X unit cell dimension is %f Angstrom\n",xlen);
    printf("map2map> Enter the same or a new value: ");
    xlen = readln_double();
    printf("map2map> The currently assigned Y unit cell dimension is %f Angstrom\n",ylen);
    printf("map2map> Enter the same or a new value: ");
    ylen = readln_double();
    printf("map2map> The currently assigned Z unit cell dimension is %f Angstrom\n",zlen);
    printf("map2map> Enter the same or a new value: ");
    zlen = readln_double();
    printf("map2map> The currently assigned unit cell angle alpha is %f degrees \n",alpha);
    printf("map2map> Enter the same or a new value: ");
    alpha = readln_double();
    printf("map2map> The currently assigned unit cell angle beta is %f degrees \n",beta);
    printf("map2map> Enter the same or a new value: ");
    beta = readln_double();
    printf("map2map> The currently assigned unit cell angle gamma is %f degrees \n",gamma);
    printf("map2map> Enter the same or a new value: ");
    gamma = readln_double();
    printf("map2map> The currently assigned MAPC field (axis order) is %ld \n",mapc);
    printf("map2map> Enter the same or a new value: ");
    mapc = readln_int();
    wmaperr=1;
    if (mapc == 1 || mapc == 2 || mapc == 3) wmaperr = 0;
    if (wmaperr) {
      mapc = 1; mapr = 2; maps = 3;
      printf("map2map> Inconsistent axis order values, assuming mapc = 1, mapr = 2, maps = 3. \n");
    } else {
      printf("map2map> The currently assigned MAPR field (axis order) is %ld \n",mapr);
      printf("map2map> Enter the same or a new value: ");
      mapr = readln_int();
      wmaperr=1;
      if (mapc==1 && mapr==2) wmaperr = 0;
      if (mapc==1 && mapr==3) wmaperr = 0;
      if (mapc==2 && mapr==1) wmaperr = 0;
      if (mapc==2 && mapr==3) wmaperr = 0;
      if (mapc==3 && mapr==1) wmaperr = 0;
      if (mapc==3 && mapr==2) wmaperr = 0;
      if (wmaperr) {
	mapc = 1; mapr = 2; maps = 3;
	printf("map2map> Inconsistent axis order values, assuming mapc = 1, mapr = 2, maps = 3. \n");
      } else {
	printf("map2map> The currently assigned MAPS field (axis order) is %ld \n",maps);
	printf("map2map> Enter the same or a new value: ");
	maps = readln_int();
	wmaperr=1;
	if (mapc==1 && mapr==2 && maps==3) wmaperr = 0;
	if (mapc==1 && mapr==3 && maps==2) wmaperr = 0;
	if (mapc==2 && mapr==1 && maps==3) wmaperr = 0;
	if (mapc==2 && mapr==3 && maps==1) wmaperr = 0;
	if (mapc==3 && mapr==1 && maps==2) wmaperr = 0;
	if (mapc==3 && mapr==2 && maps==1) wmaperr = 0;
	if (wmaperr) {
	  mapc = 1; mapr = 2; maps = 3;
	  printf("map2map> Inconsistent axis order values, assuming mapc = 1, mapr = 2, maps = 3. \n");
	}
      }
    }
    printf("map2map> The currently assigned MRC2000-style XORIGIN field is %f Angstrom\n",xorigin);
    printf("map2map> Enter the same or a new value: ");
    xorigin = readln_double();
    printf("map2map> The currently assigned MRC2000-style YORIGIN field is %f Angstrom\n",yorigin);
    printf("map2map> Enter the same or a new value: ");
    yorigin = readln_double();
    printf("map2map> The currently assigned MRC2000-style ZORIGIN field is %f Angstrom\n",zorigin);
    printf("map2map> Enter the same or a new value: ");
    zorigin = readln_double();
  }

  /* renormalize and threshold map if necessary */
  if (mode == 0) {
    if (fmax < 10.0 || fmax > 127.0) {      
      printf ("map2map> Density values will be rescaled to a max value of 127 to fit selected char data format.\n");
      printf ("map2map> Scaling factor: %f .\n", 127.0/fmax);
      normalize (pphi, pnvox, fmax/127.0);
      calc_map_info_short (pphi, pnvox, &dmax, &dmin, &dav);
      fmax = dmax; fmin = dmin; fav = dav;
    }
    if (fmin < -128.0) {
      printf ("map2map> Density values below -128 will be ignored to fit selected char data format.\n");
      threshold (pphi, pnvox, -128.0);
      calc_map_info_short (pphi, pnvox, &dmax, &dmin, &dav);
      fmax = dmax; fmin = dmin; fav = dav;
    }
  }
    
  /* write header */    
  printf("map2map>\n");
  printf("map2map> Writing MRC2000 / CCP4 (binary) volumetric map \n");
  wmaperr=0;
  if (fwrite(&nx, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&ny, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&nz, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&mode, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&nxstart, 4, 1, fout)!=1) wmaperr=1; 
  if (fwrite(&nystart, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&nzstart, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&mx, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&my, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&mz, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&xlen, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&ylen, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&zlen, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&alpha, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&beta, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&gamma, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&mapc, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&mapr, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&maps, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&fmin, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&fmax, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&fav, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&fispg, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&dummy, 4, 1, fout)!=1) wmaperr=1;
  for (i=25;i<50;++i) if (fwrite(&dummy, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&xorigin, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&yorigin, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&zorigin, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(mapstring, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&dummy, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&dummy, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&dummy, 4, 1, fout)!=1) wmaperr=1;
  for (i=0;i<800;++i) if (fwrite(&dummychar, sizeof(char), 1, fout)!=1) wmaperr=1;

  /* write actual data */  
  switch (mode) {
  case 0:
    for (count=0;count<pnvox;count++) {
      dummy = floor(*(pphi+count)+0.5);
      if (dummy < -128 || dummy > 127) wmaperr=1;
      dummychar = dummy;
      if (fwrite(&dummychar, 1, 1, fout)!=1) wmaperr=1;
    }
    break;
  case 2:
    for (count=0;count<pnvox;count++) {
      dummy = *(pphi+count);
      if (fwrite(&dummy, 4, 1, fout)!=1) wmaperr=1;
    }
    break;
  default:
    wmaperr=1;
    break;
  }
  if (wmaperr != 0) {
    error_write_filename(71430, "map2map");
  }
  fclose(fout);

  /* print some info */
  printf("map2map> Volumetric map written to file %s \n",vol_file);
  printf("map2map> Header information: \n");
  printf("map2map>       NC = %8ld  (# columns)\n",nx); 
  printf("map2map>       NR = %8ld  (# rows)\n",ny); 
  printf("map2map>       NS = %8ld  (# sections)\n",nz); 
  printf("map2map>     MODE = %8ld  (data type: 0: 1-byte char, 2: 4-byte float)\n",mode); 
  printf("map2map>  NCSTART = %8ld  (index of first column, counting from 0)\n",nxstart); 
  printf("map2map>  NRSTART = %8ld  (index of first row, counting from 0)\n",nystart);
  printf("map2map>  NSSTART = %8ld  (index of first section, counting from 0)\n",nzstart);
  printf("map2map>       MX = %8ld  (# of X intervals in unit cell)\n",mx); 
  printf("map2map>       MY = %8ld  (# of Y intervals in unit cell)\n",my);
  printf("map2map>       MZ = %8ld  (# of Z intervals in unit cell)\n",mz);
  printf("map2map> X length = %8.3f  (unit cell dimension)\n",xlen);
  printf("map2map> Y length = %8.3f  (unit cell dimension)\n",ylen);
  printf("map2map> Z length = %8.3f  (unit cell dimension)\n",zlen);
  printf("map2map>    Alpha = %8.3f  (unit cell angle)\n",alpha);
  printf("map2map>     Beta = %8.3f  (unit cell angle)\n",beta);
  printf("map2map>    Gamma = %8.3f  (unit cell angle)\n",gamma);
  printf("map2map>     MAPC = %8ld  (columns axis: 1=X,2=Y,3=Z)\n",mapc); 
  printf("map2map>     MAPR = %8ld  (rows axis: 1=X,2=Y,3=Z)\n",mapr);
  printf("map2map>     MAPS = %8ld  (sections axis: 1=X,2=Y,3=Z)\n",maps);
  printf("map2map>     DMIN = %8.3f  (minimum density value)\n",fmin);
  printf("map2map>     DMAX = %8.3f  (maximum density value)\n",fmax);
  printf("map2map>    DMEAN = %8.3f  (mean density value)\n",fav);
  printf("map2map>     ISPG = %8.3f  (space group number)\n",fispg);
  printf("map2map>   NSYMBT = %8d  (# bytes used for storing symmetry operators)\n",0);
  printf("map2map>  XORIGIN = %8.3f  (X origin - MRC2000 only)\n",xorigin);
  printf("map2map>  YORIGIN = %8.3f  (Y origin - MRC2000 only)\n",yorigin);
  printf("map2map>  ZORIGIN = %8.3f  (Z origin - MRC2000 only)\n",zorigin);
  printf("map2map> \n");

  /* if auto mode and grid is not in register with origin, issue a warning */
  if (automode==1 && test_registration(porigx,porigy,porigz,pwidth)==0) {
    printf("map2map>\n");
    printf("map2map> Input grid not in register with origin of coordinate system.\n");
    printf("map2map> The origin information was saved only in the MRC2000 format fields.\n");
    printf("map2map> The CCP4 start indexing was set to zero.\n");
    printf("map2map> To invoke a compatible (crystallographic) CCP4 indexing, you can force an interpolation \n");
    printf("map2map> as follows: Situs -> X-PLOR (forces an interpolation) -> Situs -> MRC2000/CCP4. \n");
    printf("map2map> \n");  
  } 
  printf("map2map> All done.\n");
}


/* write SPIDER map to vol_file */
static void write_spider (char *vol_file, double pwidth, double porigx, double porigy,
			  double porigz, unsigned pextx, unsigned pexty, unsigned pextz,
			  double *pphi) {
  FILE *fout; 
  unsigned long count, pnvox; 
  float nslice, nrow, iform, imami, nsam, headrec;
  double dmax, dmin, dav, dsig;
  float fmax, fmin, fav, fsig;
  float iangle, phi, theta, xoff, yoff, zoff, scale, labbyt, lenbyt, irec;
  float istack, inuse, maxim, kangle, phi1, theta1, psi1, phi2, theta2, psi2;
  float gamma;
  float dummy;
  int i, headlen, wmaperr;
  
  fout = fopen(vol_file, "wb");
  if( fout == NULL ) {
    error_open_filename(71460, "map2map", vol_file); 
  }
  
  /* set map parameters */
  pnvox = pextx * pexty * pextz;
  nsam = pextx; nrow = pexty; nslice = pextz; iform = 3.0f; imami = 1.0f; istack = 0.0f; inuse = -1.0f; maxim = 0.0f; kangle = 0.0f; 
  phi1 = 0.0f; theta1 = 0.0f; psi1 = 0.0f; phi2 = 0.0f; theta2 = 0.0f; psi2 = 0.0f; scale = 0.0f; dummy = 0.0f;
  phi = 0; theta = 0; gamma = 0; xoff = 0; yoff = 0; zoff = 0; iangle = 0;
  
  /* compute density info */
  calc_map_info (pphi, pnvox, &dmax, &dmin, &dav, &dsig);
  fmax = dmax; fmin = dmin; fav = dav; fsig = dsig;
  
  /* write map */
  headrec = ceil(256.0f/(pextx*1.0f)); 
  lenbyt = nsam * 4.0f;
  labbyt = headrec * lenbyt;
  irec = nslice * nrow + headrec;
  headlen = headrec * nsam;
  printf("map2map>\n");
  printf("map2map> Writing SPIDER (binary) volumetric map \n");
  wmaperr=0;	
  if (fwrite(&nslice, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&nrow, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&irec, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&dummy, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&iform, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&imami, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&fmax, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&fmin, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&fav, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&fsig, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&dummy, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&nsam, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&headrec, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&iangle, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&phi, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&theta, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&gamma, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&xoff, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&yoff, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&zoff, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&scale, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&labbyt, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&lenbyt, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&istack, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&inuse, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&maxim, 4, 1, fout)!=1) wmaperr=1;
  for (i=0; i<4; ++i) if (fwrite(&dummy, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&kangle, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&phi1, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&theta1, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&psi1, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&phi2, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&theta2, 4, 1, fout)!=1) wmaperr=1;
  if (fwrite(&psi2, 4, 1, fout)!=1) wmaperr=1;
  for (i=0;i<(headlen-37);++i) if (fwrite(&dummy, 4, 1, fout)!=1) wmaperr=1;
  for (count=0;count<pnvox;count++) {
    dummy = *(pphi+count);
    if (fwrite(&dummy, 4, 1, fout)!=1) wmaperr=1;
  }
  
  if (wmaperr != 0) {
    error_write_filename(71470, "map2map");
  }
  fclose(fout);
  
  /* print some info */
  printf("map2map>   SPIDER map written to file %s \n",vol_file);
  printf("map2map>   SPIDER map indexing: \n");
  printf("map2map>   NSLICE = %8.f  (# sections)\n",nslice); 
  printf("map2map>     NROW = %8.f  (# rows)\n",nrow); 
  printf("map2map>    IFORM = %8.f  (file type specifier)\n",iform); 
  printf("map2map>    IMAMI = %8.f  (flag: =1 the maximum and minimum values are computed)\n",imami); 
  printf("map2map>     FMAX = %8.3f  (maximum density value)\n",fmax); 
  printf("map2map>     FMIN = %8.3f  (minimum density value)\n",fmin);
  printf("map2map>       AV = %8.3f  (average density value)\n",fav);
  printf("map2map>      SIG = %8.3f  (standard deviation of density distribution)\n",fsig); 
  printf("map2map>     NSAM = %8.f  (# columns)\n",nsam);
  printf("map2map>  HEADREC = %8.f  (number of records in file header)\n",headrec);
  printf("map2map>   IANGLE = %8.f  (flag: =1 tilt angles filled)\n",iangle);
  printf("map2map>      PHI = %8.3f  (tilt angle)\n",phi);
  printf("map2map>    THETA = %8.3f  (tilt angle)\n",theta);
  printf("map2map>    GAMMA = %8.3f  (tilt angle)\n",gamma);
  printf("map2map>     XOFF = %8.3f  (X offset)\n",xoff);
  printf("map2map>     YOFF = %8.3f  (Y offset)\n",yoff);
  printf("map2map>     ZOFF = %8.3f  (Z offset)\n",zoff); 
  printf("map2map>    SCALE = %8.3f  (scale factor)\n",scale);
  printf("map2map>   LABBYT = %8.f  (total number of bytes in header)\n",labbyt);
  printf("map2map>   LENBYT = %8.f  (record length in bytes)\n",lenbyt);
  printf("map2map>   ISTACK = %8.f  (flag; file contains a stack of images)\n",istack);
  printf("map2map>    INUSE = %8.f  (flag; this image in stack is used)\n",inuse);
  printf("map2map>    MAXIM = %8.f  (maximum image used in stack)\n",maxim);
  printf("map2map>   KANGLE = %8.f  (flag; additional angles set)\n",kangle);
  printf("map2map>     PHI1 = %8.3f  (additional rotation)\n",phi1);
  printf("map2map>   THETA1 = %8.3f  (additional rotation)\n",theta1);
  printf("map2map>     PSI1 = %8.3f  (additional rotation)\n",psi1);
  printf("map2map>     PHI2 = %8.3f  (additional rotation)\n",phi2);
  printf("map2map>   THETA2 = %8.3f  (additional rotation)\n",theta2);
  printf("map2map>     PSI2 = %8.3f  (additional rotation)\n",psi2);
  printf("map2map> \n");
  printf("map2map> Warning: The Situs voxel spacing %f is not saved in the SPIDER map.\n",pwidth);
  printf("map2map> Warning: The Situs origin %f,%f,%f (first voxel position) is not saved in the SPIDER map.\n",porigx,porigy,porigz);
  printf("map2map> To keep maps properly aligned and scaled, avoid round-trip conversions SPIDER<->Situs\n");
  printf("map2map> \n");
  printf("map2map> All done.\n");
}


