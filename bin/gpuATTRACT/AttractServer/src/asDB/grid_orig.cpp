/*******************************************************************************
 * gpuATTRACT framework
 * Copyright (C) 2015 Uwe Ehmann
 *
 * This file is part of the gpuATTRACT framework.
 *
 * The gpuATTRACT framework is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The gpuATTRACT framework is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
#include "grid_orig.h"
//#include "nonbon.h"
//#include "state.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include <fcntl.h>
#include <unistd.h>

#ifdef TORQUEGRID
#define EnerGradX EnerGradTorque
#define energradsX energrads_torque
#define IS_TORQUEGRID 1
#else
#define EnerGradX EnerGradStd
#define energradsX energrads_std
#define IS_TORQUEGRID 0
#endif

static void get_shm_name(int shm_id, char *shm_name) {
	sprintf(shm_name, "/attract-grid%d", shm_id);
}

void Grid::init(double gridspacing0, int gridextension0, double plateaudis0,
		double neighbourdis0, bool (&alphabet0)[MAXATOMTYPES]) {
	gridspacing = gridspacing0;
	gridextension = gridextension0;
	plateaudis = plateaudis0;
	plateaudissq = plateaudis * plateaudis;
	plateaudissqinv = 1.0 / plateaudissq;
	neighbourdis = neighbourdis0;
	neighbourdissq = neighbourdis * neighbourdis;
	architecture = ARCHITECTURE;
	//Pre-compute the scale-down-distance ratios
	int size_ratio = int(10000 * plateaudissq);
	_ratio = new double[size_ratio + 1];
	for (int n = 0; n <= size_ratio; n++) {
		double dissq = ((n + 0.5) / 10000);
		_ratio[n] = sqrt(dissq / plateaudissq);
	}
	memset(alphabet, 0, MAXATOMTYPES * sizeof(bool));
	memcpy(alphabet, alphabet0, sizeof(alphabet));
	alphabetsize = 0;
	for (int n = 0; n < MAXATOMTYPES; n++) {
		if (alphabet[n]) {
			alphabetsize++;
		}
	}
}

static void error(const char *filename) {
	fprintf(stderr, "Reading error in grid file %s\n", filename);
	exit(1);
}

extern void get_shm_name(int shm_id, char *shm_name);

#ifdef TORQUEGRID
void Grid::read_torque(const char *filename) {
#else
void Grid::read_std(const char *filename) {
#endif
	int n;

	int read;

	this->torquegrid = IS_TORQUEGRID;

	FILE *f = fopen(filename, "rb");
	if (f == NULL)
		error(filename);

	bool torquegrid0;
	read = fread(&torquegrid0, sizeof(bool), 1, f);
	if (IS_TORQUEGRID && (!torquegrid0)) {
		fprintf(stderr,
				"Reading error in grid file %s, grid was computed as a normal grid, but it was specified as a torque grid\n",
				filename);
		exit(1);
	}
	if ((!IS_TORQUEGRID) && torquegrid0) {
		fprintf(stderr,
				"Reading error in grid file %s, grid was computed as a torque grid, but it was specified as a normal grid\n",
				filename);
		exit(1);
	}
	read = fread(&architecture, sizeof(short), 1, f);
	if (architecture != ARCHITECTURE) {
		fprintf(stderr,
				"Reading error in grid file %s, grid was computed on %d bit, but we are on %d bit\n",
				filename, architecture, ARCHITECTURE);
		exit(1);
	}
	if (!read)
		error(filename);
	read = fread(&gridspacing, sizeof(double), 1, f);
	if (!read)
		error(filename);
	read = fread(&gridextension, sizeof(int), 1, f);
	if (!read)
		error(filename);
	read = fread(&plateaudis, sizeof(double), 1, f);
	if (!read)
		error(filename);
	read = fread(&neighbourdis, sizeof(double), 1, f);
	if (!read)
		error(filename);
	bool alphabet0[MAXATOMTYPES];
	read = fread(&alphabet0, sizeof(alphabet0), 1, f);
	if (!read)
		error(filename);
	init(gridspacing, gridextension, plateaudis, neighbourdis, alphabet0);

	float arr1[3];
	read = fread(arr1, 3 * sizeof(float), 1, f);
	if (!read)
		error(filename);
	ori[0] = arr1[0];
	ori[1] = arr1[1];
	ori[2] = arr1[2];
	int arr2[6];
	read = fread(arr2, 6 * sizeof(int), 1, f);
	if (!read)
		error(filename);
	gridx = arr2[0];
	gridy = arr2[1];
	gridz = arr2[2];
	gridx2 = arr2[3];
	gridy2 = arr2[4];
	gridz2 = arr2[5];
	read = fread(&natoms, sizeof(int), 1, f);
	if (!read)
		error(filename);
	read = fread(&pivot, sizeof(Coor), 1, f);
	if (!read)
		error(filename);

	read = fread(&nr_energrads, sizeof(nr_energrads), 1, f);
	if (!read)
		error(filename);
	read = fread(&shm_energrads, sizeof(shm_energrads), 1, f);
	if (!read)
		error(filename);

	if (nr_energrads) {
		if (shm_energrads == -1) {
			energradsX = new EnerGradX[nr_energrads];
			read = fread(energradsX, nr_energrads * sizeof(EnerGradX), 1, f);
			if (!read)
				error(filename);
		} else {
			char shm_name[100];
			get_shm_name(shm_energrads, shm_name);
			int fshm1 = shm_open(shm_name, O_RDONLY, S_IREAD);
			if (fshm1 == -1) {
				fprintf(stderr,
						"Reading error in grid file %s: shared memory segment %d for potential list does not exist\n",
						filename, shm_energrads);
				exit(1);
			}
			int success = ftruncate(fshm1, nr_energrads * sizeof(EnerGradX));
			if (success == -1) {
				fprintf(stderr, "Error ftruncate\n");
				exit(1);
			}
			energradsX = (EnerGradX *) mmap(0, nr_energrads * sizeof(EnerGradX),
					PROT_READ, MAP_SHARED | MAP_NORESERVE, fshm1, 0);
			if (energradsX == NULL) {
				fprintf(stderr,
						"Reading error in grid file %s: Could not load shared memory segment %d\n",
						filename, shm_energrads);
				exit(1);
			}
		}
	}

	read = fread(&nr_neighbours, sizeof(nr_neighbours), 1, f);
	if (!read)
		error(filename);
	read = fread(&shm_neighbours, sizeof(shm_neighbours), 1, f);
	if (!read)
		error(filename);
	if (shm_neighbours == -1) {
		neighbours = new Neighbour[nr_neighbours];
		read = fread(neighbours, nr_neighbours * sizeof(Neighbour), 1, f);
		if (!read)
			error(filename);
	} else {
		char shm_name[100];
		get_shm_name(shm_neighbours, shm_name);
		int fshm2 = shm_open(shm_name, O_RDONLY, S_IREAD);
		if (fshm2 == -1) {
			fprintf(stderr,
					"Reading error in grid file %s: shared memory segment %d for neighbour list does not exist\n",
					filename, shm_neighbours);
			exit(1);
		}
		int success = ftruncate(fshm2, nr_neighbours * sizeof(Neighbour));
		if (success == -1) {
			fprintf(stderr, "Error ftruncate\n");
			exit(1);
		}
		neighbours = (Neighbour *) mmap(0, nr_neighbours * sizeof(Neighbour),
				PROT_READ, MAP_SHARED | MAP_NORESERVE, fshm2, 0);
		if (neighbours == NULL) {
			fprintf(stderr,
					"Reading error in grid file %s: Could not load shared memory segment %d\n",
					filename, shm_neighbours);
			exit(1);
		}
	}
	long innergridsize, biggridsize;
	read = fread(&innergridsize, sizeof(innergridsize), 1, f);
	if (!read)
		error(filename);
	innergrid = new Voxel[innergridsize];
	read = fread(innergrid, innergridsize * sizeof(Voxel), 1, f);
	if (!read)
		error(filename);
	read = fread(&biggridsize, sizeof(biggridsize), 1, f);
	if (!read)
		error(filename);
	if (biggridsize) {
		biggrid = new Potential[biggridsize];
		read = fread(biggrid, biggridsize * sizeof(Potential), 1, f);
	}
	if (!read)
		error(filename);
	fclose(f);
	for (n = 0; n < innergridsize; n++) {
		if (innergrid[n].potential[MAXATOMTYPES]) {
			for (int i = 0; i <= MAXATOMTYPES; i++) {
				if (i < MAXATOMTYPES && alphabet[i] == 0)
					continue;
				int dif = innergrid[n].potential[i] - 1;
				if (dif < 0 || dif >= nr_energrads) {
					fprintf(stderr,
							"Reading error in %s, innergrid voxel %d atom type %d: %d >= %d\n",
							filename, n + 1, i + 1, dif, nr_energrads);
					exit(1);
				}
			}
			int nl = innergrid[n].neighbourlist;
			int nr = innergrid[n].nr_neighbours;
			bool empty1 = (nl == 0);
			bool empty2 = (nr == 0);

			if (empty1 != empty2 || (nl + nr - 1 > nr_neighbours)) {
				fprintf(stderr,
						"Reading error in %s, innergrid voxel %d neighbourlist: %d + %d >= %d\n",
						filename, n + 1, nl - 1, nr, nr_neighbours);
				exit(1);
			}

		}
	}
	for (n = 0; n < biggridsize; n++) {
		if (biggrid[n][MAXATOMTYPES]) {
			for (int i = 0; i <= MAXATOMTYPES; i++) {
				if (i < MAXATOMTYPES && alphabet[i] == 0)
					continue;
				int dif = biggrid[n][i] - 1;
				if (dif < 0 || dif >= nr_energrads) {
					fprintf(stderr,
							"Reading error in %s, biggrid voxel %d atom type %d: %d >= %d\n",
							filename, n + 1, i + 1, dif, nr_energrads);
					exit(1);
				}
			}
		}
	}

	init(gridspacing, gridextension, plateaudis, neighbourdis, alphabet0);
}

#ifdef TORQUEGRID
void Grid::write_torque(const char *filename) {
#else
void Grid::write_std(const char *filename) {
#endif
	long n;
	long innergridsize = gridx * gridy * gridz;
	long biggridsize = gridx2 * gridy2 * gridz2;
	if (!nr_energrads)
		biggridsize = 0;

	FILE *f = fopen(filename, "wb");
	if (f == NULL) {
		fprintf(stderr,
				"Grid::write error for %s: Cannot open file for writing\n",
				filename);
		exit(1);
	}

	EnerGradX *shmptr1 = NULL;
	Neighbour *shmptr2 = NULL;
	if (nr_energrads && shm_energrads != -1) {
		char shm_name[100];
		get_shm_name(shm_energrads, shm_name);
		int fshm1 = shm_open(shm_name, (O_CREAT | O_RDWR),
				(S_IREAD | S_IWRITE));
		if (fshm1 == -1) {
			fprintf(stderr,
					"Grid::write error for %s: Cannot open shared memory for writing\n",
					filename);
			exit(1);
		}
		int success = ftruncate(fshm1, nr_energrads * sizeof(EnerGradX));
		if (success == -1) {
			fprintf(stderr, "Error ftruncate\n");
			exit(1);
		}
		shmptr1 = (EnerGradX *) mmap(0, nr_energrads * sizeof(EnerGradX),
				(PROT_READ | PROT_WRITE), MAP_SHARED, fshm1, 0);
		if (shmptr1 == MAP_FAILED ) {
			fprintf(stderr,
					"Grid::write error for %s: Cannot map shared memory for writing\n",
					filename);
			exit(1);
		}
		memset(shmptr1, 0, nr_energrads * sizeof(EnerGradX));
		close(fshm1);
	}

	if (shm_neighbours != -1) {
		char shm_name[100];
		get_shm_name(shm_neighbours, shm_name);
		int fshm2 = shm_open(shm_name, (O_CREAT | O_RDWR),
				(S_IREAD | S_IWRITE));
		if (fshm2 == -1) {
			fprintf(stderr,
					"Grid::write error for %s: Cannot open shared memory for writing\n",
					filename);
			exit(1);
		}
		int success = ftruncate(fshm2, nr_neighbours * sizeof(Neighbour));
		if (success == -1) {
			fprintf(stderr, "Error ftruncate\n");
			exit(1);
		}
		shmptr2 = (Neighbour *) mmap(0, nr_neighbours * sizeof(Neighbour),
				(PROT_READ | PROT_WRITE), MAP_SHARED, fshm2, 0);
		if (shmptr2 == MAP_FAILED ) {
			fprintf(stderr,
					"Grid::write error for %s: Cannot map shared memory for writing\n",
					filename);
			exit(1);
		}
		memset(shmptr2, 0, nr_neighbours * sizeof(Neighbour));
		close(fshm2);
	}

	fwrite(&torquegrid, sizeof(bool), 1, f);
	fwrite(&architecture, sizeof(unsigned short), 1, f);
	fwrite(&gridspacing, sizeof(double), 1, f);
	fwrite(&gridextension, sizeof(int), 1, f);
	fwrite(&plateaudis, sizeof(double), 1, f);
	fwrite(&neighbourdis, sizeof(double), 1, f);
	fwrite(&alphabet, sizeof(alphabet), 1, f);
	float arr1[] = { (float) ori[0], (float) ori[1], (float) ori[2] };
	fwrite(arr1, 3 * sizeof(float), 1, f);
	int arr2[] = { gridx, gridy, gridz, gridx2, gridy2, gridz2 };
	fwrite(arr2, 6 * sizeof(int), 1, f);
	fwrite(&natoms, sizeof(int), 1, f);
	fwrite(&pivot, sizeof(Coor), 1, f);

	if (nr_energrads) {
		energradsX = (EnerGradX *) realloc(energradsX,
				nr_energrads * sizeof(EnerGradX));
		EnerGradX *energrads_reordered = new EnerGradX[nr_energrads];
		memset(energrads_reordered, 0, nr_energrads * sizeof(EnerGradX));
		if (energrads_reordered) { //only re-order the energrads if we got the memory for it
			int nr_energrads2 = 0;
			for (n = 0; n < innergridsize; n++) {
				if (innergrid[n].potential[MAXATOMTYPES]) {
					for (int i = 0; i <= MAXATOMTYPES; i++) {
						if (i < MAXATOMTYPES && alphabet[i] == 0)
							continue;
						unsigned int &oldpos = innergrid[n].potential[i];
						unsigned int newpos = nr_energrads2 + 1;
						memcpy(&energrads_reordered[newpos - 1],
								&energradsX[oldpos - 1], sizeof(EnerGradX));
						oldpos = newpos;
						nr_energrads2++;
					}
				}
			}
			for (n = 0; n < biggridsize; n++) {
				if (biggrid[n][MAXATOMTYPES]) {
					for (int i = 0; i <= MAXATOMTYPES; i++) {
						if (i < MAXATOMTYPES && alphabet[i] == 0)
							continue;
						unsigned int &oldpos = biggrid[n][i];
						unsigned int newpos = nr_energrads2 + 1;
						memcpy(&energrads_reordered[newpos - 1],
								&energradsX[oldpos - 1], sizeof(EnerGradX));
						oldpos = newpos;
						nr_energrads2++;
					}
				}
			}
			if (nr_energrads != nr_energrads2) {
				fprintf(stderr, "ERR nr_energrads %d %d\n", nr_energrads,
						nr_energrads2);
			}
			free(energradsX);
			energradsX = energrads_reordered;
		}
	}

	fwrite(&nr_energrads, sizeof(nr_energrads), 1, f);
	fwrite(&shm_energrads, sizeof(shm_energrads), 1, f);

	if (nr_energrads) {
		if (shm_energrads == -1)
			fwrite(energradsX, nr_energrads * sizeof(EnerGradX), 1, f);
		else
			memcpy(shmptr1, energradsX, nr_energrads * sizeof(EnerGradX));
	}

	if (nr_neighbours) {
		Neighbour *neighbours_reordered = new Neighbour[nr_neighbours];
		memset(neighbours_reordered, 0, nr_neighbours * sizeof(Neighbour));
		if (neighbours_reordered) { //only re-order the neighbours if we got the memory for it
			int nr_neighbours2 = 0;
			for (n = 0; n < innergridsize; n++) {
				Voxel &v = innergrid[n];
				if (v.nr_neighbours) {
					memcpy(&neighbours_reordered[nr_neighbours2],
							&neighbours[v.neighbourlist - 1],
							v.nr_neighbours * sizeof(Neighbour));
					v.neighbourlist = nr_neighbours2 + 1;
					nr_neighbours2 += v.nr_neighbours;
				}
			}
			if (nr_neighbours2 != nr_neighbours) {
				fprintf(stderr, "ERR nr_neighbours %d %d\n", nr_neighbours,
						nr_neighbours2);
			}
			delete[] neighbours;
			neighbours = neighbours_reordered;
		}
	}

	fwrite(&nr_neighbours, sizeof(nr_neighbours), 1, f);
	fwrite(&shm_neighbours, sizeof(shm_neighbours), 1, f);
	if (shm_neighbours == -1)
		fwrite(neighbours, nr_neighbours * sizeof(Neighbour), 1, f);
	else
		memcpy(shmptr2, neighbours, nr_neighbours * sizeof(Neighbour));

	fwrite(&innergridsize, sizeof(innergridsize), 1, f);
	fwrite(innergrid, innergridsize * sizeof(Voxel), 1, f);
	fwrite(&biggridsize, sizeof(biggridsize), 1, f);
	if (biggridsize) {
		fwrite(biggrid, biggridsize * sizeof(Potential), 1, f);
	}
	fclose(f);
}

