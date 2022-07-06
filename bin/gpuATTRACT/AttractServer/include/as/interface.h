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

#ifndef INTERFACE_H_
#define INTERFACE_H_

#include <cstdlib>
#include <vector_types.h>

#define MAXIDL 30

/*
 ** @brief: description of a protein. Position and normal mode deformation vectors
 ** are given according to the center of mass frame.
 ** Rotations are performed along the the center of mass frame (for now).
 */
struct ProteinDesc {
	char id[MAXIDL];	/** identifier */
	unsigned nAtoms; // number of atoms of a single molecular conformer
	unsigned ntotAtoms; // total number of atoms of the protein (all conformers)

	/*
	 ** @data layout:
	 ** pos
	 ** size: 3*ntotAtoms
	 ** {x0,x1,...,xN,y0,y1,...,yN,z0,z1,...zN}
	 ** N: number of atoms/particles
	 ** the first nAtoms values are for conformer 1, then for conformer 2, etc.
	 */
	float *pos;	/** Cartesian coordinates */

	/*
	 ** @data layout:
	 ** type & charge
	 ** size: ntotAtoms
	 ** {e1,e2,...,eN}
	 ** N: ntotAtoms
	 ** the first nAtoms values are for conformer 1, then for conformer 2, etc.
	 */
	unsigned* type; 	/** atom type */
	float* charge;	/** charge of the atoms/particle */

	unsigned numModes; /** number of modes */

	/*
	 ** @data layout:
	 ** modes:
	 ** size: 3*numModes*nAtoms
	 ** {x0_0,x0_1,...,x0_(M-1),  x1_0,x1_1,...,x1_(M-1),...,  x(N-1)_0,x(N-1)_1,...,xN_(M-1),
	 **  y0_0,y0_1,...,y0_(M-1),  y1_0,y1_1,...,y1_(M-1),...,  y(N-1)_0,y(N-1)_1,...,yN_(M-1),
	 **  z0_0,z0_1,...,z0_(M-1),  z1_0,z1_1,...,z1_(M-1),...,  z(N-1)_0,z(N-1)_1,...,zN_(M-1)}
	 ** x0_1: x-component of mode 1 of particle 0
	 ** M: number of modes
	 ** N: number of atoms
	 */
	float* modes; /** normal mode deformation vectors */

};


/*
 ** @brief: Descriptions of a gradient-energy-grid
 */
struct GradEnGridDesc {
	// ToDo: Remove typemask since it is very application/client specific! Instead use numGrids.
	bool typemask[99]; /** mask of atom types which are supported in that grid */
	unsigned numGrids;	/** number of grids. The last Grid is the charge grid */
	unsigned width;		/** number of elements along x */
	unsigned height;	/** number of elements along y */
	unsigned depth;		/** number of elements along z */

	float gridSpacing; 	/** grid spacing */
	float posMin[3];	/** lower coordinate bounds == position of grid[0] */


	/*
	 ** @data layout:
	 ** grid
	 ** size: numGrids*width*height*depth
	 ** {grid_0,grid_1,...,grid_(K-1)} with each grid defined as slices of 2D grids:
	 ** {s0,...,s(Nz-1)} with each slice defined like
	 ** {e00,...,e(Nx-1)0, e01,...,e(Nx-1)1,...,e(Nx-1)(Ny-1)}
	 ** K: number of grids
	 ** Nz: depth
	 ** Ny: height
	 ** Nx: width
	 */
	float4* grid; /** grid data. The first 3 entries of float4 contain the forces.
	 	 	 	 	  The last one contains the energy */
};

/*
 ** @brief: data that is contained by each grid element of a neighbourlist
 ** grid.
 */

typedef enum {
	constant,
	variable
} dielec_t;

struct NeighbourDesc {
	unsigned numEl;	/** number of elements in the neighbourlist.
	 	 	 	 If numEl == 0 then idx is ignored  */
	unsigned idx;	/** starting index in NLGridDesc.neighbourArray */
};


/*
 ** @brief: Descriptions of a neighbour list grid
 */
struct NLGridDesc {
	unsigned width;		/** number of elements along x */
	unsigned height;	/** number of elements along y */
	unsigned depth;		/** number of elements along z */

	float gridSpacing; 	/** grid spacing */
	float dPlateau; 	/** Plateau distance: cutOff */
	float posMin[3];	/** lower coordinate bounds == position of grid[0] */

	/*
	 ** @data layout:
	 ** grid
	 ** size: width*height*depth
	 ** slices of 2D grids {s0,...,s(Nz-1)} with each slice defined like
	 ** {e00,...,e(Nx-1)0, e01,...,e(Nx-1)1,...,e(Nx-1)(Ny-1)}
	 ** Nz: depth
	 ** Ny: height
	 ** Nx: width
	 */
	NeighbourDesc* grid;	/** grid of neighbour descriptions */

	unsigned numEl;				/** number of total elements in neighbourArray */

	/*
	 ** @data layout:
	 ** neighbourArray
	 ** size: numEl
	 ** {list_0,...,list_(N-1)} with each list
	 ** {idx_0,...,idx_Li-1} where Li is the length of the i-th list
	 ** The sum of all Li's equals numEl
	 ** N = width*height*depth
	 */
	unsigned* neighborArray;	/** contains receptor atom indices for all neighbours.
								It is basically a concatenation of subsequent neighbour lists.
								The exact sequence is the same as for the grid */
};


#endif /* INTERFACE_H_ */
