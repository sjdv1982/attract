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

#ifndef CORE_TYPES_H_
#define CORE_TYPES_H_

#include "cuda_runtime.h"
#include <iostream>
#include <iomanip>
#include "config.h"
#include "as/interface.h"
#include "as/ParamTable.h"

namespace as {

/*
 ** @brief: Represents the volume spanned by 8 voxels.
 ** This is the basic entity interpolation is performed on.
 */
struct VoxelOctet {
	float4 data[2][2][2];	/** function values at the voxels */
	float3 min;				/** lower bound of the voxel coordinates */
	float3 max;				/** upper bound of the voxel coordinates */
};


struct DevFillLevel {
	unsigned id;
	unsigned level;

	inline bool operator< (const DevFillLevel& rhs) const {
		return level < rhs.level;
	}
};

typedef struct deviceIntrplGridDesc {
	unsigned width;		/** width of each grid in elements */
	unsigned height;   	/** height of each grid in elements */
	unsigned depth;   	/** depth of each grid in elements */
	float dVox;			/** voxel distance */
	float dVox_inv;		/** inverse voxel distance */
	float voxelVol_inv;	/** inverse voxel volume */
 	float3 minDim;		/** lower bound of grid dimensions */
 	float3 maxDim;		/** upper bound of grid dimensions */
 	cudaTextureObject_t*
		texArrayLin; 		/** texture objects for built-in interpolation  */
 	cudaTextureObject_t*
		texArrayPt; 		/** texture objects for manual interpolation  */

// 	deviceIntrplGridDesc () :
// 		width(0),height(0),	depth(0), dVox(0), dVox_inv(0),	voxelVol_inv(0),
// 		minDim(make_float3(0,0,0)),	maxDim(make_float3(0,0,0)),
// 		texArrayLin(NULL), texArrayPt(NULL) {}
} deviceIntrplGridDesc;

typedef struct hostIntrplGridResource {
	unsigned numArrays; 	/** number of cuArrays */
	/** Host pointers */
 	cudaArray** cuArrayPtr; /** Array of cudaArrays (host pointer)*/
 	cudaTextureObject_t* h_texArrayLin;
 	cudaTextureObject_t* h_texArrayPt;

 	/** Device pointers */
 	cudaTextureObject_t* d_texArrayLin;
 	cudaTextureObject_t* d_texArrayPt;


 	hostIntrplGridResource() : numArrays(0), cuArrayPtr(NULL),
 			h_texArrayLin(NULL), h_texArrayPt(NULL),
 			d_texArrayLin(NULL), d_texArrayPt(NULL){}
} hostIntrplGridResource;

typedef struct cudaIntrplGridDesc {
	deviceIntrplGridDesc deviceDesc;
	hostIntrplGridResource hostResc;

 	cudaIntrplGridDesc () : deviceDesc(), hostResc() {}

} cudaIntrplGridDesc;

typedef struct deviceNLGridDesc {
	unsigned width;		/** width of the grid in elements */
	unsigned height;   	/** height of the grid in elements */
	unsigned depth;   	/** depth of the grid in elements */
	float dVox_inv;		/** inverse voxel distance */
	float dPlateau2;     /** Plateau distance */
 	float3 minDim;		/** lower bound of grid dimensions */
 	float3 maxDim;		/** upper bound of grid dimensions */
 	cudaTextureObject_t tex; 	/** texture object */
 	unsigned* neighborList; /** device NL pointer*/

// 	deviceNLGridDesc () :
// 		width(0), height(0), depth(0), dVox_inv(0), dPlateau2(0),
// 		minDim(make_float3(0,0,0)),	maxDim(make_float3(0,0,0)),
// 		tex(), neighborList(NULL) {}

} deviceNLGridDesc;

typedef struct hostNLGridResource {
	/** Host pointers / objects*/
	cudaTextureObject_t tex;
	cudaArray* cuArray; /** cudaArray */

	hostNLGridResource() : tex(), cuArray(nullptr) {}
} hostNLGridResource;

typedef struct cudaNLGridDesc {
	deviceNLGridDesc deviceDesc;
	hostNLGridResource hostResc;

	cudaNLGridDesc() : deviceDesc(), hostResc() {}
} cudaNLGridDesc;

typedef struct deviceGridUnionDesc {
	deviceIntrplGridDesc inner;
	deviceIntrplGridDesc outer;
	deviceNLGridDesc NL;

} deviceGridUnionDesc;

typedef struct hostGridUnionResource {
	hostIntrplGridResource inner;
	hostIntrplGridResource outer;
	hostNLGridResource NL;

	hostGridUnionResource () : inner(), outer(), NL() {}
} hostGridUnionResource;

typedef struct cudaGridUnionDesc {
	deviceGridUnionDesc deviceDesc;
	hostGridUnionResource hostResc;

	cudaGridUnionDesc () : deviceDesc(), hostResc() {}
} cudaGridUnionDesc;


typedef struct deviceProteinDesc {
	unsigned nAtoms; // number of atoms of a single molecular conformer
	unsigned ntotAtoms; // total number of atoms of the protein (all conformers)

	float *xPos;	/** Cartesian coordinates in cm-frame*/
	float *yPos;
	float *zPos;

	unsigned* type; 	/** atom type */
	unsigned* mappedType;
	float* charge;	/** charge of the atoms/particle */

	unsigned numModes; /** number of modes */
	float* xModes; /** normal mode deformation vectors */
	float* yModes;
	float* zModes;

} deviceProteinDesc;

typedef deviceProteinDesc hostProteinResource;

typedef struct cudaProteinDesc {
	deviceProteinDesc deviceDesc;
	hostProteinResource hostResc;
} cudaProteinDesc;




struct deviceParamTableDesc {
	uint numTypes;  					/** number of particle/atom types */
	AttrParamTable::PotShape shape;		/** potential shape that is supported by the table */
 	AttrParamTable::type* paramTable; 	/** texture object */

 	inline __device__ const AttrParamTable::type getParams(int typeA, int typeB) const {
		return paramTable[numTypes*typeA + typeB];
	}
};


typedef deviceParamTableDesc hostParamTableResource;


typedef struct cudaParamTableDesc {
	deviceParamTableDesc deviceDesc;
	hostParamTableResource hostResc;
} cudaParamTableDesc;

template <unsigned T>
struct DOF_t {
	float3 pos;
	float3 ang;
	float modes[T];
	ushort conf;

	template <unsigned S>
	friend std::ostream& operator <<(std::ostream& outStream,
		const struct DOF_t<S> &dof);

	DOF_t() = default;
	DOF_t(float value) {
		conf = 0;
		pos = make_float3(value, value, value);
		ang = make_float3(value, value, value);
	}
};

template <unsigned T>
std::ostream& operator <<(std::ostream& outStream,
		const struct DOF_t<T> &dof)
{
	using namespace std;
	int precisionSetting = outStream.precision( );
	ios::fmtflags flagSettings = outStream.flags();
	outStream.setf(ios::scientific);
	outStream.precision(3);

	int w = 13;
	outStream 	<< setw(w) << "DOF"
				<< setw(w) << dof.pos.x << setw(w) << dof.pos.y << setw(w) << dof.pos.z
				<< setw(w) << dof.ang.x << setw(w) << dof.ang.y << setw(w) << dof.ang.z;

//	outStream 	<< dof.ang.x << "\t" << dof.ang.y << "\t" << dof.ang.z << "\t"
//				<< dof.pos.x << "\t" << dof.pos.y << "\t" << dof.pos.z;

	if(T > 0) {outStream << endl;}
	for (int i = 0; i < T; ++i) {
		outStream << setw(w) << dof.modes[i];
	}

	outStream.precision(precisionSetting);
	outStream.flags(flagSettings);

	return outStream;
}

template <unsigned T>
struct ServerDOF_t {
	int id; 		/** memory location for in order storage in result Buffer. */
	int locGridId;
	int locProtId;
	DOF_t<T> dof;
};

template <unsigned T>
std::ostream& operator <<(std::ostream& outStream,
		const struct ServerDOF_t<T> &dof)
{
	using namespace std;
	int precisionSetting = outStream.precision( );
	ios::fmtflags flagSettings = outStream.flags();
	outStream.precision(4);

	outStream << dof.id << " " << dof.locGridId << " " << dof.locProtId << endl;
	outStream << dof.dof << endl;

	outStream.precision(precisionSetting);
	outStream.flags(flagSettings);

	return outStream;
}

template <unsigned T>
struct ServerResult_t {
	float E_VdW;
	float E_El;
	float3 pos;
	float3 ang;
//	const unsigned numModes = T;
	float modes[T];
};

template <unsigned T>
std::ostream& operator <<(std::ostream& outStream,
		const struct ServerResult_t<T> &dof)
{
	using namespace std;
	int precisionSetting = outStream.precision( );
	ios::fmtflags flagSettings = outStream.flags();
	outStream.setf(ios::scientific | ios::showpos);
	outStream.precision(3);

	int w = 13;
	outStream << setw(w) << "Energy" << setw(w) << dof.E_VdW + dof.E_El  << setw(w) << dof.E_VdW << setw(w) << dof.E_El << endl;
	outStream 	<< setw(w) << "Gradients"
				<< setw(w) << dof.ang.x << setw(w) << dof.ang.y << setw(w) << dof.ang.z
				<< setw(w) << dof.pos.x << setw(w) << dof.pos.y << setw(w) << dof.pos.z;

//	outStream << dof.E_VdW + dof.E_El << "\t" << dof.E_VdW << "\t" << dof.E_El << "\t" << dof.ang.x << "\t" << dof.ang.y << "\t" << dof.ang.z << "\t"
//				<< dof.pos.x << "\t" << dof.pos.y << "\t" << dof.pos.z;

//	outStream << dof.E_VdW << "\t" << dof.E_El << "\t" << dof.ang.x << "\t" << dof.ang.y << "\t" << dof.ang.z << "\t"
//				<< dof.pos.x << "\t" << dof.pos.y << "\t" << dof.pos.z;

	if(T > 0) {outStream << endl;}
	for (int i = 0; i < T; ++i) {
		outStream << "ModeGrads\t"
				<< "\t" << dof.modes[i];
	}

	outStream.precision(precisionSetting);
	outStream.flags(flagSettings);

	return outStream;
}

typedef struct DOF_t<MAXMODES> DOF;
typedef struct ServerResult_t<MAXMODES> EnGrad;

struct deviceSimParam {
	dielec_t dielec;				/** type of dielectric constant */
	float epsilon;				/** dielectric constant */
	float ffelec;	/** precomputed factor felec/epsilon */
//	bool  useSwi;					/** using switching potential */
//	float swiOn;					/** min. switching potential distance */
//	float swiOff;					/** max. switching potential distance */
	bool useRecGrad;		/** using Receptor gradients */
	bool usePot;				/** use Potential grid */
};


} // namespace as

#endif /* CORE_TYPES_H_ */
