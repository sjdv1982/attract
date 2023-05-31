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

#include "cuda_runtime.h"
#include <cstring>
#include <cassert>

#include "as/DeviceDataFactory.h"
#include "asUtils/macros.h"


/****************************
 * public member functions
 ****************************/
as::cudaGridUnionDesc as::DeviceDataFactory::initDeviceGridUnion(const GridUnion* gridUnion, int deviceId)
{
	ASSERT(deviceId >= 0);
	CUDA_CHECK(cudaSetDevice(deviceId));


	cudaIntrplGridDesc inner = DeviceDataFactory::initIntrpl(gridUnion->innerGrid());
	cudaIntrplGridDesc outer = DeviceDataFactory::initIntrpl(gridUnion->outerGrid());
	cudaNLGridDesc NL = DeviceDataFactory::initNL(gridUnion->NLgrid());

	deviceGridUnionDesc deviceDesc;
	deviceDesc.inner = inner.deviceDesc;
	deviceDesc.outer = outer.deviceDesc;
	deviceDesc.NL = NL.deviceDesc;

	hostGridUnionResource hostResc;
	hostResc.inner = inner.hostResc;
	hostResc.outer = outer.hostResc;
	hostResc.NL = NL.hostResc;

	cudaGridUnionDesc cudaDesc;
	cudaDesc.deviceDesc = deviceDesc;
	cudaDesc.hostResc = hostResc;
	return cudaDesc;
}

as::cudaProteinDesc as::DeviceDataFactory::initDeviceProtein(const Protein* protein, int deviceId)
{
	ASSERT(deviceId >= 0);
	CUDA_CHECK(cudaSetDevice(deviceId));

	deviceProteinDesc deviceDesc;
	unsigned ntotAtoms = protein->ntotAtoms();
	unsigned nAtoms = protein->nAtoms();

	float *d_xPos;
	CUDA_CHECK(cudaMalloc((void**) &d_xPos, ntotAtoms * sizeof(float)));
	CUDA_CHECK(cudaMemcpy(d_xPos, protein->xPos(0), ntotAtoms * sizeof(float), cudaMemcpyHostToDevice));
	float *d_yPos;
	CUDA_CHECK(cudaMalloc((void**) &d_yPos, ntotAtoms * sizeof(float)));
	CUDA_CHECK(cudaMemcpy(d_yPos, protein->yPos(0), ntotAtoms * sizeof(float), cudaMemcpyHostToDevice));
	float *d_zPos;
	CUDA_CHECK(cudaMalloc((void**) &d_zPos, ntotAtoms * sizeof(float)));
	CUDA_CHECK(cudaMemcpy(d_zPos, protein->zPos(0), ntotAtoms * sizeof(float), cudaMemcpyHostToDevice));

	unsigned* d_type;
	CUDA_CHECK(cudaMalloc((void**) &d_type, nAtoms * sizeof(unsigned)));
	CUDA_CHECK(cudaMemcpy(d_type, protein->type(), nAtoms * sizeof(unsigned), cudaMemcpyHostToDevice));
	unsigned* d_mappedType;
	CUDA_CHECK(cudaMalloc((void**) &d_mappedType, nAtoms * sizeof(unsigned)));
	CUDA_CHECK(cudaMemcpy(d_mappedType, protein->mappedType(), nAtoms * sizeof(unsigned), cudaMemcpyHostToDevice));
	float* d_charge;
	CUDA_CHECK(cudaMalloc((void**) &d_charge, nAtoms * sizeof(float)));
	CUDA_CHECK(cudaMemcpy(d_charge, protein->charge(), nAtoms * sizeof(float), cudaMemcpyHostToDevice));

	unsigned numModes = protein->numModes();
	float* d_xModes = NULL;
	float* d_yModes = NULL;
	float* d_zModes = NULL;
	if (numModes != 0) {
		CUDA_CHECK(cudaMalloc((void**) d_xModes, nAtoms * numModes * sizeof(float)));
		CUDA_CHECK(cudaMemcpy(d_xModes, protein->xModes(), nAtoms * numModes * sizeof(float), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMalloc((void**) d_yModes, nAtoms * numModes * sizeof(float)));
		CUDA_CHECK(cudaMemcpy(d_yModes, protein->yModes(), nAtoms * numModes * sizeof(float), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMalloc((void**) d_zModes, nAtoms * numModes * sizeof(float)));
		CUDA_CHECK(cudaMemcpy(d_zModes, protein->zModes(), nAtoms * numModes * sizeof(float), cudaMemcpyHostToDevice));
	}

	deviceDesc.ntotAtoms = ntotAtoms;
	deviceDesc.nAtoms = nAtoms;
	deviceDesc.xPos = d_xPos;
	deviceDesc.yPos = d_yPos;
	deviceDesc.zPos = d_zPos;
	deviceDesc.type = d_type;
	deviceDesc.mappedType = d_mappedType;
	deviceDesc.charge = d_charge;
	deviceDesc.numModes = numModes;
	deviceDesc.xModes = d_xModes;
	deviceDesc.yModes = d_yModes;
	deviceDesc.zModes = d_zModes;

	cudaProteinDesc cudaDesc;
	cudaDesc.deviceDesc = deviceDesc;
	cudaDesc.hostResc = deviceDesc;

	return cudaDesc;
}

as::cudaParamTableDesc as::DeviceDataFactory::initDeviceParamTable (
		const AttrParamTable* table, int deviceId)
{
	ASSERT(deviceId >= 0);
	CUDA_CHECK(cudaSetDevice(deviceId));

	deviceParamTableDesc deviceDesc;

	const unsigned numTypes = table->numTypes();

	AttrParamTable::type* d_table;
	CUDA_CHECK(cudaMalloc((void**)&d_table, numTypes*numTypes*sizeof(AttrParamTable::type)));
	CUDA_CHECK(cudaMemcpy(d_table, table->table(), numTypes*numTypes*sizeof(AttrParamTable::type), cudaMemcpyHostToDevice));

	deviceDesc.numTypes = numTypes;
	deviceDesc.shape = table->potShape();
	deviceDesc.paramTable = d_table;
	cudaParamTableDesc cudaDesc;
	cudaDesc.deviceDesc = deviceDesc;
	cudaDesc.hostResc = deviceDesc;

	return cudaDesc;
}

void as::DeviceDataFactory::disposeDeviceGridUnion(hostGridUnionResource resc, int deviceId)
{
	ASSERT(deviceId >= 0);
	CUDA_CHECK(cudaSetDevice(deviceId));

	hostIntrplGridResource& inner = resc.inner;
	hostIntrplGridResource& outer = resc.outer;
	hostNLGridResource& NL = resc.NL;

	/* Free interpolation grid resources */
	for (uint i = 0; i<inner.numArrays; i++) {
		CUDA_CHECK(cudaFreeArray(inner.cuArrayPtr[i]));
		CUDA_CHECK(cudaDestroyTextureObject(inner.h_texArrayLin[i]));
		CUDA_CHECK(cudaDestroyTextureObject(inner.h_texArrayPt[i]));
	}
	delete[] inner.cuArrayPtr;
	delete[] inner.h_texArrayLin;
	delete[] inner.h_texArrayPt;
	CUDA_CHECK(cudaFree(inner.d_texArrayLin));
	CUDA_CHECK(cudaFree(inner.d_texArrayPt));

	for (uint i = 0; i<outer.numArrays; i++) {
		CUDA_CHECK(cudaFreeArray(outer.cuArrayPtr[i]));
		CUDA_CHECK(cudaDestroyTextureObject(outer.h_texArrayLin[i]));
		CUDA_CHECK(cudaDestroyTextureObject(outer.h_texArrayPt[i]));
	}
	delete[] outer.cuArrayPtr;
	delete[] outer.h_texArrayLin;
	delete[] outer.h_texArrayPt;
	CUDA_CHECK(cudaFree(outer.d_texArrayLin));
	CUDA_CHECK(cudaFree(outer.d_texArrayPt));

	/* Free NL grid resources */
	CUDA_CHECK(cudaFreeArray(NL.cuArray));
	CUDA_CHECK(cudaDestroyTextureObject(NL.tex));
}
void as::DeviceDataFactory::disposeDeviceProtein(hostProteinResource resc, int deviceId)
{
	ASSERT(deviceId >= 0);
	CUDA_CHECK(cudaSetDevice(deviceId));

	CUDA_CHECK(cudaFree(resc.xPos));
	CUDA_CHECK(cudaFree(resc.yPos));
	CUDA_CHECK(cudaFree(resc.zPos));
	CUDA_CHECK(cudaFree(resc.charge));
	CUDA_CHECK(cudaFree(resc.type));
	CUDA_CHECK(cudaFree(resc.mappedType));
	if (resc.numModes != 0) {
		CUDA_CHECK(cudaFree(resc.xModes));
		CUDA_CHECK(cudaFree(resc.yModes));
		CUDA_CHECK(cudaFree(resc.zModes));
	}
}

void as::DeviceDataFactory::disposeDeviceParamTable(hostParamTableResource resc, int deviceId) {
	ASSERT(deviceId >= 0);
	CUDA_CHECK(cudaSetDevice(deviceId));

	CUDA_CHECK(cudaFree(resc.paramTable));

}

/****************************
 * protected member functions
 ****************************/

/****************************
 * private member functions
 ****************************/
as::cudaIntrplGridDesc as::DeviceDataFactory::initIntrpl(const IntrplGrid* grid)
{
	/** array of CUDA texture objects for built-in interpolation */
	cudaTextureObject_t* h_texArrayLin; /** located on host */
	cudaTextureObject_t* d_texArrayLin; /** located on device */
	/** array of CUDA texture objects for manual interpolation */
	cudaTextureObject_t* h_texArrayPt; 	/** located on host */
	cudaTextureObject_t* d_texArrayPt; 	/** located on device */
	cudaArray** h_cuArrayPtr;

	/** The following pointers are deleted when the grid gets detached from the device.
	 * They are needed to destroy device recourses (textures) and are kept in the hostResc object
	 * below.*/
	h_cuArrayPtr = new cudaArray*[grid->numTypes()];
	h_texArrayLin = new cudaTextureObject_t[grid->numTypes()];
	h_texArrayPt = new cudaTextureObject_t[grid->numTypes()];

	// Allocate CUDA array in device memory
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float4>();
	struct cudaExtent cuExtent = make_cudaExtent(grid->width(), grid->height(), grid->depth());
	// For cudaMalloc3DArray: width range in elements.
	// For cudaMalloc3D(...): width range in bytes.

	// Specify texture object parameters
	cudaTextureDesc texDescLin; // for built in device interpolation
	memset(&texDescLin, 0, sizeof(cudaTextureDesc));
	texDescLin.addressMode[0] = cudaAddressModeBorder; // return 0 if out of bounds
	texDescLin.addressMode[1] = cudaAddressModeBorder;
	texDescLin.addressMode[2] = cudaAddressModeBorder;
	texDescLin.filterMode = cudaFilterModeLinear;
	texDescLin.readMode = cudaReadModeElementType;
	texDescLin.normalizedCoords = false;

	cudaTextureDesc texDescPt; // for manual interpolation kernel
	memset(&texDescPt, 0, sizeof(cudaTextureDesc));
	texDescPt.addressMode[0] = cudaAddressModeBorder; // return 0 if out of bounds
	texDescPt.addressMode[1] = cudaAddressModeBorder;
	texDescPt.addressMode[2] = cudaAddressModeBorder;
	texDescPt.filterMode = cudaFilterModePoint;
	texDescPt.readMode = cudaReadModeElementType;
	texDescPt.normalizedCoords = false;

	for (unsigned i = 0; i<grid->numTypes(); i++) {
		cudaArray* &cuArray = h_cuArrayPtr[i];
		CUDA_CHECK(cudaMalloc3DArray(&cuArray, &channelDesc, cuExtent, cudaChannelFormatKindFloat));

		// copy data to 3D array
		cudaMemcpy3DParms copyParams;
		memset(&copyParams, 0, sizeof(copyParams));
		void* gridPtr = (void*)grid->getHostGridPtr(i);
		copyParams.srcPtr   = make_cudaPitchedPtr(gridPtr, cuExtent.width*sizeof(float4), cuExtent.width, cuExtent.height);
		copyParams.dstArray = cuArray;
		copyParams.extent   = cuExtent;
		copyParams.kind     = cudaMemcpyHostToDevice;
		CUDA_CHECK(cudaMemcpy3D(&copyParams));

		// Specify resource
		cudaResourceDesc resDesc;
		memset(&resDesc, 0, sizeof(resDesc));
		resDesc.resType = cudaResourceTypeArray;
		resDesc.res.array.array = cuArray;

		// Create texture objects
		cudaTextureObject_t &texObjLin = h_texArrayLin[i];
		texObjLin = (long long)NULL;
		CUDA_CHECK(cudaCreateTextureObject(&texObjLin, &resDesc, &texDescLin, NULL));

		cudaTextureObject_t &texObjPt = h_texArrayPt[i];
		texObjPt = (long long)NULL;
		CUDA_CHECK(cudaCreateTextureObject(&texObjPt, &resDesc, &texDescPt, NULL));
	}

	CUDA_CHECK(cudaMalloc((void**)&d_texArrayLin, grid->numTypes()*sizeof(cudaTextureObject_t)));
	CUDA_CHECK(cudaMemcpy(d_texArrayLin, h_texArrayLin,grid->numTypes()*sizeof(cudaTextureObject_t), cudaMemcpyHostToDevice));

	CUDA_CHECK(cudaMalloc((void**)&d_texArrayPt, grid->numTypes()*sizeof(cudaTextureObject_t)));
	CUDA_CHECK(cudaMemcpy(d_texArrayPt, h_texArrayPt,grid->numTypes()*sizeof(cudaTextureObject_t), cudaMemcpyHostToDevice));

	/* create deviceIntrplGridDesc */
	deviceIntrplGridDesc deviceDesc;
	deviceDesc.width = grid->width();
	deviceDesc.height = grid->height();
	deviceDesc.depth = grid->depth();
	deviceDesc.dVox = grid->dVox();
	deviceDesc.dVox_inv = grid->dVox_inv();
	deviceDesc.voxelVol_inv = grid->voxelVol_inv();
 	deviceDesc.minDim = grid->pos();
 	deviceDesc.maxDim = grid->maxDim();
 	deviceDesc.texArrayLin = d_texArrayLin;
 	deviceDesc.texArrayPt = d_texArrayPt;

 	/* create hostResc (for device resource deletion) */
 	hostIntrplGridResource hostResc;
 	hostResc.numArrays = grid->numTypes();
 	hostResc.cuArrayPtr = h_cuArrayPtr;
 	hostResc.h_texArrayLin = h_texArrayLin;
    hostResc.h_texArrayPt = h_texArrayPt;
    hostResc.d_texArrayLin = d_texArrayLin;
    hostResc.d_texArrayPt = d_texArrayPt;

    /* create cudaIntrplGridDesc */
    cudaIntrplGridDesc cudaDesc;
    cudaDesc.deviceDesc = deviceDesc;
    cudaDesc.hostResc = hostResc;

    return cudaDesc;
}

as::cudaNLGridDesc as::DeviceDataFactory::initNL(const NLGrid* grid) {

	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 32, 0, 0, cudaChannelFormatKindUnsigned);

	struct cudaExtent cuExtent = make_cudaExtent(grid->width(), grid->height(), grid->depth());
	// For cudaMalloc3DArray: width range in elements.
	// For cudaMalloc3D(...): width range in bytes.

	// Specify texture object parameters
	cudaTextureDesc texDesc;
	memset(&texDesc, 0, sizeof(texDesc));
	texDesc.addressMode[0] = cudaAddressModeBorder; // return 0 if out of bounds
	texDesc.addressMode[1] = cudaAddressModeBorder;
	texDesc.addressMode[2] = cudaAddressModeBorder;
	texDesc.filterMode = cudaFilterModePoint;
	texDesc.readMode = cudaReadModeElementType;
	texDesc.normalizedCoords = false;

	cudaArray* cuArray;
	CUDA_CHECK(cudaMalloc3DArray(&cuArray, &channelDesc, cuExtent));

	// copy data to 3D array
	cudaMemcpy3DParms copyParams;
	memset(&copyParams, 0, sizeof(copyParams));
	copyParams.srcPtr = make_cudaPitchedPtr((void*) grid->grid(), cuExtent.width*sizeof(uint2), cuExtent.width, cuExtent.height);
	copyParams.dstArray = cuArray;
	copyParams.extent   = cuExtent;
	copyParams.kind     = cudaMemcpyHostToDevice;
	CUDA_CHECK(cudaMemcpy3D(&copyParams));

	// Specify resource
	cudaResourceDesc resDesc;
	memset(&resDesc, 0, sizeof(resDesc));
	resDesc.resType = cudaResourceTypeArray;
	resDesc.res.array.array = cuArray;

	// Create texture object
	cudaTextureObject_t texObj;
	texObj = (long long)NULL;
	CUDA_CHECK(cudaCreateTextureObject(&texObj, &resDesc, &texDesc, NULL));

	// Create device neighbor list
	unsigned* d_neighborList;
	CUDA_CHECK(cudaMalloc((void**)&d_neighborList, grid->neighborListSize()*sizeof(unsigned)));
	CUDA_CHECK(cudaMemcpy(d_neighborList, grid->neighborList(),grid->neighborListSize()*sizeof(unsigned), cudaMemcpyHostToDevice));

	/* create deviceNLGridDesc */
	deviceNLGridDesc deviceDesc;
	deviceDesc.width 	= grid->width();
	deviceDesc.height = grid->height();
	deviceDesc.depth	= grid->depth();
	deviceDesc.dVox_inv  = grid->dVox_inv();
	deviceDesc.dPlateau2 = grid->dPlateau2();
	deviceDesc.minDim	= grid->minDim();
	deviceDesc.maxDim	= grid->maxDim();
	deviceDesc.tex = texObj;
	deviceDesc.neighborList = d_neighborList;

	hostNLGridResource hostResc;
	hostResc.tex = texObj;
	hostResc.cuArray = cuArray;

	cudaNLGridDesc cudaDesc;
	cudaDesc.deviceDesc = deviceDesc;
    cudaDesc.hostResc = hostResc;

	return cudaDesc;

}




