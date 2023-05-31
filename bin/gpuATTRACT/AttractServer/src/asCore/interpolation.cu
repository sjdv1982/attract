

#include <stdio.h>

#include "asCore/interpolation.h"
#include "asUtils/cudaMath.h"


extern __constant__ as::deviceProteinDesc c_Proteins[DEVICE_MAXPROTEINS];
extern __constant__ as::deviceGridUnionDesc c_Grids[DEVICE_MAXGRIDS];

/*
 ** @brief: returns the 8 voxels for device interpolation given a set of coordinates
 */
__forceinline__ __device__ void getVoxelDevice(const as::deviceIntrplGridDesc& grid, const unsigned &type,
		const float &x, const float &y,	const float &z,  as::VoxelOctet& voxelOct)
{
	unsigned idxX = (unsigned) floor(
			(x - grid.minDim.x) * grid.dVox_inv);
	unsigned idxY = (unsigned) floor(
			(y - grid.minDim.y) * grid.dVox_inv);
	unsigned idxZ = (unsigned) floor(
			(z - grid.minDim.z) * grid.dVox_inv);

	// compute absolute position of the vertices
	voxelOct.min.x = idxX * grid.dVox + grid.minDim.x;
	voxelOct.min.y = idxY * grid.dVox + grid.minDim.y;
	voxelOct.min.z = idxZ * grid.dVox + grid.minDim.z;
	voxelOct.max.x = voxelOct.min.x + grid.dVox;
	voxelOct.max.y = voxelOct.min.y + grid.dVox;
	voxelOct.max.z = voxelOct.min.z + grid.dVox;

	/** use of non-normalized coordinates */

	float idxNx = 0.5 + (float) idxX;
	float idxNy = 0.5 + (float) idxY;
	float idxNz = 0.5 + (float) idxZ;

	voxelOct.data[0][0][0] = tex3D<float4>(grid.texArrayPt[type], idxNx, idxNy, idxNz);
	voxelOct.data[1][0][0] = tex3D<float4>(grid.texArrayPt[type], idxNx + 1, idxNy, idxNz);
	voxelOct.data[0][1][0] = tex3D<float4>(grid.texArrayPt[type], idxNx, idxNy + 1, idxNz);
	voxelOct.data[1][1][0] = tex3D<float4>(grid.texArrayPt[type], idxNx + 1, idxNy + 1, idxNz);

	voxelOct.data[0][0][1] = tex3D<float4>(grid.texArrayPt[type], idxNx, idxNy, idxNz + 1);
	voxelOct.data[1][0][1] = tex3D<float4>(grid.texArrayPt[type], idxNx + 1, idxNy, idxNz + 1);
	voxelOct.data[0][1][1] = tex3D<float4>(grid.texArrayPt[type], idxNx, idxNy + 1, idxNz + 1);
	voxelOct.data[1][1][1] = tex3D<float4>(grid.texArrayPt[type], idxNx + 1, idxNy + 1, idxNz + 1);

}


/*
 ** @brief: function body for a trilinear interpolation.
 */
__forceinline__ __host__ __device__ void trilinearInterpolation(const float &x,
		const float &y, const float &z, const as::VoxelOctet &voxelOct,
		const float &voxelVol_inv, float4 &V)
{
	/* for operator overloading of *,+,-,/ for cuda types (float4)
	 * they are defined in asUtils/cudaMath*/
	using namespace asUtils;

	float3 pos = make_float3(x, y, z);
	float3 pos_m_posMin = pos - voxelOct.min;
	float3 posMax_m_pos = voxelOct.max - pos;

	float tmpMax_xy = (posMax_m_pos.x) * (posMax_m_pos.y);
	float tmpMin_xy = (pos_m_posMin.x) * (pos_m_posMin.y);

	V = voxelOct.data[0][0][0] * (tmpMax_xy * (posMax_m_pos.z))
			+ voxelOct.data[1][0][0]
					* ((pos_m_posMin.x) * (posMax_m_pos.y) * (posMax_m_pos.z))
			+ voxelOct.data[0][1][0]
					* ((posMax_m_pos.x) * (pos_m_posMin.y) * (posMax_m_pos.z))
			+ voxelOct.data[0][0][1] * (tmpMax_xy * (pos_m_posMin.z))
			+ voxelOct.data[1][0][1]
					* ((pos_m_posMin.x) * (posMax_m_pos.y) * (pos_m_posMin.z))
			+ voxelOct.data[0][1][1]
					* ((posMax_m_pos.x) * (pos_m_posMin.y) * (pos_m_posMin.z))
			+ voxelOct.data[1][1][0] * (tmpMin_xy * (posMax_m_pos.z))
			+ voxelOct.data[1][1][1] * (tmpMin_xy * (pos_m_posMin.z));

	V = V * voxelVol_inv;
	return;
}

__forceinline__ __device__ float4 Intrpl3D(const as::deviceIntrplGridDesc& grid, const unsigned& type, const float &x, const float &y,
		const float &z)
{

	as::VoxelOctet voxelOct;
	getVoxelDevice(grid, type, x, y, z, voxelOct);
	float4 V;
	trilinearInterpolation(x, y, z, voxelOct, grid.voxelVol_inv, V);
	return V;
}




template<asCore::IntrplType T>
__global__ void asCore::d_InnerPotForce(
		const unsigned gridId, const unsigned protId,
		const unsigned numDOFs,
		const float* data_in_x, const float* data_in_y, const float* data_in_z,
		float* data_out_x, float* data_out_y, float* data_out_z,
		float* data_out_eEl, float* data_out_eVdW)
{
	/* for operator overloading of *,+,-,/ for cuda types (float4)
	 * they are defined in asUtils/cudaMath*/
	using namespace asUtils;

	/* calculate element index that is to be processed */
	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;
	const unsigned nAtoms = c_Proteins[protId].nAtoms;
	if (idx < nAtoms*numDOFs) {
		unsigned type = c_Proteins[protId].mappedType[idx % nAtoms];
		float4 V_el = {0};
		float4 V_VdW = {0};
		if (type != 0) {

			float x = data_in_x[idx];
			float y = data_in_y[idx];
			float z = data_in_z[idx];

			if ((x >= c_Grids[gridId].inner.minDim.x && x <= c_Grids[gridId].inner.maxDim.x)
					&& (y >= c_Grids[gridId].inner.minDim.y && y <= c_Grids[gridId].inner.maxDim.y)
					&& (z >= c_Grids[gridId].inner.minDim.z && z <= c_Grids[gridId].inner.maxDim.z)) {

				// interpolate
				/* switch is optimized out by the compiler. */
				switch (T) {
				case built_in:
					x = (x - c_Grids[gridId].inner.minDim.x) * c_Grids[gridId].inner.dVox_inv + 0.5f;
					y = (y - c_Grids[gridId].inner.minDim.y) * c_Grids[gridId].inner.dVox_inv + 0.5f;
					z = (z - c_Grids[gridId].inner.minDim.z) * c_Grids[gridId].inner.dVox_inv + 0.5f;

					V_VdW = tex3D<float4>(c_Grids[gridId].inner.texArrayLin[type], x, y, z); /** Interpolated value */
					break;
				case manual:
					V_VdW = Intrpl3D(c_Grids[gridId].inner, type, x, y, z); /** Interpolated value */
				}

				float charge = c_Proteins[protId].charge[idx % nAtoms];
				if (fabs(charge) > 0.001f) {

					switch (T) {
					case built_in:
						V_el = tex3D<float4>(c_Grids[gridId].inner.texArrayLin[0], x, y, z); /** Interpolated value */
						break;
					case manual:
						V_el = Intrpl3D(c_Grids[gridId].inner, 0, x, y, z); /** Interpolated value */
					}

					V_el = V_el * charge;
				}
			}
		}
		data_out_x[idx] = V_VdW.x + V_el.x;
		data_out_y[idx] = V_VdW.y + V_el.y;
		data_out_z[idx] = V_VdW.z + V_el.z;
		data_out_eVdW[idx] = V_VdW.w;
		data_out_eEl[idx] = V_el.w;
	}
}

template<asCore::IntrplType T>
__global__ void asCore::d_OuterPotForce(
		const unsigned gridId, const unsigned protId,
		const unsigned numDOFs,
		const float* data_in_x, const float* data_in_y, const float* data_in_z,
		float* data_out_x, float* data_out_y, float* data_out_z,
		float* data_out_eEl, float* data_out_eVdw)
{
	/* for operator overloading of *,+,-,/ for cuda types (float4)
	 * they are defined in asUtils/cudaMath*/
	using namespace asUtils;

	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;
	const unsigned nAtoms = c_Proteins[protId].nAtoms;
	if (idx < nAtoms*numDOFs) {
		unsigned type = c_Proteins[protId].mappedType[idx % nAtoms];
		if (type != 0) {

			float x = data_in_x[idx];
			float y = data_in_y[idx];
			float z = data_in_z[idx];

			//DEBUG
//			if (idx >703 && idx < 713) {
//				printf("ATOM %u %f %f %f\n", idx, x, y, z);
//
//			}

			//DEBUG
//			if (idx == 708) {
//				printf("ATOM %u %f %f %f\n", idx, x, y, z);
//				printf("%f %f %f %f %f %f\n%f %f %f %f %f %f" ,
//						c_Grids[gridId].inner.minDim.x, c_Grids[gridId].inner.minDim.y, c_Grids[gridId].inner.minDim.z,
//						c_Grids[gridId].inner.maxDim.x, c_Grids[gridId].inner.maxDim.y, c_Grids[gridId].inner.maxDim.z,
//						c_Grids[gridId].outer.minDim.x, c_Grids[gridId].outer.minDim.y, c_Grids[gridId].outer.minDim.z,
//						c_Grids[gridId].outer.maxDim.x, c_Grids[gridId].outer.maxDim.y, c_Grids[gridId].outer.maxDim.z);
//			}

			if (      ((x < c_Grids[gridId].inner.minDim.x || x > c_Grids[gridId].inner.maxDim.x)
					|| (y < c_Grids[gridId].inner.minDim.y || y > c_Grids[gridId].inner.maxDim.y)
					|| (z < c_Grids[gridId].inner.minDim.z || z > c_Grids[gridId].inner.maxDim.z))
					&&
					  ((x >= c_Grids[gridId].outer.minDim.x && x <= c_Grids[gridId].outer.maxDim.x)
					&& (y >= c_Grids[gridId].outer.minDim.y && y <= c_Grids[gridId].outer.maxDim.y)
					&& (z >= c_Grids[gridId].outer.minDim.z && z <= c_Grids[gridId].outer.maxDim.z))) {

				//DEBUG
//				if (idx == 708)
//					printf("ATOM %u outer", idx);

				// interpolate
				/* switch is optimized out by the compiler. */
				float4 V_VdW;

				switch (T) {
				case built_in:
					x = (x - c_Grids[gridId].outer.minDim.x) * c_Grids[gridId].outer.dVox_inv + 0.5f;
					y = (y - c_Grids[gridId].outer.minDim.y) * c_Grids[gridId].outer.dVox_inv + 0.5f;
					z = (z - c_Grids[gridId].outer.minDim.z) * c_Grids[gridId].outer.dVox_inv + 0.5f;

					V_VdW = tex3D<float4>(c_Grids[gridId].outer.texArrayLin[type], x, y, z); /** Interpolated value */
					break;
				case manual:
					V_VdW = Intrpl3D(c_Grids[gridId].outer, type, x, y, z); /** Interpolated value */
				}

				float charge = c_Proteins[protId].charge[idx % nAtoms];
				float4 V_el = {0};
				if (fabs(charge) > 0.001f) {

					switch (T) {
					case built_in:
						V_el = tex3D<float4>(c_Grids[gridId].outer.texArrayLin[0], x, y, z); /** Interpolated value */
						break;
					case manual:
						V_el = Intrpl3D(c_Grids[gridId].outer, 0, x, y, z); /** Interpolated value */
					}

					V_el = V_el * charge;
					data_out_eEl[idx] = V_el.w;
				}

				data_out_x[idx] = V_VdW.x + V_el.x;
				data_out_y[idx] = V_VdW.y + V_el.y;
				data_out_z[idx] = V_VdW.z + V_el.z;
				data_out_eVdw[idx] = V_VdW.w;
			}
		}
	}
}

void asCore::h_PotForce(const as::IntrplGrid* innerGrid,
		const as::IntrplGrid* outerGrid, const as::Protein* prot,
		const float* LigPosX,
		const float* LigPosY,
		const float* LigPosZ,
		float* data_out_x, float* data_out_y, float* data_out_z,
		float* data_out_eEl, float* data_out_eVdw)
{
	/* for operator overloading of *,+,-,/ for cuda types (float4) */
	using namespace asUtils;

	const unsigned nAtoms = prot->nAtoms();
	/* loop over all elements in LigPos/output */
	for (unsigned i = 0; i < nAtoms; ++i) {
		const unsigned type = prot->mappedType()[i];
		if (type == 0)
			continue;
		float3 pos;
		pos.x = LigPosX[i];
		pos.y = LigPosY[i];
		pos.z = LigPosZ[i];
		const float charge = prot->charge()[i];

		//DEBUG

		//DEBUG
//		if (i >703 && i < 713) {
//			printf("ATOM %u %f %f %f\n", i, pos.x, pos.y, pos.z);
//		}

//		if (i == 36) {
////		if (true) {
//			printf("ATOM %u %f %f %f\n", i, pos.x, pos.y, pos.z);
//			printf("%f %f %f %f %f %f\n%f %f %f %f %f %f" ,
//						innerGrid->minDim().x, innerGrid->minDim().y, innerGrid->minDim().z,
//						innerGrid->maxDim().x, innerGrid->maxDim().y, innerGrid->maxDim().z,
//						outerGrid->minDim().x, outerGrid->minDim().y, outerGrid->minDim().z,
//						outerGrid->maxDim().x, outerGrid->maxDim().y, outerGrid->maxDim().z);
//		}


		// fetch voxel from grid
		as::VoxelOctet voxel;
		float4 V_VdW, V_el; //interpolated value for Van-der-Waals & electorstatic interactions
		V_el = V_VdW = make_float4(0.0f, 0.0f, 0.0f, 0.0f);

		if (innerGrid->outOfBounds_byPos(pos.x, pos.y, pos.z)) {

			if (!outerGrid->outOfBounds_byPos(pos.x, pos.y, pos.z)) {

				int idxX, idxY, idxZ;
				outerGrid->getIndex(pos.x, pos.y, pos.z, idxX, idxY, idxZ);

				//DEBUG
//				if (i == 708)
//					std::cout << "outer" << std::endl;

				/* VdW - Forces/Energy */
				outerGrid->host_getVoxelByIndex(idxX, idxY, idxZ, type, voxel);

				trilinearInterpolation(pos.x, pos.y, pos.z, voxel,
						outerGrid->voxelVol_inv(), V_VdW);
				/* El.stat. - Forces/Energy */
				if (fabs(charge) > 0.001f) {
					outerGrid->host_getVoxelByIndex(idxX, idxY, idxZ, 0, voxel);
					trilinearInterpolation(pos.x, pos.y, pos.z, voxel,
							outerGrid->voxelVol_inv(), V_el);
					V_el = V_el * charge;
				}

			}

		} else {
			//DEBUG
//			if (i == 708)
//				std::cout << "inner" << std::endl;

			int idxX, idxY, idxZ;
			innerGrid->getIndex(pos.x, pos.y, pos.z, idxX, idxY, idxZ);

			/* We are using the inner grid and can fetch the voxel by index */
			innerGrid->host_getVoxelByIndex(idxX, idxY, idxZ, type, voxel);

			trilinearInterpolation(pos.x, pos.y, pos.z, voxel,
					innerGrid->voxelVol_inv(), V_VdW);

			/* El.stat. - Forces/Energy */
			if (fabs(charge) > 0.001f) {
				innerGrid->host_getVoxelByIndex(idxX, idxY, idxZ, 0, voxel);
				trilinearInterpolation(pos.x, pos.y, pos.z, voxel,
						innerGrid->voxelVol_inv(), V_el);
				V_el = V_el * charge;
			}

		}

		data_out_x[i] = V_VdW.x + V_el.x;
		data_out_y[i] = V_VdW.y + V_el.y;
		data_out_z[i] = V_VdW.z + V_el.z;
		data_out_eVdw[i] = V_VdW.w;
		data_out_eEl[i] = V_el.w;

	}

	return;
}






/* explicit instantiation */
template __global__ void
asCore::d_InnerPotForce<asCore::built_in> (const unsigned gridId, const unsigned protId,
		const unsigned numDOFs,
		const float* data_in_x, const float* data_in_y, const float* data_in_z,
		float* data_out_x, float* data_out_y, float* data_out_z, float* data_out_eEl,
		float* data_out_eVdW);

template __global__ void
asCore::d_InnerPotForce<asCore::manual> (const unsigned gridId, const unsigned protId,
		const unsigned numDOFs,
		const float* data_in_x, const float* data_in_y, const float* data_in_z,
		float* data_out_x, float* data_out_y, float* data_out_z, float* data_out_eEl,
		float* data_out_eVdW);

template __global__ void
asCore::d_OuterPotForce<asCore::built_in>(
		const unsigned gridId, const unsigned protId,
		const unsigned numDOFs,
		const float* data_in_x, const float* data_in_y, const float* data_in_z,
		float* data_out_x, float* data_out_y, float* data_out_z, float* data_out_eEl,
		float* data_out_eVdW);

template __global__ void
asCore::d_OuterPotForce<asCore::manual>(
		const unsigned gridId, const unsigned protId,
		const unsigned numDOFs,
		const float* data_in_x, const float* data_in_y, const float* data_in_z,
		float* data_out_x, float* data_out_y, float* data_out_z, float* data_out_eEl,
		float* data_out_eVdW);

