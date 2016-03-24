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

#ifndef REDUCTION_H_
#define REDUCTION_H_

#include "asUtils/macros.h"
#include "asUtils/helper.h"


namespace asCore {
/*
 ** @remark: in order to support normal modes it is not enough to supply the protId.
 ** With the protId you can gather only undeformed coordinates. However, for reduction, deformed
 ** coordinates are required.
 ** To integrate normal modes into the pipeline means that one needs to separate deformation
 ** and rotation+translation in terms of execution and memory storage. Otherwise deformation
 ** needs to be calculated twice.
 */
//template <class T>
//__global__ void reduce(unsigned protId,
//		T *d_fx, T *d_fy, T *d_fz, T *d_eVdW, T *d_eEl,
//		T *d_odata);

template <class T>
__host__ void reduceAll(const unsigned& threads, const unsigned& blocks,
		const unsigned& protId, const unsigned& size,
		T *d_fx, T *d_fy, T *d_fz, T *d_eVdW, T *d_eEl,
		T *d_odata,
		const cudaStream_t& stream);

/*
 ** @brief: ToDo: Implementation
 */
template <class T>
__global__ void reduce_modes(unsigned protId,
		float* xDef, float* yDef, float* zDef,
		T *d_fx, T *d_fy, T *d_fz, T *d_eVdW, T *d_eEl,
		T *d_odata);


inline __host__ void euler2torquemat(const float& phi, const float& ssi, const float& rot,
		float torquemat[3][3][3])
{
	float cSSI = cos(ssi);
	float cPHI = cos(phi);
	float sSSI = sin(ssi);
	float sPHI = sin(phi);

	float cSSI_cPHI = cSSI  *  cPHI;
	float cSSI_sPHI = cSSI  *  sPHI;
	float sSSI_cPHI = sSSI  *  cPHI;
	float sSSI_sPHI = sSSI  *  sPHI;
	float cROT = cos(rot);
	float sROT = sin(rot);

	torquemat[0][0][0] = -cROT * cSSI_sPHI + sROT * cPHI;
	torquemat[0][0][1] = -sROT * cSSI_sPHI - cROT * cPHI; // -rotmat[4]
	torquemat[0][0][2] = -sSSI_sPHI; // -rotmat[5]
	torquemat[1][0][0] = cROT * cSSI_cPHI + sROT * sPHI; // rotmat[0]
	torquemat[1][0][1] = sROT * cSSI_cPHI - cROT * sPHI; // rotmat[1]
	torquemat[1][0][2] = sSSI_cPHI; // rotmat[2]
	torquemat[2][0][0] = 0.0;
	torquemat[2][0][1] = 0.0;
	torquemat[2][0][2] = 0.0;

	torquemat[0][1][0] = -cROT * sSSI_cPHI;
	torquemat[0][1][1] = -sROT * sSSI_cPHI;
	torquemat[0][1][2] = cSSI_cPHI;
	torquemat[1][1][0] = -cROT * sSSI_sPHI;
	torquemat[1][1][1] = -sROT * sSSI_sPHI;
	torquemat[1][1][2] = cSSI_sPHI;
	torquemat[2][1][0] = -cROT * cSSI;
	torquemat[2][1][1] = -sROT * cSSI; // rotmat[7]
	torquemat[2][1][2] = -sSSI;

	torquemat[0][2][0] = -sROT * cSSI_cPHI + cROT * sPHI;
	torquemat[0][2][1] = cROT * cSSI_cPHI + sROT * sPHI;
	torquemat[0][2][2] = 0.0;
	torquemat[1][2][0] = -sROT * cSSI_sPHI - cROT * cPHI;
	torquemat[1][2][1] = cROT * cSSI_sPHI - sROT * cPHI; // rotmat[3]
	torquemat[1][2][2] = 0.0;
	torquemat[2][2][0] = sROT * sSSI;
	torquemat[2][2][1] = -cROT * sSSI; // rotmat[6]
	torquemat[2][2][2] = 0.0;
}


/*
 ** @brief: reduction/accumulation of energies and forces to calc. the translational gradients
 */
const float ForceLim = 1.0e18;
inline void redPotForce(const float* fx, const float* fy, const float* fz,
		const float* energies0, const float* energies1,
		const unsigned &numAtoms,
		float &accFx, float &accFy, float &accFz, float &energy0, float &energy1)
{
	accFx = 0;
	accFy = 0;
	accFz = 0;
	energy0 = 0;
	energy1 = 0;
	for(unsigned i = 0; i < numAtoms; ++i) {
		accFx += fx[i];
		accFy += fy[i];
		accFz += fz[i];
		energy0 += energies0[i];
		energy1 += energies1[i];
	}

	// force reduction, some times helps in case of very "bad" start structure
	// taken from original ATTRACT code in trans.f
	for(unsigned i = 0; i < 3; ++i) {
		float magn2 = accFx*accFx + accFy*accFy + accFz*accFz;
		if(magn2 > ForceLim) {
			accFx *= 0.01;
			accFy *= 0.01;
			accFz *= 0.01;
		}
	}
}

/*
 ** @brief: reduction/transformation of forces to calc. the rotational gradients
 */
inline void redTorque(const float* x, const float* y, const float* z,
		const float* fx, const float* fy, const float* fz ,
		const unsigned &numAtoms, const float torqueMat[3][3][3],
		float &torPhi, float &torSsi, float &torRot)
{
	torPhi = 0;
	torSsi = 0;
	torRot = 0;

	float torque[3][3] = {0};
	for (unsigned i = 0; i < numAtoms; ++i) {
		torque[0][0] += x[i]*fx[i];
		torque[0][1] += y[i]*fx[i];
		torque[0][2] += z[i]*fx[i];
		torque[1][0] += x[i]*fy[i];
		torque[1][1] += y[i]*fy[i];
		torque[1][2] += z[i]*fy[i];
		torque[2][0] += x[i]*fz[i];
		torque[2][1] += y[i]*fz[i];
		torque[2][2] += z[i]*fz[i];
	}

	for (unsigned k = 0; k < 3; ++k) {
		for (unsigned l = 0; l < 3; ++l) {
			torPhi += torqueMat[k][0][l] * torque[k][l];
			torSsi += torqueMat[k][1][l] * torque[k][l];
			torRot += torqueMat[k][2][l] * torque[k][l];
		}
	}
}

/*
 ** @remark: needs verification!!!
 */
inline void redModes(const float* xModes, const float* yModes, const float* zModes,
		const float* fx, const float* fy, const float* fz ,
		const unsigned &numAtoms, const unsigned &numModes,
		float* modes)
{
	for (unsigned i = 0; i < numAtoms; ++i) {
		unsigned idx = i*numModes;
		for (unsigned modeId = 0; modeId < numModes; ++modeId) {
			unsigned mode = idx + modeId;
			float grad = fx[i]*xModes[mode] + fy[i]*yModes[mode] + fz[i]*zModes[mode];

			/* -=: taken from the original ATTRACT code in ligmin.f*/
			modes[modeId] -= grad;
		}
	}
}

}  // namespace asCore

#endif /* REDUCTION_H_ */
