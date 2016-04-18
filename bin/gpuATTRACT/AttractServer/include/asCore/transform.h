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

#ifndef COORDTRANSFORM_H_
#define COORDTRANSFORM_H_

#include <cassert>

#include "asUtils/RotMat.h"
#include "as/asTypes.h"


/*
 ** @brief: This file provides basic transformation routines
 */

namespace asCore {


/*
 ** @brief: Used to transform the ligand.
 */
__global__ void d_DOF2Pos(const unsigned protId,
		const unsigned numDOFs, const as::DOF* dofs,
		float* xTr, float* yTr, float* zTr);


__global__ void d_DOF2Pos_modes(const unsigned protId,
		const unsigned numDOFs, const as::DOF* dofs,
		float* xTr, float* yTr, float* zTr,
		float* xDef, float* yDef, float* zDef);

/*
 ** @brief: Used to transform the receptor.
 */
__global__ void d_DOF2Deform(const unsigned protId,
	const as::DOF* dofs, const unsigned numDOFs,
	float* xDef, float* yDef, float* zDef);


inline void  euler2rotmat(const float& phi, const float& ssi, const float& rot, asUtils::RotMatf& rotmat) {
	//first rotate using rot
	//		this is a rotation around the Z axis
	//			rotating Y into X and X into -Y
	//then rotate using ssi
	//		this is a rotation around the (new) Y axis
	//			rotating X into -Z and Z into X
	//finally, rotate using phi
	//		this is a rotation around the (new) Z axis
	//			rotating X into Y and Y into -X

	float cSSI = cos(ssi);
	float cPHI = cos(phi);
	float sSSI = sin(ssi);
	float sPHI = sin(phi);

	float cSSI_cPHI = cSSI * cPHI;
	float cSSI_sPHI = cSSI * sPHI;
	float sSSI_cPHI = sSSI * cPHI;
	float sSSI_sPHI = sSSI * sPHI;
	float cROT = cos(rot);
	float sROT = sin(rot);

	rotmat.mat[0] = cROT * cSSI_cPHI + sROT * sPHI;
	rotmat.mat[1] = sROT * cSSI_cPHI - cROT * sPHI;
	rotmat.mat[2] = sSSI_cPHI;

	rotmat.mat[3] = cROT * cSSI_sPHI - sROT * cPHI;
	rotmat.mat[4] = sROT * cSSI_sPHI + cROT * cPHI;
	rotmat.mat[5] = sSSI_sPHI;

	rotmat.mat[6] = -cROT * sSSI;
	rotmat.mat[7] = -sROT * sSSI;
	rotmat.mat[8] = cSSI;
}

/*
 ** @brief: rotates the array components in-place = overwriting existing coordinates
 ** @remark: needs verification!!!
 */
inline void rotate(const unsigned& numAtoms, float* x, float* y, float* z, const asUtils::RotMatf& rotmat) {
	for (unsigned i = 0; i < numAtoms; ++i) {
		float x_tmp = x[i];
		float y_tmp = y[i];
		float z_tmp = z[i];
		x[i] = rotmat.mat[0] * x_tmp + rotmat.mat[1] * y_tmp + rotmat.mat[2] * z_tmp;
		y[i] = rotmat.mat[3] * x_tmp + rotmat.mat[4] * y_tmp + rotmat.mat[5] * z_tmp;
		z[i] = rotmat.mat[6] * x_tmp + rotmat.mat[7] * y_tmp + rotmat.mat[8] * z_tmp;
	}
}

inline void h_DOF2Pos(const float* x, const float* y, const float* z,
		const float &xPos, const float &yPos, const float &zPos,
		const float* defVecX, const float* defVecY, const float* defVecZ,
		const asUtils::RotMatf& rotmat, const float* scaleFacMode,
		const unsigned &numAtoms, const unsigned &numModes,
		float* xTr, float* yTr, float* zTr,
		float* xTrDef, float* yTrDef, float* zTrDef)
{
	for (unsigned i = 0; i < numAtoms; ++i) {
		float pos[3];

		pos[0] = x[i];
		pos[1] = y[i];
		pos[2] = z[i];

		for (unsigned mode = 0; mode < numModes; ++mode) {
			float scale = scaleFacMode[mode];
			unsigned idx = i* numModes + mode;
			pos[0] += scale * defVecX[idx];
			pos[1] += scale * defVecY[idx];
			pos[2] += scale * defVecZ[idx];

		}
		if (numModes > 0) {
			xTrDef[i] = pos[0];
			yTrDef[i] = pos[1];
			zTrDef[i] = pos[2];
		}

		float x_tmp = pos[0];
		float y_tmp = pos[1];
		float z_tmp = pos[2];
		pos[0] = rotmat.mat[0] * x_tmp + rotmat.mat[1] * y_tmp + rotmat.mat[2] * z_tmp;
		pos[1] = rotmat.mat[3] * x_tmp + rotmat.mat[4] * y_tmp + rotmat.mat[5] * z_tmp;
		pos[2] = rotmat.mat[6] * x_tmp + rotmat.mat[7] * y_tmp + rotmat.mat[8] * z_tmp;

		pos[0] += xPos;
		pos[1] += yPos;
		pos[2] += zPos;

		xTr[i] = pos[0];
		yTr[i] = pos[1];
		zTr[i] = pos[2];
	}
}

/*
 ** @brief: Used to deform the receptor
 */
inline void h_DOF2Deform(const float* x, const float* y, const float* z,
		const float* defVecX, const float* defVecY, const float* defVecZ,
		const float* scaleFacMode,
		const unsigned &numAtoms, const unsigned &numModes,
		float* xTrDef, float* yTrDef, float* zTrDef)
{
	assert(numModes > 0);
	for (unsigned i = 0; i < numAtoms; ++i) {
		float pos[3];
		pos[0] = x[i];
		pos[1] = y[i];
		pos[2] = z[i];

		for (unsigned mode = 0; mode < numModes; ++mode) {
			float scale = scaleFacMode[mode];
			unsigned idx = i* numModes + mode;
			pos[0] += scale * defVecX[idx];
			pos[1] += scale * defVecY[idx];
			pos[2] += scale * defVecZ[idx];

		}
		xTrDef[i] = pos[0];
		yTrDef[i] = pos[1];
		zTrDef[i] = pos[2];
	}
}

}  // namespace asCore

#endif /* COORDTRANSFORM_H_ */
