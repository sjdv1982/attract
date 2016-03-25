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

#ifndef NEIGHBORLIST_H_
#define NEIGHBORLIST_H_

#include "as/NLGrid.h"
#include "as/Protein.h"
#include "as/ParamTable.h"
#include "as/SimParam.h"

namespace asCore {


template<bool NLOnly>
__global__ void d_NLPotForce(const unsigned gridId,
		const unsigned RecId, const unsigned LigId,
		const unsigned numDOFs,
		const float* LigPosX,
		const float* LigPosY,
		const float* LigPosZ,
		float* outLig_fx,
		float* outLig_fy,
		float* outLig_fz,
		float* outLigand_eEl,
		float* outLigand_eVdW);

void h_NLPotForce(const as::NLGrid *grid,
		const as::Protein* rec, const as::Protein* lig,
		const as::SimParam* simParam,
		const as::AttrParamTable *table,
		const float* LigPosX,
		const float* LigPosY,
		const float* LigPosZ,
		const float* RecPosX,
		const float* RecPosY,
		const float* RecPosZ,
		float* outLig_fx,
		float* outLig_fy,
		float* outLig_fz,
		float* outLig_eEl,
		float* outLig_eVdW,
		float* outRec_fx,
		float* outRec_fy,
		float* outRec_fz);

/*
 ** @brief: ToDo: implementation
 */
template<bool NLOnly>
__global__ void d_NLPotForce_modes(const unsigned gridId,
		const unsigned RecId, const unsigned LigId,
		const unsigned tableId,
		const float* LigPosX,
		const float* LigPosY,
		const float* LigPosZ,
		const float* RecPosX,
		const float* RecPosY,
		const float* RecPosZ,
		float* outLigandX,
		float* outLigandY,
		float* outLigandZ,
		float* outLigand_eEl,
		float* outLigand_eVdW);

}  // namespace asCore




#endif /* NEIGHBORLIST_H_ */
