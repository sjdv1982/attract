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

#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

#include "cuda_runtime.h"

#include "as/IntrplGrid.h"
#include "as/Protein.h"


namespace asCore {

/*
 ** @brief: Interpolation type.
 ** built_in: 	use cuda built_in interpolation routine
 ** manual: 	use "home-made" routine
 */
enum IntrplType {
	built_in,
	manual
};

template<IntrplType T>
__global__ void d_InnerPotForce(
		const unsigned gridId, const unsigned protId,
		const unsigned numDOFs,
		const float* data_in_x, const float* data_in_y, const float* data_in_z,
		float* data_out_x, float* data_out_y, float* data_out_z,
		float* data_out_eEl, float* data_out_eVdw);

template<IntrplType T>
__global__ void d_OuterPotForce(
		const unsigned gridId, const unsigned protId,
		const unsigned numDOFs,
		const float* data_in_x, const float* data_in_y, const float* data_in_z,
		float* data_out_x, float* data_out_y, float* data_out_z,
		float* data_out_eEl, float* data_out_eVdw);


void h_PotForce(const as::IntrplGrid* innerGrid,
		const as::IntrplGrid* outerGrid, const as::Protein* prot,
		const float* LigPosX,
		const float* LigPosY,
		const float* LigPosZ,
		float* data_out_x, float* data_out_y, float* data_out_z,
		float* data_out_eEl, float* data_out_eVdw);



}  // namespace asCore


#endif /* INTERPOLATION_H_ */
