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

#ifndef INTERPOLATOR_H_
#define INTERPOLATOR_H_

#include "config.h"
#include "asCore/interpolation.h"
#include "asCore/neighborList.h"
#include "as/GridUnion.h"
#include "as/Protein.h"
#include "as/SimParam.h"
#include "as/Comp1_HD.h"
#include "as/Comp3_HD.h"
#include "as/Comp5_HD.h"

namespace asCore {

class Interpolator {
public:
	/* Constructor */
	Interpolator() : _gridSize(0) {};

	/* Destructor */
	~Interpolator() {};

	/***************
	* G E T T E R
	***************/

	/***************
	* S E T T E R
	***************/

	/****************************
	 * public member functions
	 ****************************/
	/*
	 ** @brief: shift > 0 if multiple dofs are processed per stage.
	 */
	inline static void h_PotForce(const as::GridUnion* grid,
		const as::Protein* prot,
		const as::Comp3_HD<float, as::HOSTONLY>* posTr,
		as::Comp5_HD<float, as::HOSTONLY>* outPotForce,
		const unsigned& shift = 0)
	{
		asCore::h_PotForce(grid->innerGrid(), grid->outerGrid(), prot,
				posTr->h_x() + shift,
				posTr->h_y() + shift,
				posTr->h_z() + shift,
				outPotForce->h_x() + shift,
				outPotForce->h_y() + shift,
				outPotForce->h_z() + shift,
				outPotForce->h_v() + shift,
				outPotForce->h_w() + shift);
	}

	inline static void h_NLPotForce(const as::GridUnion* grid,
		const as::Protein* rec, const as::Protein* lig, const as::SimParam* simParam,
		const as::AttrParamTable *table,
		const as::Comp3_HD<float, as::HOSTONLY>* LigPosTr,
		as::Comp5_HD<float, as::HOSTONLY>* outLigPotForce,
		const unsigned& shift = 0)
	{
		asCore::h_NLPotForce(grid->NLgrid(), rec, lig, simParam, table,
				LigPosTr->h_x() + shift,
				LigPosTr->h_y() + shift,
				LigPosTr->h_z() + shift,
				NULL,
				NULL,
				NULL,
				outLigPotForce->h_x() + shift,
				outLigPotForce->h_y() + shift,
				outLigPotForce->h_z() + shift,
				outLigPotForce->h_v() + shift,
				outLigPotForce->h_w() + shift,
				NULL,
				NULL,
				NULL);
	}

	inline static void h_NLPotForce_modes(const as::GridUnion* grid,
		const as::Protein* rec, const as::Protein* lig, const as::SimParam* simParam,
		const as::AttrParamTable *table,
		const as::Comp3_HD<float, as::HOSTONLY>* LigPosTr,
		const as::Comp3_HD<float, as::HOSTONLY>* RecPosTr,
		as::Comp5_HD<float, as::HOSTONLY>* outLigPotForce,
		as::Comp3_HD<float, as::HOSTONLY>* outRecPotForce,
		const unsigned& shift = 0)
	{
		asCore::h_NLPotForce(grid->NLgrid(), rec, lig, simParam, table,
				LigPosTr->h_x() + shift,
				LigPosTr->h_y() + shift,
				LigPosTr->h_z() + shift,
				RecPosTr->h_x() + shift,
				RecPosTr->h_y() + shift,
				RecPosTr->h_z() + shift,
				outLigPotForce->h_x() + shift,
				outLigPotForce->h_y() + shift,
				outLigPotForce->h_z() + shift,
				outLigPotForce->h_v() + shift,
				outLigPotForce->h_w() + shift,
				outRecPotForce->h_x() + shift,
				outRecPotForce->h_y() + shift,
				outRecPotForce->h_z() + shift);
	}

	template<asCore::IntrplType T>
	void d_PotForce (const unsigned& gridId,
		const unsigned& protId,
		const unsigned& numDOFs,
		const as::Comp3_HD<float, as::DEVONLY>* posTr,
		as::Comp5_HD<float, as::DEVONLY>* outPotForce,
		const cudaStream_t &stream = 0);

	// Receptor Gradients not (yet) supported
	template<bool NLOnly>
	void d_NLPotForce (const unsigned& gridId,
		const unsigned& recId, const unsigned& ligId,
		const unsigned& numDOFs,
		const as::Comp3_HD<float, as::DEVONLY>* LigPosTr,
		as::Comp5_HD<float, as::DEVONLY>* outLigPotForce,
		const cudaStream_t &stream = 0);


	/*
	 ** @brief: ToDo: Implementation
	 */
	// Receptor Gradients not (yet) supported
	template<bool NLOnly>
	void d_NLPotForce_modes (const unsigned& NLGridId,
		const unsigned& recId, const unsigned& ligId,
		const unsigned& numDOFs,
		const as::SimParam& simParam,
		const as::Comp3_HD<float, as::DEVONLY>* LigPosTr,
		const as::Comp3_HD<float, as::DEVONLY>* RecPosTr,
		as::Comp5_HD<float, as::DEVONLY>* outLigPotForce,
		as::Comp3_HD<float, as::DEVONLY>* outRecPotForce,
		const cudaStream_t &stream = 0);

	inline void calcAndSetGridSize(const unsigned& numEl) {
		_gridSize = ( numEl + BLSZ_INTRPL -1) / BLSZ_INTRPL;
	}

	/****************************
	 * public member variables
	 ****************************/

protected:
	/****************************
	 * protected member functions
	 ****************************/

	/****************************
	 * protected member variables
	 ****************************/

private:
	/****************************
	 * private member functions
	 ****************************/

	/****************************
	 * private member variables
	 ****************************/

	static const unsigned _blockSize = BLSZ_INTRPL;
	unsigned _gridSize;
};

} // namespace

#endif /* INTERPOLATOR_H_ */
