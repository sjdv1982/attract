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

#ifndef TRANSFORMER_H_
#define TRANSFORMER_H_

#include "config.h"
#include "as/asTypes.h"
#include "as/Protein.h"
#include "asUtils/RotMat.h"
#include "as/asBuffer.h"
#include "asCore/transform.h"
#include "asCore/reduction.h"

namespace asCore {


/*
 ** @remark: in order to support normal modes it is not enough to supply the protId.
 ** With the protId you can gather only undeformed coordinates. However, for reduction, deformed
 ** coordinates are required.
 ** To integrate normal modes into the pipeline means that one needs to separate deformation
 ** and rotation+translation in terms of execution and memory storage. Otherwise deformation
 ** needs to be calculated twice.
 */


class Transformer {

public:

	/* Constructor */
	Transformer() : _gridSize(0) {}

	/* Destructor */
	~Transformer() {}

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
	static void h_DOF2Pos(const as::Protein* prot,
			const as::DOF& dof, const asUtils::RotMatf& rotmat,
			as::Comp3_HD<float, as::HOSTONLY>* posTr,
			const unsigned& shift = 0)
	{
		*(posTr->h_conf() + shift) = dof.conf;
		asCore::h_DOF2Pos(prot->xPos(dof.conf),prot->yPos(dof.conf), prot->zPos(dof.conf),
				dof.pos.x, dof.pos.y, dof.pos.z,
				prot->xModes(),prot->yModes(), prot->zModes(),
				rotmat, dof.modes,
				prot->nAtoms(), prot->numModes(),
				posTr->h_x() + shift, posTr->h_y() + shift, posTr->h_z() + shift,
				NULL, NULL, NULL);
	}

	static void h_DOF2Pos_modes(const as::Protein* prot,
			const as::DOF& dof, const asUtils::RotMatf& rotmat,
			as::Comp3_HD<float, as::HOSTONLY>* posTr,
			as::Comp3_HD<float, as::HOSTONLY>* posDef,
			const unsigned& shift = 0)
	{
		*(posTr->h_conf() + shift) = dof.conf;
		*(posDef->h_conf() + shift) = dof.conf;
		asCore::h_DOF2Pos(prot->xPos(dof.conf),prot->yPos(dof.conf), prot->zPos(dof.conf),
				dof.pos.x, dof.pos.y, dof.pos.z,
				prot->xModes(),prot->yModes(), prot->zModes(),
				rotmat, dof.modes,
				prot->nAtoms(), prot->numModes(),
				posTr->h_x() + shift, posTr->h_y() + shift, posTr->h_z() + shift,
				posDef->h_x() + shift, posDef->h_y() + shift, posDef->h_z() + shift);
	}

	/*
	 ** @brief: Used to transform receptor coordinates.
	 */
	static void h_DOF2Deform(const as::Protein* prot,
			const as::DOF& dof,
			as::Comp3_HD<float, as::HOSTONLY>* posDef,
			const unsigned& shift = 0)
	{
		asCore::h_DOF2Deform(prot->xPos(dof.conf),prot->yPos(dof.conf), prot->zPos(dof.conf),
				prot->xModes(),prot->yModes(), prot->zModes(),
				dof.modes,
				prot->nAtoms(), prot->numModes(),
				posDef->h_x() + shift, posDef->h_x() + shift, posDef->h_x() + shift);

	}


	/*
	 ** @brief: need verification for mode reduction!!!
	 */
	static void h_RedEnForce(const as::Protein* prot,
			const as::DOF& dof, const asUtils::RotMatf& rotmat,
			const as::Comp5_HD<float, as::HOSTONLY>* outPotForce,
			const as::Comp3_HD<float, as::HOSTONLY>* posDef,
			as::EnGrad& gradEn,
			const unsigned& shift = 0)
	{
		const unsigned& nAtoms = prot->nAtoms();
		asCore::redPotForce(
				outPotForce->h_x() + shift,
				outPotForce->h_y() + shift,
				outPotForce->h_z() + shift,
				outPotForce->h_w() + shift,
				outPotForce->h_v() + shift,
				nAtoms,
				gradEn.pos.x, gradEn.pos.y, gradEn.pos.z, gradEn.E_VdW, gradEn.E_El);

		float torqueMat[3][3][3];
		asCore::euler2torquemat(dof.ang.x, dof.ang.y, dof.ang.z, torqueMat);

		const unsigned& numModes = prot->numModes();
		if (numModes <= 0) {
			asCore::redTorque(prot->xPos(dof.conf), prot->yPos(dof.conf), prot->zPos(dof.conf),
				outPotForce->h_x() + shift,
				outPotForce->h_y() + shift,
				outPotForce->h_z() + shift,
				nAtoms, torqueMat,
				gradEn.ang.x, gradEn.ang.y, gradEn.ang.z);
		} else {
			// verification needed! Was never tested!

			asCore::redTorque(posDef->h_x(), posDef->h_y(), posDef->h_z(),
				outPotForce->h_x() + shift,
				outPotForce->h_y() + shift,
				outPotForce->h_z() + shift,
				nAtoms, torqueMat,
				gradEn.ang.x, gradEn.ang.y, gradEn.ang.z);

			/* invert rotmat */
			asUtils::RotMatf invMat = rotmat.getInv();

			/* in-place Rotation! */
			asCore::rotate(nAtoms, posDef->h_x(), posDef->h_y(), posDef->h_z(), invMat);

			asCore::redModes(prot->xModes(), prot->yModes(), prot->zModes(),
				outPotForce->h_x() + shift,
				outPotForce->h_y() + shift,
				outPotForce->h_z() + shift,
				nAtoms, numModes,
				gradEn.modes);
		}
	}

	void d_DOF2Pos(const unsigned& protId,
			const unsigned &numDOFs, const as::Comp1_HD<as::DOF, as::DEVONLY>* dofs,
			as::Comp3_HD<float, as::DEVONLY>* posTr,
			const cudaStream_t &stream = 0);

	void d_DOF2Pos_modes(const unsigned& protId,
			const unsigned &numDOFs, const as::Comp1_HD<as::DOF, as::DEVONLY>* dofs,
			as::Comp3_HD<float, as::DEVONLY>* posTr,
			as::Comp3_HD<float, as::DEVONLY>* posDef,
			const cudaStream_t &stream = 0);


	void d_partForce2Grad(const unsigned& protId,
			const unsigned &numDOFs,
			const as::Comp5_HD<float, as::DEVONLY>* outPotForce,
			as::Comp1_HD<float, as::HOST_PINNED>* reduce_res,
			const cudaStream_t &stream = 0);

	void d_partForce2GradAll(const unsigned& protId,
			const unsigned &numDOFs,
			const unsigned &sizeLigand,
			const as::Comp5_HD<float, as::DEVONLY>* outPotForce,
			as::Comp1_HD<float, as::HOST_PINNED>* reduce_res,
			const cudaStream_t &stream = 0);

	/*
	 ** @brief: ToDo: Implementation
	 */
	void d_partForce2Grad_modes(const unsigned& protId,
			const unsigned& numDOFs,
			const as::Comp1_HD<as::DOF, as::DEVONLY>* dofs,
			const as::Comp5_HD<float, as::DEVONLY>* outPotForce,
			const as::Comp3_HD<float, as::DEVONLY>* posDef,
			as::Comp1_HD<float, as::HOST_PINNED>* reduce_res,
			const cudaStream_t &stream = 0);

	static void h_finalForce2Grad(const as::Protein* prot,
			const unsigned& numDOFs,
			const as::DOF* dofs,
			as::Comp1_HD<float, as::HOST_PINNED>* reduce_res,
			as::EnGrad* enGrads);

	/*
	 ** @brief: ToDo: Implementation
	 */
	static void h_finalForce2Grad_modes(const as::Protein* prot,
			const unsigned& numDOFs,
			const as::Comp1_HD<as::DOF, as::DEVONLY>* dofs,
			as::Comp1_HD<float, as::HOST_PINNED>* reduce_res,
			as::Comp1_HD<as::EnGrad, as::HOSTONLY>* enGrads);

	void calcAndSetGridSize(const unsigned& numEl) {
		_gridSize = ( numEl + BLSZ_TRAFO -1) / BLSZ_TRAFO;
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
	static const unsigned _blockSize = BLSZ_TRAFO;
	static const unsigned _blockSizeRed = BLSZ_REDUCE;
	unsigned _gridSize;

};

} // namespace

#endif /* TRANSFORMER_H_ */
