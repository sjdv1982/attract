
#include "asCore/Transformer.h"
#include "asCore/reduction.h"
#include "asUtils/macros.h"
#include "asClient/DOFTransform.h"
/* Constructor */

/* Destructor */

/****************************
 * public member functions
 ****************************/
using namespace as;

void asCore::Transformer::d_DOF2Pos(const unsigned& protId,
			const unsigned &numDOFs, const Comp1_HD<DOF, DEVONLY>* dofs,
			Comp3_HD<float, DEVONLY>* posTr,
			const cudaStream_t &stream)
{
	cudaVerifyKernel((
			asCore::d_DOF2Pos<<<_gridSize, _blockSize, 0, stream>>>(protId,
				numDOFs, dofs->d_data(),
				posTr->d_x(), posTr->d_y(), posTr->d_z())
			));
}

void asCore::Transformer::d_DOF2Pos_modes(const unsigned& protId,
			const unsigned &numDOFs, const Comp1_HD<DOF, DEVONLY>* dofs,
			Comp3_HD<float, DEVONLY>* posTr,
			Comp3_HD<float, DEVONLY>* posDef,
			const cudaStream_t &stream)
{
	cudaVerifyKernel((
			asCore::d_DOF2Pos_modes<<<_gridSize, _blockSize, 0, stream>>>(protId,
				numDOFs, dofs->d_data(),
				posTr->d_x(), posTr->d_y(), posTr->d_z(),
				posDef->d_x(), posDef->d_y(), posDef->d_z())
	));
}

//void asCore::Transformer::d_partForce2Grad(const unsigned& protId,
//			const unsigned &numDOFs,
//			const Comp5_HD<float, DEVONLY>* outPotForce,
//			Comp1_HD<float, HOST_PINNED>* reduce_res,
//			const cudaStream_t &stream)
//{
//	cudaVerifyKernel((
//			asCore::reduce<float><<<numDOFs, 32, 0, stream>>>(protId,
//				   outPotForce->d_x(),
//				   outPotForce->d_y(),
//				   outPotForce->d_z(),
//				   outPotForce->d_w(),
//				   outPotForce->d_v(),
//				   reduce_res->d_data())
//	));
//
//}

void asCore::Transformer::d_partForce2GradAll(const unsigned& protId,
			const unsigned &numDOFs,
			const unsigned &sizeLigand,
			const Comp5_HD<float, DEVONLY>* outPotForce,
			Comp1_HD<float, HOST_PINNED>* reduce_res,
			const cudaStream_t &stream)
{
	unsigned threads  = (sizeLigand < _blockSizeRed*2) ? asUtils::nextPow2((sizeLigand + 1)/ 2) : _blockSizeRed;
	unsigned blocks = numDOFs;
	asCore::reduceAll(threads, blocks, protId, sizeLigand,
			outPotForce->d_x(),
		    outPotForce->d_y(),
		    outPotForce->d_z(),
		    outPotForce->d_w(),
		    outPotForce->d_v(),
		    reduce_res->d_data(), stream);
}

void asCore::Transformer::h_finalForce2Grad(
			const Protein* prot, // unused
			const unsigned& numDOFs,
			const DOF* dofs,
			Comp1_HD<float, HOST_PINNED>* reduce_res,
			EnGrad* enGrads)
{
	float* redResPtr = reduce_res->h_data();
	for (unsigned i = 0; i < numDOFs; ++i)
	{
		EnGrad &enGrad = enGrads[i];
		enGrad.pos.x = redResPtr[i*14 + 0];
		enGrad.pos.y = redResPtr[i*14 + 1];
		enGrad.pos.z = redResPtr[i*14 + 2];

		for(unsigned j = 0; j < 3; ++j) {
			float magn2 = enGrad.pos.x*enGrad.pos.x + enGrad.pos.y*enGrad.pos.y + enGrad.pos.z*enGrad.pos.z;
			if(magn2 > asCore::ForceLim) {
				enGrad.pos.x *= 0.01;
				enGrad.pos.y *= 0.01;
				enGrad.pos.z *= 0.01;
			}
		}

		enGrad.E_VdW = redResPtr[i*14 + 3];
		enGrad.E_El  = redResPtr[i*14 + 4];

		float torque[3][3] = {0};
		torque[0][0] = redResPtr[i*14 + 5 ];
		torque[0][1] = redResPtr[i*14 + 6 ];
		torque[0][2] = redResPtr[i*14 + 7 ];
		torque[1][0] = redResPtr[i*14 + 8 ];
		torque[1][1] = redResPtr[i*14 + 9 ];
		torque[1][2] = redResPtr[i*14 + 10];
		torque[2][0] = redResPtr[i*14 + 11];
		torque[2][1] = redResPtr[i*14 + 12];
		torque[2][2] = redResPtr[i*14 + 13];

		const DOF &dof = dofs[i];
		float torqueMat[3][3][3];
		asCore::euler2torquemat(dof.ang.x, dof.ang.y, dof.ang.z, torqueMat);

//		using namespace asUtils;
//
//		RotMatf mat_rec;
//		asClient::euler2rotmat(0.0f, 1.0f, 0.0f, mat_rec);
//
//		RotMatf mat_lig;
//		asClient::euler2rotmat(dof.ang.x, dof.ang.y, dof.ang.z, mat_lig);
//
//		RotMatf mat1(0.0f); mat1[0] = 1.0f; mat1[4] = 1.0f; mat1[8] = 1.0f;
//		RotMatf mat2(0.0f); mat2[0] = 1.0f; mat2[4] = 1.0f; mat2[8] = 1.0f;


//		float torqueMat0[3][3][3];
//		std::copy(&torqueMat[0][0][0], &torqueMat[0][0][0] + 27, &torqueMat0[0][0][0]);
//
//
//		for (int j = 0; j < 3; ++j) {
//			for (int l = 0; l < 3; ++l) {
//				torqueMat[0][j][l] = mat1[0]*torqueMat0[0][j][l]+mat1[1]*torqueMat0[1][j][l] + mat1[2]*torqueMat0[2][j][l];
//				torqueMat[1][j][l] = mat1[3]*torqueMat0[0][j][l]+mat1[4]*torqueMat0[1][j][l] + mat1[5]*torqueMat0[2][j][l];
//				torqueMat[2][j][l] = mat1[6]*torqueMat0[0][j][l]+mat1[7]*torqueMat0[1][j][l] + mat1[8]*torqueMat0[2][j][l];
//			}
//		}
//


		/* taken from the original ATTRACT code in rota.f */
		float torPhi = 0;
		float torSsi = 0;
		float torRot = 0;
		for (unsigned k = 0; k < 3; ++k) {
			for (unsigned l = 0; l < 3; ++l) {
				torPhi += torqueMat[k][0][l] * torque[k][l];
				torSsi += torqueMat[k][1][l] * torque[k][l];
				torRot += torqueMat[k][2][l] * torque[k][l];
			}
		}

		enGrad.ang.x = torPhi;
		enGrad.ang.y = torSsi;
		enGrad.ang.z = torRot;
	}
}

/****************************
 * protected member functions
 ****************************/

/****************************
 * private member functions
 ****************************/


