
#include "asCore/Interpolator.h"
#include "asUtils/macros.h"

using namespace as;

/* Constructor */


/* Destructor */



/****************************
 * public member functions
 ****************************/
template<asCore::IntrplType T>
void asCore::Interpolator::d_PotForce(const unsigned& gridId,
		const unsigned& protId,
		const unsigned& numDOFs,
		const Comp3_HD<float, DEVONLY>* posTr,
		Comp5_HD<float, DEVONLY>* outPotForce,
		const cudaStream_t &stream)
		{
	cudaVerifyKernel((
			asCore::d_InnerPotForce<T><<<_gridSize, _blockSize, 0, stream>>>(
					gridId, protId, numDOFs,
					posTr->d_x(),
					posTr->d_y(),
			        posTr->d_z(),
			        outPotForce->d_x(),
			        outPotForce->d_y(),
			        outPotForce->d_z(),
			        outPotForce->d_v(),
			        outPotForce->d_w())
			));
	cudaVerifyKernel((
			asCore::d_OuterPotForce<T><<<_gridSize, _blockSize, 0, stream>>>(
					gridId, protId, numDOFs,
					posTr->d_x(),
					posTr->d_y(),
			        posTr->d_z(),
			        outPotForce->d_x(),
			        outPotForce->d_y(),
			        outPotForce->d_z(),
			        outPotForce->d_v(),
			        outPotForce->d_w())
			));
}

template<bool NLOnly>
void asCore::Interpolator::d_NLPotForce(const unsigned& gridId,
		const unsigned& recId, const unsigned& ligId,
		const unsigned& numDOFs,
		const Comp3_HD<float, DEVONLY>* LigPosTr,
		Comp5_HD<float, DEVONLY>* outLigPotForce,
		const cudaStream_t &stream)
{
	cudaVerifyKernel((
		asCore::d_NLPotForce<NLOnly><<<_gridSize, _blockSize, 0, stream>>>(gridId,
				recId, ligId, numDOFs,
				LigPosTr->d_x(),
				LigPosTr->d_y(),
				LigPosTr->d_z(),
				outLigPotForce->d_x(),
				outLigPotForce->d_y(),
				outLigPotForce->d_z(),
				outLigPotForce->d_v(),
				outLigPotForce->d_w())
	));
}

/* explicit instantiation */
template
void asCore::Interpolator::d_PotForce<asCore::built_in>(const unsigned& gridId,
		const unsigned& protId,
		const unsigned& numDOFs,
		const Comp3_HD<float, DEVONLY>* posTr,
		Comp5_HD<float, DEVONLY>* outPotForce,
		const cudaStream_t &stream);

template
void asCore::Interpolator::d_PotForce<asCore::manual>(const unsigned& gridId,
		const unsigned& protId,
		const unsigned& numDOFs,
		const Comp3_HD<float, DEVONLY>* posTr,
		Comp5_HD<float, DEVONLY>* outPotForce,
		const cudaStream_t &stream);

template
void asCore::Interpolator::d_NLPotForce<true>(const unsigned& gridId,
		const unsigned& recId, const unsigned& ligId,
		const unsigned& numDOFs,
		const Comp3_HD<float, DEVONLY>* LigPosTr,
		Comp5_HD<float, DEVONLY>* outLigPotForce,
		const cudaStream_t &stream);

template
void asCore::Interpolator::d_NLPotForce<false>(const unsigned& gridId,
		const unsigned& recId, const unsigned& ligId,
		const unsigned& numDOFs,
		const Comp3_HD<float, DEVONLY>* LigPosTr,
		Comp5_HD<float, DEVONLY>* outLigPotForce,
		const cudaStream_t &stream);

/****************************
 * protected member functions
 ****************************/

/****************************
 * private member functions
 ****************************/

