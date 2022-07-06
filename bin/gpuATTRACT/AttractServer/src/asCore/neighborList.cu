
#include "asCore/neighborList.h"
#include "as/asTypes.h"

extern __constant__ as::deviceProteinDesc c_Proteins[DEVICE_MAXPROTEINS];
extern __constant__ as::deviceGridUnionDesc c_Grids[DEVICE_MAXGRIDS];
extern __constant__ as::deviceParamTableDesc c_ParamTable;
extern __constant__ as::deviceSimParam c_SimParam;

__forceinline__ __host__ __device__ void LJPotForce(const float &dr2, const float &dr2_inv, const float &dx, const float &dy,	const float &dz,
		const as::AttrParamTable::type &LJParams, const float &swi, const as::AttrParamTable::PotShape &potShape,
		float &fx, float &fy, float &fz,
		float &energy
		) {

	const float &alen = LJParams.ac; // 24*eps
	const float &rlen = LJParams.rc; // sigma squared
	const float &ivor = LJParams.ipon;
	const float &rmin2 = LJParams.rmin2;
	const float &emin = LJParams.emin;

	const float rr23 = dr2_inv*dr2_inv*dr2_inv;

	float rrd, shapedelta;

	switch (potShape) {
	case as::AttrParamTable::PotShape::_8_6:
		rrd = dr2_inv;
		shapedelta = 2;
		break;
	case as::AttrParamTable::PotShape::_12_6:
		rrd = rr23;
		shapedelta = 6;
		break;
	default:
		break;
	}
	float rep = rlen * rrd;
	float vlj = (rep - alen) * rr23;
	if (dr2 < rmin2) {
		energy = swi * (vlj + (ivor - 1) * emin);
		float fb_swi = (6.0 * vlj + shapedelta * (rep * rr23))*swi;
		fx = fb_swi * dx;
		fy = fb_swi * dy;
		fz = fb_swi * dz;

	} else {
		float swi_ivor = swi * ivor;
		energy = swi_ivor * vlj;
		float ivor_fb_swi = swi_ivor*(6.0 * vlj + shapedelta * (rep * rr23));
		fx = ivor_fb_swi * dx;
		fy = ivor_fb_swi * dy;
		fz = ivor_fb_swi * dz;
	}
}

__forceinline__ __host__ __device__ void ChargePotForce(const float &dr2_inv, const float &dx, const float &dy,	const float &dz,
		const float &chargeLigRec,
		const float &swi,
		const dielec_t &dielec,
		float &fx, float &fy, float &fz, float &energy) {


	float dd;

	switch(dielec) {
	case constant:
		dd = sqrt(dr2_inv) - 1.0 / 50.0;
		break;
	case variable:
		dd = dr2_inv - 1.0 / (50.0 * 50.0);
		break;
	}

	/* (cap all distances at 50 A) */
	if (dd < 0)
		dd = 0;

//	std::cout <<"CHARGE "<< chargeLigRec << " " << dd << "\n";
	energy = swi * chargeLigRec * dd;

	switch (dielec) {
	case constant:
		if (dd <= 0) {
			fx = 0;
			fy = 0;
			fz = 0;
		} else {
			double et2;
			et2 = swi * chargeLigRec * sqrt(dr2_inv);
			fx = et2 * dx;
			fy = et2 * dy;
			fz = et2 * dz;
		}
		break;
	case variable:
		if (dd <= 0) {
			fx = 0;
			fy = 0;
			fz = 0;
		} else {
			double et2;
			et2 = swi * chargeLigRec * dr2_inv;
			fx = 2 * et2 * dx;
			fy = 2 * et2 * dy;
			fz = 2 * et2 * dz;
		}
		break;
	}

}


/*
 ** @brief: calculation of the switching function;
 */
inline __host__ __device__ float calcSwi(const float &dr2, const float &swiOn,
		const float &swiOff) {

	float swi;
	if (dr2 > swiOff * swiOff) {
		swi = 0;
	} else {
		float distance = sqrt(dr2);
		swi = 1 - (distance - swiOn) / (swiOff - swiOn);
	}
	return swi;
}

// Receptor Gradients not (yet) supported
template<bool NLOnly>
__global__ void asCore::d_NLPotForce(const unsigned gridId,
		const unsigned RecId, const unsigned LigId,
		const unsigned numDOFs,
		const float* LigPosX,
		const float* LigPosY,
		const float* LigPosZ,
		float* outLig_fx,
		float* outLig_fy,
		float* outLig_fz,
		float* outLigand_eEl,
		float* outLigand_eVdW)
{
	const unsigned i = blockDim.x * blockIdx.x + threadIdx.x;
	const unsigned LigNumEl = c_Proteins[LigId].nAtoms;
	if (i < LigNumEl*numDOFs) {
		const unsigned LigAttrIdx = i % LigNumEl;

		const unsigned atomTypeLig = c_Proteins[LigId].type[LigAttrIdx];

		if (atomTypeLig != 0) {


			const float posLigX = LigPosX[i];
			const float posLigY = LigPosY[i];
			const float posLigZ = LigPosZ[i];

			//DEBUG
//			if (i == 814)
//				printf("ATOM %u here\n", i);

			/* test if particle is out of bounds and perform data fetch and neigbourlist calculations */
			if (!(     (posLigX < c_Grids[gridId].NL.minDim.x || posLigX > c_Grids[gridId].NL.maxDim.x)
					|| (posLigY < c_Grids[gridId].NL.minDim.y || posLigY > c_Grids[gridId].NL.maxDim.y)
					|| (posLigZ < c_Grids[gridId].NL.minDim.z || posLigZ > c_Grids[gridId].NL.maxDim.z) )) {

				//DEBUG
//				if (i == 814)
//					printf("ATOM %u out of bounds test passed\n", i);

				const uint2 nDesc = tex3D<uint2>(c_Grids[gridId].NL.tex,
						(posLigX - c_Grids[gridId].NL.minDim.x) * c_Grids[gridId].NL.dVox_inv + 0.5f,
						(posLigY - c_Grids[gridId].NL.minDim.y) * c_Grids[gridId].NL.dVox_inv + 0.5f,
						(posLigZ - c_Grids[gridId].NL.minDim.z) * c_Grids[gridId].NL.dVox_inv + 0.5f);
				/* numEl = x; idx = y */

				//DEBUG
//				if (i == 814)
//					printf("ATOM %u List numEl startIdx %u %u \n", i, nDesc.x, nDesc.y);

				float3 fAcc = {0.0, 0.0, 0.0};
				float eVdWAcc = 0;
				float eElAcc = 0;
				for (unsigned j = 0; j < nDesc.x; ++j) {
					const unsigned nIdx = c_Grids[gridId].NL.neighborList[nDesc.y + j];

//					printf("ATOM %u RecAtomIdx %u \n", i, nIdx);


					float dx = posLigX - c_Proteins[RecId].xPos[nIdx];
					float dy = posLigY - c_Proteins[RecId].yPos[nIdx];
					float dz = posLigZ - c_Proteins[RecId].zPos[nIdx];
					const float dr2 = dx * dx + dy * dy + dz * dz;
					const float dPlateau2 = c_Grids[gridId].NL.dPlateau2;
					if ((dr2) > dPlateau2) {
						continue;
					}
//					printf("ATOM %u RecAtomIdx %u  under Plateau\n", i, nIdx);

					const float dr2_inv = 1.0/dr2; // inverse of dr2

					// Scale distances
					dx *= dr2_inv;
					dy *= dr2_inv;
					dz *= dr2_inv;

					float3 fVdW;
					float eVdW;

					const size_t atomTypeRec = c_Proteins[RecId].type[nIdx];

					// Calculate Switching Function -> not supported
//					float swi = 1;
//					float swiPlateau = 1;
//					if (c_SimParam.useSwi) {
//						swi = calcSwi(dr2, c_SimParam.swiOn, c_SimParam.swiOff);
//						swiPlateau = calcSwi(grid->dPlateau2(), c_SimParam.swiOn, c_SimParam.swiOff);
//					}


					// calculate energy and potential/energy of LJ/VdW potential

					as::AttrParamTable::type params = c_ParamTable.getParams(atomTypeRec-1, atomTypeLig-1);
					LJPotForce(dr2, dr2_inv, dx, dy, dz,
							params,
							1, c_ParamTable.shape,
							fVdW.x, fVdW.y, fVdW.z, eVdW);

					fAcc.x  += fVdW.x;
					fAcc.y  += fVdW.y;
					fAcc.z  += fVdW.z;
					eVdWAcc += eVdW;

					const float chargeLig = c_Proteins[LigId].charge[LigAttrIdx];
					const float chargeRec = c_Proteins[RecId].charge[nIdx];
					const float chargeLigRec = chargeLig * chargeRec * c_SimParam.ffelec;

					const bool calc_elec = abs(chargeLigRec) > 0.001; // evaluate electric potential

					float rdx, rdy, rdz;
					float dPlateau2_inv;
					if (!NLOnly || calc_elec) {
						dPlateau2_inv = 1/c_Grids[gridId].NL.dPlateau2;
						float ratio = sqrt(dr2*dPlateau2_inv);
						rdx = ratio*dx;
						rdy = ratio*dy;
						rdz = ratio*dz;
					}

					if (!NLOnly) {
						LJPotForce(dPlateau2, dPlateau2_inv, rdx, rdy, rdz,
							params,
							1, c_ParamTable.shape,
							fVdW.x, fVdW.y, fVdW.z, eVdW);
						fAcc.x  -= fVdW.x;
						fAcc.y  -= fVdW.y;
						fAcc.z  -= fVdW.z;
						eVdWAcc -= eVdW;
					}

					if (calc_elec) {
						float eEl;
						float3 fEl;

						// calculate energy and potential/energy of charge potential

						ChargePotForce(dr2_inv, dx, dy, dz,
								chargeLigRec,
								1, c_SimParam.dielec,
								fEl.x, fEl.y, fEl.z, eEl);

						fAcc.x += fEl.x;
						fAcc.y += fEl.y;
						fAcc.z += fEl.z;
						eElAcc += eEl;

						ChargePotForce(dPlateau2_inv, rdx, rdy, rdz,
								chargeLigRec,
								1, c_SimParam.dielec,
								fEl.x, fEl.y, fEl.z, eEl);
						fAcc.x -= fEl.x;
						fAcc.y -= fEl.y;
						fAcc.z -= fEl.z;
						eElAcc -= eEl;

					}
				}

				/* store results back to global memory */
				if (nDesc.x > 0) {
					outLig_fx[i] += fAcc.x;
					outLig_fy[i] += fAcc.y;
					outLig_fz[i] += fAcc.z;
					outLigand_eEl[i] += eElAcc;
					outLigand_eVdW[i] += eVdWAcc;
				}
			}
		} // if (atomtype != 0)
	}
}

void asCore::h_NLPotForce(const as::NLGrid *grid,
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
		float* outRec_fz)
{
	//DEBUG
//	nlTotCount = grid->neighborListSize();

	const unsigned nAtomsLig = lig->nAtoms();
	/* loop over all elements in input/output */
	for (unsigned i = 0; i < nAtomsLig; ++i) {
		const unsigned atomTypeLig = lig->type()[i];

		if (atomTypeLig == 0)
			continue;

		const float posLigX = LigPosX[i];
		const float posLigY = LigPosY[i];
		const float posLigZ = LigPosZ[i];

//		//DEBUG
//		if (i == 814)
//			printf("ATOM %u here\n", i);

		/* test if particle is out of bounds and perform data fetch and neigbourlist calculations */
		if (!(grid->outOfBounds(posLigX, posLigY, posLigZ))) {

//			//DEBUG
//			if (i == 814)
//				printf("ATOM %u %f %f %f out of bounds test passed\n", i, posLigX, posLigY, posLigZ);

			int idxX, idxY, idxZ;
			grid->getIndex(posLigX, posLigY, posLigZ, idxX, idxY, idxZ);
			const NeighbourDesc &nDesc = grid->getNeighbourDesc(idxX, idxY,
					idxZ);

			//DEBUG
//			nlCount += nDesc.numEl;
//			//DEBUG
//			if (i == 814)
//				printf("ATOM %u List numEl startIdx %u %u \n", i, nDesc.numEl, nDesc.idx);

			for (unsigned j = 0; j < nDesc.numEl; ++j) {
				const unsigned nIdx = grid->getNeighbor(nDesc.idx + j);

				float dx, dy, dz;
				if (rec->numModes() > 0) {
					dx = posLigX - RecPosX[nIdx];
					dy = posLigY - RecPosY[nIdx];
					dz = posLigZ - RecPosZ[nIdx];
				} else {
					dx = posLigX - rec->xPos(0)[nIdx];
					dy = posLigY - rec->yPos(0)[nIdx];
					dz = posLigZ - rec->zPos(0)[nIdx];
				}

				const float dr2 = dx * dx + dy * dy + dz * dz;

				if (grid->outOfPlateau(dr2)) {
					continue;
				}

				const float dr2_inv = 1.0/dr2; // inverse of dr2

				// Scale distances
				dx *= dr2_inv;
				dy *= dr2_inv;
				dz *= dr2_inv;

				float3 fVdW;
				float3 fAcc = {0.0,0.0,0.0};
				float eVdW;
				float eVdWAcc = 0;

				const size_t atomTypeRec = rec->type()[nIdx];

				// Calculate Switching Function
				float swi = 1;
				float swiPlateau = 1;
				// ToDo: this is not equivalent to chargenonzero > 0.001 in ATTRACT grid.cpp! Why only if charge of Ligand > 0.001 in attract!
//				if (simParam->useSwi) {
//					swi = calcSwi(dr2, simParam->swiOn, simParam->swiOff);
//					swiPlateau = calcSwi(grid->dPlateau2(), simParam->swiOn, simParam->swiOff);
//				}

				// calculate energy and potential/energy of LJ/VdW potential
				assert(atomTypeRec > 0);
				assert(atomTypeRec < 99);
				assert(atomTypeLig > 0);
				assert(atomTypeLig < 99);
				LJPotForce(dr2, dr2_inv, dx, dy, dz,
						table->getParams(atomTypeRec-1, atomTypeLig-1),
						swi, table->potShape(),
						fVdW.x, fVdW.y, fVdW.z, eVdW);

				fAcc.x  += fVdW.x;
				fAcc.y  += fVdW.y;
				fAcc.z  += fVdW.z;
				eVdWAcc += eVdW;

				if (simParam->useRecGrad) {
					outRec_fx[nIdx] -= fVdW.x;
					outRec_fy[nIdx] -= fVdW.y;
					outRec_fz[nIdx] -= fVdW.z;
				}

				const float chargeLig = lig->charge()[i];
				const float chargeRec = rec->charge()[nIdx];

				const float chargeLigRec = chargeLig * chargeRec * simParam->ffelec;

				const bool calc_elec = abs(chargeLigRec) > 0.001; // evaluate electric potential

//				std::cout << calc_elec << std::endl;
				float rdx, rdy, rdz;
				if (simParam->usePot || calc_elec) {

					float ratio = grid->getRatio(dr2);
					rdx = ratio*dx;
					rdy = ratio*dy;
					rdz = ratio*dz;
				}

				if (simParam->usePot) {
					LJPotForce(grid->dPlateau2(), grid->dPlateau2_inv(), rdx, rdy, rdz,
						table->getParams(atomTypeRec-1, atomTypeLig-1),
						swiPlateau, table->potShape(),
						fVdW.x, fVdW.y, fVdW.z, eVdW);
					fAcc.x  -= fVdW.x;
					fAcc.y  -= fVdW.y;
					fAcc.z  -= fVdW.z;
					eVdWAcc -= eVdW;
					if (simParam->useRecGrad) {
						outRec_fx[nIdx] += fVdW.x;
						outRec_fy[nIdx] += fVdW.y;
						outRec_fz[nIdx] += fVdW.z;
					}
				}

				float eElAcc = 0;
				if (calc_elec) {

					float eEl;
					float4 fEl;

					// calculate energy and potential/energy of charge potential
					ChargePotForce(dr2_inv, dx, dy, dz,
							chargeLigRec,
							swi, simParam->dielec,
							fEl.x, fEl.y, fEl.z, eEl);

					fAcc.x += fEl.x;
					fAcc.y += fEl.y;
					fAcc.z += fEl.z;
					eElAcc += eEl;

					if (simParam->useRecGrad) {
						outRec_fx[nIdx] -= fEl.x;
						outRec_fy[nIdx] -= fEl.y;
						outRec_fz[nIdx] -= fEl.z;
					}

					ChargePotForce(grid->dPlateau2_inv(), rdx, rdy, rdz,
							chargeLigRec,
							swiPlateau, simParam->dielec,
							fEl.x, fEl.y, fEl.z, eEl);
					fAcc.x -= fEl.x;
					fAcc.y -= fEl.y;
					fAcc.z -= fEl.z;
					eElAcc -= eEl;

					if (simParam->useRecGrad) {
						outRec_fx[nIdx] += fEl.x;
						outRec_fy[nIdx] += fEl.y;
						outRec_fz[nIdx] += fEl.z;
					}
				}

				// TODO: put store outside loop

				outLig_fx[i]   +=  fAcc.x;
				outLig_fy[i]   +=  fAcc.y;
				outLig_fz[i]   +=  fAcc.z;
				outLig_eEl[i]  +=  eElAcc;
				outLig_eVdW[i] +=  eVdWAcc;
			}
		}
	}
}

template __global__ void
asCore::d_NLPotForce<true>(const unsigned gridId,
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

template __global__ void
asCore::d_NLPotForce<false>(const unsigned gridId,
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


