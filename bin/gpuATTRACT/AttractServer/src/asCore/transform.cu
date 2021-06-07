
#include "asCore/transform.h"
#include "config.h"

extern __constant__ as::deviceProteinDesc c_Proteins[DEVICE_MAXPROTEINS];

__global__ void asCore::d_DOF2Pos(const unsigned protId,
		const unsigned numDOFs, const as::DOF* dofs,
		unsigned short *confTr,
		float* xTr, float* yTr, float* zTr)

{
	/* calculate element index that is to be prcessed */
	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;
	/* get number of ligand atoms from constant memory */
	const unsigned nAtoms = c_Proteins[protId].nAtoms;	

	if (idx < nAtoms*numDOFs) {
		/* load DOF from global memory */
		unsigned DOFidx = idx / nAtoms;
		as::DOF dof = dofs[DOFidx];

		if (idx % nAtoms == 0) {
			confTr[DOFidx] = dof.conf;
		}

		const unsigned int conf_base = dof.conf * c_Proteins[protId].nAtoms;

		/* load original position from global memory */
		float pos[3];
		pos[0] = c_Proteins[protId].xPos[conf_base + idx % nAtoms];
		pos[1] = c_Proteins[protId].yPos[conf_base + idx % nAtoms];
		pos[2] = c_Proteins[protId].zPos[conf_base + idx % nAtoms];

		/* calculate rotation matrix */
		//first rotate using rot
		//		this is a rotation around the Z axis
		//			rotating Y into X and X into -Y
		//then rotate using ssi
		//		this is a rotation around the (new) Y axis
		//			rotating X into -Z and Z into X
		//finally, rotate using phi
		//		this is a rotation around the (new) Z axis
		//			rotating X into Y and Y into -X
		const float cSSI = cos(dof.ang.y);
		const float cPHI = cos(dof.ang.x);
		const float sSSI = sin(dof.ang.y);
		const float sPHI = sin(dof.ang.x);
		const float cROT = cos(dof.ang.z);
		const float sROT = sin(dof.ang.z);
		float rotmat[9];


		rotmat[0] = cROT * cSSI * cPHI + sROT * sPHI;
		rotmat[1] = sROT * cSSI * cPHI - cROT * sPHI;
		rotmat[2] = sSSI * cPHI;

		rotmat[3] = cROT * cSSI * sPHI - sROT * cPHI;
		rotmat[4] = sROT * cSSI * sPHI + cROT * cPHI;
		rotmat[5] = sSSI * sPHI;

		rotmat[6] = -cROT * sSSI;
		rotmat[7] = -sROT * sSSI;
		rotmat[8] = cSSI;

		/* apply rotation matrix */
		const float x_tmp = pos[0];
		const float y_tmp = pos[1];
		const float z_tmp = pos[2];
		pos[0] = rotmat[0] * x_tmp + rotmat[1] * y_tmp + rotmat[2] * z_tmp;
		pos[1] = rotmat[3] * x_tmp + rotmat[4] * y_tmp + rotmat[5] * z_tmp;
		pos[2] = rotmat[6] * x_tmp + rotmat[7] * y_tmp + rotmat[8] * z_tmp;

		/* translate position */
		pos[0] += dof.pos.x;
		pos[1] += dof.pos.y;
		pos[2] += dof.pos.z;

		/* store position to global memory */
		xTr[idx] = pos[0];
		yTr[idx] = pos[1];
		zTr[idx] = pos[2];
	}
}

__global__ void asCore::d_DOF2Pos_modes(const unsigned protId,
		const unsigned numDOFs, const as::DOF* dofs,
		unsigned short *confTr,
		float* xTr, float* yTr, float* zTr,
		unsigned short *confDef,
		float* xDef, float* yDef, float* zDef)

{

	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;
	const unsigned nAtoms = c_Proteins[protId].nAtoms;
	if (idx < nAtoms*numDOFs) {

		unsigned DOFidx = idx / nAtoms;
		const as::DOF dof = dofs[DOFidx];

		if (idx % nAtoms == 0) {
			confTr[DOFidx] = dof.conf;
			confDef[DOFidx] = dof.conf;
		}
		const unsigned int conf_base = dof.conf * c_Proteins[protId].nAtoms;

		/* load original position from global memory */
		float pos[3];
		pos[0] = c_Proteins[protId].xPos[conf_base + idx % nAtoms];
		pos[1] = c_Proteins[protId].yPos[conf_base + idx % nAtoms];
		pos[2] = c_Proteins[protId].zPos[conf_base + idx % nAtoms];
		
		/* deform protein if necessary */
		const unsigned numModes = c_Proteins[protId].numModes;
		for (uint mode = 0; mode < numModes; ++mode) {
			const float scale = dof.modes[mode];
			const unsigned modeIdx = idx* numModes + mode;
			pos[0] += scale * c_Proteins[protId].xModes[modeIdx];
			pos[1] += scale * c_Proteins[protId].yModes[modeIdx];
			pos[2] += scale * c_Proteins[protId].zModes[modeIdx];
		}

		xDef[idx] = pos[0];
		yDef[idx] = pos[1];
        zDef[idx] = pos[2];

		//first rotate using rot
		//		this is a rotation around the Z axis
		//			rotating Y into X and X into -Y
		//then rotate using ssi
		//		this is a rotation around the (new) Y axis
		//			rotating X into -Z and Z into X
		//finally, rotate using phi
		//		this is a rotation around the (new) Z axis
		//			rotating X into Y and Y into -X
		const float cSSI = cos(dof.ang.y);
		const float cPHI = cos(dof.ang.x);
		const float sSSI = sin(dof.ang.y);
		const float sPHI = sin(dof.ang.x);
		const float cROT = cos(dof.ang.z);
		const float sROT = sin(dof.ang.z);

		float rotmat[9];
		rotmat[0] = cROT * cSSI * cPHI + sROT * sPHI;
		rotmat[1] = sROT * cSSI * cPHI - cROT * sPHI;
		rotmat[2] = sSSI * cPHI;

		rotmat[3] = cROT * cSSI * sPHI - sROT * cPHI;
		rotmat[4] = sROT * cSSI * sPHI + cROT * cPHI;
		rotmat[5] = sSSI * sPHI;

		rotmat[6] = -cROT * sSSI;
		rotmat[7] = -sROT * sSSI;
		rotmat[8] = cSSI;

		const float x_tmp = pos[0];
		const float y_tmp = pos[1];
		const float z_tmp = pos[2];
		pos[0] = rotmat[0] * x_tmp + rotmat[1] * y_tmp + rotmat[2] * z_tmp;
		pos[1] = rotmat[3] * x_tmp + rotmat[4] * y_tmp + rotmat[5] * z_tmp;
		pos[2] = rotmat[6] * x_tmp + rotmat[7] * y_tmp + rotmat[8] * z_tmp;

		pos[0] += dof.pos.x;
		pos[1] += dof.pos.y;
		pos[2] += dof.pos.z;

		xTr[idx] = pos[0];
		yTr[idx] = pos[1];
		zTr[idx] = pos[2];
	}
}

__global__ void asCore::d_DOF2Deform(const unsigned protId,
		const as::DOF* dofs, const unsigned numDOFs,
		float* xTrDef, float* yTrDef, float* zTrDef)
{
	const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;
	const unsigned nAtoms = c_Proteins[protId].nAtoms;
	if (idx < nAtoms*numDOFs) {

		unsigned DOFidx = idx / nAtoms;
		const as::DOF dof = dofs[DOFidx];

		const unsigned int conf_base = dof.conf * c_Proteins[protId].nAtoms;

		/* load original position from global memory */
		float pos[3];
		pos[0] = c_Proteins[protId].xPos[conf_base + idx % nAtoms];
		pos[1] = c_Proteins[protId].yPos[conf_base + idx % nAtoms];
		pos[2] = c_Proteins[protId].zPos[conf_base + idx % nAtoms];


		/* deform protein */
		unsigned numModes = c_Proteins[protId].numModes;
		for (uint mode = 0; mode < numModes; ++mode) {
			float scale = dof.modes[mode];
			unsigned modeIdx = idx* numModes + mode;
			pos[0] += scale * c_Proteins[protId].xModes[modeIdx];
			pos[1] += scale * c_Proteins[protId].yModes[modeIdx];
			pos[2] += scale * c_Proteins[protId].zModes[modeIdx];
		}

		xTrDef[idx] = pos[0];
		yTrDef[idx] = pos[1];
        zTrDef[idx] = pos[2];
		
	}
}
