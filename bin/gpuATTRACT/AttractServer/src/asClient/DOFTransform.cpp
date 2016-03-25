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

#include <cmath>
#include <cassert>

#include <asClient/DOFTransform.h>
#include <asUtils/cudaMath.h>

void asClient::transformDOF_glob2rec(const std::vector<as::DOF>& dof_rec, std::vector<as::DOF>& dof_lig,
			const asUtils::Vec3f& pivot_rec, const asUtils::Vec3f& pivot_lig,
			bool centered_rec, bool centered_lig) {
	using namespace std;
	using namespace as;
	using namespace asUtils;

	assert(dof_rec.size() == dof_lig.size());

	// center ligand dofs by the respective ligand pivot
	if (centered_lig == false) {
		for (auto& dof : dof_lig) {
			dof.pos.x += pivot_lig[0];
			dof.pos.y += pivot_lig[1];
			dof.pos.z += pivot_lig[2];
		}
	}

	/* shift ligand dofs by receptor pivot */
	if (centered_rec == false) {
		if (!(pivot_rec == Vec3f(0.0f))) {
			for (auto& dof : dof_lig) {
				dof.pos.x -= pivot_rec[0];
				dof.pos.y -= pivot_rec[1];
				dof.pos.z -= pivot_rec[2];
			}
		}
	}

	/* rotate ligand into receptor frame and shift ligand by receptor dofs*/
	for (unsigned j = 0; j < dof_lig.size(); ++j) {
		const float3& pos_rec = dof_rec[j].pos;
		if (pos_rec.x != 0.0f || pos_rec.y != 0.0f || pos_rec.z != 0.0f ) {
			using namespace asUtils;
			float3& pos_lig = dof_lig[j].pos;
			pos_lig = pos_lig - pos_rec;
		}
		const float3& ang_rec = dof_rec[j].ang;
		if (ang_rec.x != 0.0f || ang_rec.y != 0.0f || ang_rec.z != 0.0f ) {

			float3& ang_lig = dof_lig[j].ang;
			float3& pos_lig = dof_lig[j].pos;

			RotMatf mat_rec_inv;
			euler2rotmat(-ang_rec.x, -ang_rec.y, -ang_rec.z, mat_rec_inv);
			RotMatf mat_lig;
			euler2rotmat(ang_lig.x, ang_lig.y, ang_lig.z, mat_lig);

//				RotMatf mat_rec_inv = mat_rec.getInv();
			Vec3f pos_lig_v(pos_lig.x, pos_lig.y, pos_lig.z);
			pos_lig_v = mat_rec_inv * pos_lig_v;
			pos_lig = make_float3(pos_lig_v[0], pos_lig_v[1], pos_lig_v[2]);
			RotMatf mat =  mat_rec_inv *mat_lig ;
			rotmat2euler(mat, ang_lig.x, ang_lig.y, ang_lig.z);
		}
	}
	/* the receptor dofs can now considered to be zero */

}

void asClient::transformEnGrad_rec2glob(const std::vector<as::DOF>& dof_rec, std::vector<as::EnGrad>& enGrad_lig) {
	using namespace std;
	using namespace as;
	using namespace asUtils;

	assert(dof_rec.size() == enGrad_lig.size());

	/* rotate forces and torques in global frame */
	for (unsigned j = 0; j < dof_rec.size(); ++j) {
		const float3& ang_rec = dof_rec[j].ang;
		if (ang_rec.x != 0.0f || ang_rec.y != 0.0f || ang_rec.z != 0.0f ) {
			float3& force_lig = enGrad_lig[j].pos;
			Vec3f force_lig_v(force_lig.x, force_lig.y, force_lig.z);

			/* rotate forces */
			RotMatf mat_rec;
			euler2rotmat(ang_rec.x, ang_rec.y, ang_rec.z, mat_rec);

			force_lig_v = mat_rec * force_lig_v;
			force_lig = make_float3(force_lig_v[0], force_lig_v[1], force_lig_v[2]);

			/* rotate torques */


		}
	}

}

template<typename T>
void asClient::rotmat2euler(const asUtils::RotMat<T>& rotmat, T& phi, T& ssi, T& rot) {
	phi = atan2(rotmat[5], rotmat[2]);
	ssi = acos(fmin(fmax(rotmat[8],-1.0), 1.0));
	rot = atan2(-rotmat[7], -rotmat[6]);
	/* handel gimbal lock */
	if (abs(rotmat[8] >= 0.9999)) {
		phi = 0.0;
		if(abs(rotmat[0]) >= 0.9999) {
			if(rotmat[0] < 0.0) {
				rot = M_PI;
			} else {
				rot = 0.0;
			}

			if(rotmat[8] < 0.0) {
				ssi = M_PI;
			} else {
				ssi = 0.0;
			}

		} else {
			if(rotmat[8] < 0.0) {
				ssi = M_PI;
				rot = -acos(-rotmat[0]);
			} else {
				ssi = 0.0;
				rot = acos(rotmat[0]);
			}

		}

		if (rotmat[1] < 0) {
			rot *= -1.0;
		}
	}
}

template<typename T>
void asClient::euler2rotmat(const T& phi, const T& ssi, const T& rot, asUtils::RotMat<T>& rotmat) {
	//first rotate using rot
	//		this is a rotation around the Z axis
	//			rotating Y into X and X into -Y
	//then rotate using ssi
	//		this is a rotation around the (new) Y axis
	//			rotating X into -Z and Z into X
	//finally, rotate using phi
	//		this is a rotation around the (new) Z axis
	//			rotating X into Y and Y into -X

	T cSSI = cos(ssi);
	T cPHI = cos(phi);
	T sSSI = sin(ssi);
	T sPHI = sin(phi);
	T cSSI_cPHI = cSSI * cPHI;
	T cSSI_sPHI = cSSI * sPHI;
	T sSSI_cPHI = sSSI * cPHI;
	T sSSI_sPHI = sSSI * sPHI;
	T cROT = cos(rot);
	T sROT = sin(rot);

	rotmat[0] = cROT * cSSI_cPHI + sROT * sPHI;
	rotmat[1] = sROT * cSSI_cPHI - cROT * sPHI;
	rotmat[2] = sSSI_cPHI;

	rotmat[3] = cROT * cSSI_sPHI - sROT * cPHI;
	rotmat[4] = sROT * cSSI_sPHI + cROT * cPHI;
	rotmat[5] = sSSI_sPHI;

	rotmat[6] = -cROT * sSSI;
	rotmat[7] = -sROT * sSSI;
	rotmat[8] = cSSI;
}

/* explicit instantiation */

template void asClient::rotmat2euler<float>(const asUtils::RotMat<float>& rotmat, float& phi, float& ssi, float& rot);
template void asClient::rotmat2euler<double>(const asUtils::RotMat<double>& rotmat, double& phi, double& ssi, double& rot);

template void asClient::euler2rotmat<float>(const float& phi, const float& ssi, const float& rot, asUtils::RotMat<float>& rotmat);
template void asClient::euler2rotmat<double>(const double& phi, const double& ssi, const double& rot, asUtils::RotMat<double>& rotmat);


