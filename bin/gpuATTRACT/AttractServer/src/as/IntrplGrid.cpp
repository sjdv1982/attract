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
#include "cuda_runtime.h"

#include "as/IntrplGrid.h"
#include "asUtils/helper.h"

/* Constructor */
as::IntrplGrid::IntrplGrid(GradEnGridDesc desc):
		Grid::Grid(desc.width, desc.height,	desc.depth,
				make_float3(desc.posMin[0], desc.posMin[1], desc.posMin[2]),
				desc.gridSpacing),
		_grid(desc.grid),
		_numGrids(desc.numGrids),
		_voxelVol(_dVox * _dVox * _dVox),
		_voxelVol_inv(1.f/_voxelVol),
		_dVox_inv(1.f/_dVox)

{
	_maxDim.x = pos().x + (_width  - 1) * _dVox;
	_maxDim.y = pos().y + (_height - 1) * _dVox;
	_maxDim.z = pos().z + (_depth  - 1) * _dVox;
}

/* Destructor */
as::IntrplGrid::~IntrplGrid() {
	freeHost();
}


/****************************
 * public member functions
 ****************************/

/*
 ** @brief: returns the size of the grid in terms of memory consumption
 */
double as::IntrplGrid::mySize(asUtils::sizeType_t type) const {
	unsigned long long sizeInByte;
	sizeInByte = (unsigned long long) _numGrids*_width*_height*_depth*sizeof(float4);

	return sizeConvert(sizeInByte, type);
}

void as::IntrplGrid::freeHost() {
	delete[] _grid;
}

/****************************
 * protected member functions
 ****************************/

/****************************
 * private member functions
 ****************************/


/* associated functions */
void as::print(as::IntrplGrid* desc, uint grid0, uint gridM,
		uint comp0, uint compM,
		uint x0, uint xM, uint y0, uint yM, uint z0, uint zM)
{
	// Store io flags
	using namespace std;
	int precisionSetting = cout.precision( );
	ios::fmtflags flagSettings = cout.flags();

	cout.setf(ios::fixed | ios::showpos | ios::showpoint);

	cout << "numGrids " << desc->numTypes() << endl;
	cout << "width " << desc->width() << endl;
	cout << "height " << desc->height() << endl;
	cout << "depth " << desc->depth() << endl;
	cout << "dVox " << desc->dVox() << endl;
	cout << "minPos " << desc->pos().x << " " << desc->pos().y << " " << desc->pos().z << endl;

//	cout << "Debug"<< " " <<   grid0 << " " <<  gridM << " " <<  comp0 << " " <<  compM << " " <<  x0 << " " <<  xM << " " <<  y0 << " " <<  yM << " " <<  z0 << " " <<  zM  << endl;

	cout.precision(6);
	uint xOff = 5;
	uint yOff = 2;
	uint width = 13;
	uint height = 2; // > 0 !!!
	for (uint grid = grid0; grid <= gridM; ++grid) {
		float4* grid_ptr = desc->getHostGridPtr(grid);
		cout << "###GRID " << grid <<"###" << endl;
		for (uint comp = comp0; comp <= compM; ++comp) {
			cout << "###COMP " << comp <<"###" << endl;
			for (uint z = z0; z <= zM; ++z) {
				cout << "---Slice" << z << "---" << endl;
				cout << setw(xOff) << " ";
				for (uint x = x0; x <= xM; ++x)
					cout << setw(width)<<  x;
				for (uint i = 0; i < yOff; ++i)
					cout << endl;

				for (uint y = y0; y <= yM; ++y) {
					cout << setw(xOff) <<  y;
					for (uint x = x0; x <= xM; ++x) {
						uint idx =  desc->width()*(z*desc->height() + y) + x;
						float* ptr = (float*)&grid_ptr[idx];
						cout << setw(width) << ptr[comp];
					}
					for (uint i = 0; i < height; ++i)
						cout << endl;
				}
				cout << endl;
			}
		}
	}

	// Restore io flags
	cout.precision(precisionSetting);
	cout.flags(flagSettings);

}

void as::print(const GradEnGridDesc& desc, uint grid0, uint gridM, uint comp0, uint compM, uint x0, uint xM, uint y0, uint yM, uint z0, uint zM)
{
	// Store io flags
	using namespace std;
	int precisionSetting = cout.precision( );
	ios::fmtflags flagSettings = cout.flags();

	cout.setf(ios::fixed | ios::showpos | ios::showpoint);

	cout << "numGrids " << desc.numGrids << endl;
	cout << "width " << desc.width << endl;
	cout << "height " << desc.height << endl;
	cout << "depth " << desc.depth << endl;
	cout << "dVox " << desc.gridSpacing << endl;
	cout << "minPos " << desc.posMin[0] << " " << desc.posMin[1] << " " << desc.posMin[2] << endl;

//	cout << "Debug"<< " " <<   grid0 << " " <<  gridM << " " <<  comp0 << " " <<  compM << " " <<  x0 << " " <<  xM << " " <<  y0 << " " <<  yM << " " <<  z0 << " " <<  zM  << endl;

	cout.precision(6);
	uint xOff = 5;
	uint yOff = 2;
	uint width = 13;
	uint height = 2; // > 0 !!!
	for (uint grid = grid0; grid <= gridM; ++grid) {
		float4* grid_ptr = desc.grid + grid*desc.width*desc.height*desc.depth;
		cout << "###GRID " << grid <<"###" << endl;
		for (uint comp = comp0; comp <= compM; ++comp) {
			cout << "###COMP " << comp <<"###" << endl;
			for (uint z = z0; z <= zM; ++z) {
				cout << "---Slice" << z << "---" << endl;
				cout << setw(xOff) << " ";
				for (uint x = x0; x <= xM; ++x)
					cout << setw(width)<<  x;
				for (uint i = 0; i < yOff; ++i)
					cout << endl;

				for (uint y = y0; y <= yM; ++y) {
					cout << setw(xOff) <<  y;
					for (uint x = x0; x <= xM; ++x) {
						uint idx =  desc.width*(z*desc.height + y) + x;
						float* ptr = (float*)&grid_ptr[idx];
						cout << setw(width) << ptr[comp];
					}
					for (uint i = 0; i < height; ++i)
						cout << endl;
				}
				cout << endl;
			}
		}
	}

	// Restore io flags
	cout.precision(precisionSetting);
	cout.flags(flagSettings);

}



