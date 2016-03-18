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

#include "as/NLGrid.h"
#include "asUtils/helper.h"

/* Constructor */
as::NLGrid::NLGrid(NLGridDesc desc):
		Grid::Grid(desc.width, desc.height,	desc.depth,
				make_float3(desc.posMin[0], desc.posMin[1], desc.posMin[2]),
				desc.gridSpacing),
		_grid(desc.grid),
		_numElInLists(desc.numEl),
		_neighborList(desc.neighborArray),
		_dPlateau2(desc.dPlateau*desc.dPlateau),
		_dPlateau2_inv(1/_dPlateau2),
		_dVox_inv(1/_dVox)
{
	_maxDim.x = pos().x + (_width  - 1) * _dVox;
	_maxDim.y = pos().y + (_height - 1) * _dVox;
	_maxDim.z = pos().z + (_depth  - 1) * _dVox;

	/* init ratios (taken from the original attract code) */
	int size  = int(10000*_dPlateau2);
	_ratio = new float[size+1];

	for (int n = 0; n <= size; n++) {
		double d2 = ((n + 0.5) / 10000);
		_ratio[n] = sqrt(d2 / _dPlateau2);
	}
}

/* Destructor */
as::NLGrid::~NLGrid() {
	freeHost();
}


/****************************
 * public member functions
 ****************************/
void as::NLGrid::freeHost() {
	delete[] _neighborList;
	delete[] _grid;
	delete[] _ratio;
}

double as::NLGrid::mySize(asUtils::sizeType_t type) const {

	unsigned long long sizeInByte;
	sizeInByte = (unsigned long long) (_numElInLists*sizeof(uint) +
			_width*_height*_depth*sizeof(NeighbourDesc));

	return sizeConvert(sizeInByte, type);
}

/****************************
 * protected member functions
 ****************************/

/****************************
 * private member functions
 ****************************/


