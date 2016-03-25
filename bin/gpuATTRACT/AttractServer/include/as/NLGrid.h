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

#ifndef NLGRID_H_
#define NLGRID_H_

#include <cmath>
#include <cassert>
#include <iostream>
#include "as/interface.h"
#include "as/Grid.h"

namespace as {


class NLGrid : public as::Grid {
public:
	/* Constructor */
	__host__ NLGrid(NLGridDesc desc);

	/* Destructor */
	__host__ virtual ~NLGrid();

	/***************
	* G E T T E R
	***************/

	float dPlateau2() const {
		return _dPlateau2;
	}

	float dPlateau2_inv() const {
		return _dPlateau2_inv;
	}

	unsigned* neighborList() const {
		return _neighborList;
	}

	unsigned neighborListSize() const {
		return _numElInLists;
	}

	float dVox_inv() const {
		return _dVox_inv;
	}

	float3 minDim() const {
		return Grid::pos();
	}

	float3 maxDim() const {
		return _maxDim;
	}

	NeighbourDesc* grid() const {
		return _grid;
	}

	/***************
	* S E T T E R
	***************/
	virtual void setPos(float3 pos) {
		_pos = pos;
		_maxDim.x = _pos.x + (_width  - 1) * _dVox;
		_maxDim.y = _pos.y + (_height - 1) * _dVox;
		_maxDim.z = _pos.z + (_depth  - 1) * _dVox;
	}

	/****************************
	 * public member functions
	 ****************************/

	/*
	 ** @brief: return the min grid indices according to the position
	 */
	__host__ void getIndex(const float &x, const float &y, const float &z, int &idxX, int &idxY, int &idxZ) const {

		idxX = round((x - _pos.x)*_dVox_inv);
		idxY = round((y - _pos.y)*_dVox_inv);
		idxZ = round((z - _pos.z)*_dVox_inv);
		assert(idxX >= 0 && idxX <= static_cast<int>(_width)-1);
		assert(idxY >= 0 && idxY <= static_cast<int>(_height)-1);
		assert(idxZ >= 0 && idxZ <= static_cast<int>(_depth)-1);
	}

	/*
	 ** @brief: returns the size of the grid in terms of memory consumption
	 */
	double mySize(asUtils::sizeType_t type = asUtils::byte) const;

	/*
	 ** @brief: check out of bounds according to indices
	 */
	bool outOfBounds(const int &idxX, const int &idxY, const int &idxZ) const {
		return (idxX < 0 || (idxX > static_cast<int>(_width))) ||
				(idxY < 0 || (idxY > static_cast<int>(_height))) ||
				(idxZ < 0 || (idxZ > static_cast<int>(_depth)));
	}

	/*
	 ** @brief: check out of bounds according to position
	 */
	bool outOfBounds(const float &x, const float &y, const float &z) const {
		return ((x < minDim().x || x > maxDim().x) ||
				(y < minDim().y || y > maxDim().y) ||
				(z < minDim().z || z > maxDim().z));
	}

	 const NeighbourDesc& getNeighbourDesc(const int &idxX, const int &idxY, const int &idxZ) const {

		return _grid[_width*(idxZ*_height + idxY) + idxX];
	}

	unsigned getNeighbor(const unsigned &idx) const {
		return _neighborList[idx];
	}

	bool outOfPlateau(const float &dr2) const {
		return dr2 > _dPlateau2;
	}

	float getRatio (float d2) const {
		return _ratio[int(d2*10000)];
	}

	/*
	 ** @brief: Free resources on the host.
	 */
	void freeHost();

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

	void initTexture();

	/****************************
	 * private member variables
	 ****************************/

	NeighbourDesc* _grid;	/** host grid of neighbourlists */

	unsigned _numElInLists;	/** number of total elements in Lists */
	unsigned* _neighborList;	/** contains receptor atom indices for all neighbours.
								It is basically a concatenation of subsequent neighbour lists.
								The exact sequence is the same as for the grid */

	float3 _maxDim;			/** max. grid dimension */
	float _dPlateau2;		/** plateau distance squared */
	float _dPlateau2_inv;	/** plateau distance squared inverse*/
	float _dVox_inv;		/** pre-calculated inverse of the voxel distance */

	float *_ratio;	/** precomputed square roots of the distance squared */
};

} // namespace

#endif /* NLGRID_H_ */
