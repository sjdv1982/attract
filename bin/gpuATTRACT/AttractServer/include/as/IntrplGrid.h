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

#ifndef INTRPLGRID_H_
#define INTRPLGRID_H_

#include <cassert>
#include <cmath>
#include <iostream>

#include "as/Grid.h"
#include "as/interface.h"
#include "as/asTypes.h"
#include "asUtils/macros.h"

namespace as {
	class IntrplGrid;

	void print(IntrplGrid* desc, uint grid0, uint gridM,
			uint comp0, uint compM,
			uint x0, uint xM, uint y0, uint yM, uint z0, uint zM);

	void print(const GradEnGridDesc& desc, uint grid0, uint gridM,
			uint comp0, uint compM,
			uint x0, uint xM, uint y0, uint yM, uint z0, uint zM);


class IntrplGrid : public as::Grid {
public:
	/* Constructor */
	__host__ IntrplGrid(GradEnGridDesc desc);

	/* Destructor */
	__host__ virtual ~IntrplGrid();

	/* Types */
	enum gridType_t {
		InnerGrid = 0,
		OuterGrid = 1
	};

	/***************
	* G E T T E R
	***************/
	float voxelVol() const {
		return _voxelVol;
	}

	float voxelVol_inv() const {
		return _voxelVol_inv;
	}

	unsigned numTypes() const {
		return _numGrids;
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

	/*
	 ** @brief: returns a pointer to the respective grid.
	 */
	float4* getHostGridPtr(uint i) const {
		return _grid + i*_width*_height*_depth;
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
		/* in cases where x is place exactly at the boundary floor does not evaluate to "dimSize" - 2
		 * --> MIN (...)*/
		idxX = MIN(floor((x - _pos.x)*_dVox_inv), _width - 2);
		idxY = MIN(floor((y - _pos.y)*_dVox_inv), _height - 2);
		idxZ = MIN(floor((z - _pos.z)*_dVox_inv), _depth -2);
	}

	/*
	 ** @brief: check out of bounds indexing
	 */
	__host__ bool outOfBounds_byIndex(const int &idxX, const int &idxY, const int &idxZ) const {
		return ((idxX < 0 || (idxX >= static_cast<int>(_width)-1)) || (idxY < 0 || (idxY >= static_cast<int>(_height)-1)) || (idxZ < 0 || (idxZ >= static_cast<int>(_depth)-1)));
	}
	
	__host__ bool outOfBounds_byPos(const float &x, const float &y, const float &z) const {
		return (( (x < minDim().x) || (x > maxDim().x) ) ||
				( (y < minDim().y) || (y > maxDim().y) ) ||
				( (z < minDim().z) || (z > maxDim().z) ) );
	}

	__host__ bool notOutOfBounds_byPos(const float &x, const float &y, const float &z) const {
		return (( (x >= minDim().x) && (x <= maxDim().x) ) &&
				( (y >= minDim().y) && (y <= maxDim().y) ) &&
				( (z >= minDim().z) && (z <= maxDim().z) ) );
	}

	/*
	 ** @brief: returns the voxelOct for host interpolation given a set of coordinates
	 ** This method depends highly on the implementation of the Grid. Therefore this
	 ** method is part of the Grid.
	 */
	__host__ void host_getVoxelByIndex(const int &idxX, const int &idxY, const int &idxZ, const uint &type, VoxelOctet &voxelOct) const {

		assert(idxX >= 0 && idxX < static_cast<int>(_width)-1);
		assert(idxY >= 0 && idxY < static_cast<int>(_height)-1);
		assert(idxZ >= 0 && idxZ < static_cast<int>(_depth)-1);

		// compute absolute position of vertices
		voxelOct.min.x = idxX * _dVox + _pos.x;
		voxelOct.min.y = idxY * _dVox + _pos.y;
		voxelOct.min.z = idxZ * _dVox + _pos.z;
		voxelOct.max.x = voxelOct.min.x + _dVox;
		voxelOct.max.y = voxelOct.min.y + _dVox;
		voxelOct.max.z = voxelOct.min.z + _dVox;

		// fetch data from the grid
		assert(type < 32);
		uint idx = type*_depth*_height*_width;

		idx += _width*(idxZ*_height + idxY) + idxX;
		voxelOct.data[0][0][0] = _grid[idx];
		voxelOct.data[1][0][0] = _grid[idx + 1];
		voxelOct.data[0][1][0] = _grid[idx + _width];
		voxelOct.data[1][1][0] = _grid[idx + _width + 1];

		idx += _width*_height;

#ifndef NDEBUG
		if (idx >= 32*width()*height()*depth()) {
			std::cout << type << " " << idxX << " " << idxY << " "  << idxZ  << std::endl;
			std::cout << type << " " << width()<< " " << height()<< " "  << depth()  << std::endl;
		}
#endif
		assert(idx < 32*width()*height()*depth());
		voxelOct.data[0][0][1] = _grid[idx];
		voxelOct.data[1][0][1] = _grid[idx + 1];
		voxelOct.data[0][1][1] = _grid[idx + _width];
		voxelOct.data[1][1][1] = _grid[idx + _width + 1];

	}

	/*
	 ** @brief: returns the size of the grid in terms of memory consumption
	 */
	double mySize(asUtils::sizeType_t type = asUtils::byte) const;

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

	/****************************
	 * private member variables
	 ****************************/

	float4* _grid; 	/** data of host grids */

	unsigned _numGrids; 	/** number of grids */

	float3 _maxDim;			/** max. grid dimension */
	float _voxelVol;		/** voxel volume */
	float _voxelVol_inv; 	/** pre-computed inverse of the volume */
	float _dVox_inv; 		/** pre-computed inverse of the voxel distance */
};

} // namespace

#endif /* INTRPLGRID_H_ */
