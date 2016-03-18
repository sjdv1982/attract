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

#ifndef OCCUPANCYLAYOUT_H_
#define OCCUPANCYLAYOUT_H_

#include <vector>
#include <unordered_map>
#include <cassert>
#include <iostream>

namespace as {

/*
 ** @brief: Book keeping class for populating an external array with
 ** arbitrary objects identified by a unique object-ID.
 ** It keeps track which value is stored at which array location or,
 ** vice versa, which array location contains which value.
 ** An unordered map is used for fast position lock up according to a value.
 */
class OccupancyLayout {
public:
	/* Constructor */
	OccupancyLayout(int size) :
		_layout(size, -1), _locationLookUp(size), _capacity(size) {};

	/* Destructor */
	~OccupancyLayout() {};

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
	 ** @brief: returns the number of valid objects
	 */
	int size() {
		return _locationLookUp.size();
	}

	int getFirstEmptyLocation () const {
		for (int i = 0; i < static_cast<int>(_layout.size()); ++i) {
			if (_layout[i] == -1) {
				return i;
			}
		}
		return -1;
	}

	void addObject (int objID, int location) {
		assert(_locationLookUp.find(objID) == _locationLookUp.end()); // assert that objID(key) is not yet contained
		_layout[location] = objID;
		_locationLookUp[objID] = location;
	}

	int removeObj (int objID) {
		assert(_locationLookUp.find(objID) != _locationLookUp.end()); // assert that objID(key) is contained
		int location = _locationLookUp[objID];
		_locationLookUp.erase(objID);
		_layout[location] = -1;
		return location;
	}

	int removeObjAtLocation (int location) {
		int objID = _layout[location];
		_layout[location] = -1;
		_locationLookUp.erase(objID);
		return objID;
	}

	inline int getLocation (const int &objID) const {
		assert(_locationLookUp.find(objID) != _locationLookUp.end()); // assert that objID(key) is contained
		// Don't want to add a new element that does not exist ==> use .at() instead of [] operator
		return _locationLookUp.at(objID);
//		return _locationLookUp[objID];
	}

	inline int getObject (const int &location) const {
		assert(_layout[location] > -1);
		return _layout[location];
	}

	inline bool isfull() {
		return static_cast<int>(_locationLookUp.size()) == _capacity;
	}

	inline bool contains(int objID) const {
		return _locationLookUp.find(objID) != _locationLookUp.end();
	}

	inline std::set<int> locationSet() {
		std::set<int> locSet;
		for (auto kv : _locationLookUp) {
			locSet.insert(kv.second); // values
		}
		return locSet;
	}

	inline std::set<int> objectSet() {
		std::set<int> objSet;
		for (auto kv : _locationLookUp) {
			objSet.insert(kv.first); // keys
		}
		return objSet;
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

	std::vector<int> _layout;
	std::unordered_map<int, int> _locationLookUp;

	int _capacity;
};

} // namespace


#endif /* OCCUPANCYLAYOUT_H_ */
