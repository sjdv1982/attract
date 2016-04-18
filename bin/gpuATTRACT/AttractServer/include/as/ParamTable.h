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

#ifndef PARAMTABLE_H_
#define PARAMTABLE_H_


namespace as {

class AttrParamTable {
public:

	struct attractFFParams_t {
		float rc; // A = rc
		float ac; // B = ac
		int ipon;
		float rmin2;
		float emin;
	};

	/* Potential shape type of the ATTRACT force field */
	enum PotShape {
		_12_6,
		_8_6,
		undefined
	};

	typedef attractFFParams_t type;

	/* Constructor */
	AttrParamTable() : _paramTable(nullptr), _numTypes(0),
			_shape(undefined), _swiOn(0), _swiOff(0) {}



	/* Destructor */
	~AttrParamTable() {
		delete[] _paramTable;
	}

	/***************
	* G E T T E R
	***************/
	unsigned numTypes() const {
		return _numTypes;
	}

	/* read only access */
	const type* table() const {
		return _paramTable;
	}

	PotShape potShape() const {
		return _shape;
	}

	float swiOn() const {
		return _swiOn;
	}

	float swiOff() const {
		return _swiOff;
	}

	/***************
	* S E T T E R
	***************/

	void setParamTable(type* table) {
		/* avoid memory leak */
		if (_paramTable != nullptr) delete _paramTable;
		_paramTable = table;
	}

	void setNumTypes(unsigned numTypes) {
		_numTypes = numTypes;
	}

	void setPotShape(PotShape shape) {
		_shape = shape;
	}

	void setSwiOn(float swiOn) {
		_swiOn = swiOn;
	}

	void setSwiOff(float swiOff) {
		_swiOff = swiOff;
	}

	/****************************
	 * public member functions
	 ****************************/
	inline const type& getParams(const int& typeA, const int& typeB) const {
		return _paramTable[_numTypes*typeA + typeB];
	}

	/*
	 * Read and write access.
	 * Should be used for initialization
	 */
	type* getOrCreateTable() {
		if (_paramTable == nullptr) {
			if (_numTypes == 0) {
				std::cerr << "Error: getOrCreateTypePtr(): the number of atoms must be set before" << std::endl;
				exit(EXIT_FAILURE);
			}
			_paramTable = new AttrParamTable::type[_numTypes * _numTypes];
		}
		return _paramTable;
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

	type* _paramTable;
	unsigned _numTypes; /** number of particle/atom types */

	PotShape _shape; /** potential shape 12:6 or 8:6 or undefined */
	float _swiOn;
	float _swiOff;
};

}


#endif /* PARAMTABLE_H_ */
