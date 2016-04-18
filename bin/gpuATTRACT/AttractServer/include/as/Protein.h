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

#ifndef PROTEIN_H_
#define PROTEIN_H_

#include <string>
#include <map>
#include <algorithm>
#include <cassert>

#include "as/interface.h"
#include "asUtils/Vec3.h"
#include "as/TypeMap.h"


namespace as {

/*
 ** @brief: Describes the properties of a Protein used in the AttractServer
 */
class Protein {
public:

	/*
	 ** @brief: maps (atom-) type ids to ids defined by the type_map
	 ** ToDo: implement and design the usage:
	 ** 	- maybe move to extern file
	 ** 	- store mapped types associated with gridId for each mapping
	 ** 	- for GPU: concatenate type arrays
	 */
	typedef std::map<int,int> type_map;

	/* Constructor */
	Protein();

	Protein(ProteinDesc);

	/* Destructor */
	~Protein();

	/***************
	* G E T T E R
	***************/
	unsigned numAtoms() const{
		return _numAtoms;
	}

	unsigned numModes() const{
		return _numModes;
	}

	float* xPos() const{
		return _pos;
	}

	float* yPos() const{
		return _pos + _numAtoms;
	}

	float* zPos() const{
		return _pos + 2*_numAtoms;
	}

	float* charge() const{
		return _charge;
	}

	unsigned* type() const{
		return _type;
	}

	const unsigned* mappedType(unsigned gridId = 0) const {
		return _mappedTypes.data();
	}

	unsigned* mappedTypes(unsigned gridId = 0) {
		return _mappedTypes.data();
	}

	float* xModes() const {
		return _modes;
	}

	float* yModes() const{
		return _modes + _numModes*_numAtoms;
	}

	float* zModes() const{
		return _modes + 2*_numModes*_numAtoms;
	}

	std::string tag() const{
		return _tag;
	}

	asUtils::Vec3f pivot() {
		return _pivot;
	}

	/*
	 ** @brief: initializes and returns pointers, if numAtoms and/or numModes
	 ** is not set, the methods fail
	 */
	float* getOrCreatePosPtr();

	unsigned* getOrCreateTypePtr();

	float* getOrCreateChargePtr();

	float* getOrCreateModePtr();


	/***************
	* S E T T E R
	***************/
	void setTag(std::string tag) {
		_tag = tag;
	}

	void setNumAtoms(unsigned numAtoms) {
		_numAtoms = numAtoms;
		_mappedTypes.resize(_numAtoms);
	}

	void setNumModes(unsigned numModes) {
		_numModes = numModes;
	}

	/****************************
	 * public member functions
	 ****************************/

	/*
	 ** @brief: centers the coordinates according to _pivot.
	 ** always uses orig. pdb-coordinates.
	 */
	void pivotize(asUtils::Vec3f pivot);

	/*
	 ** @brief: calculates and sets the pivot as center or mass
	 ** and calls pivotize().
	 ** always uses orig. pdb-coordinates.
	 */
	void auto_pivotize();

	/*
	 ** @brief: prints the contents of the protein.
	 */
	void print();

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

	std::string _tag;	/** identifier: filename (default) */
	unsigned _numAtoms; /** number of atoms/particles */


	asUtils::Vec3f _pivot;	/** rotation pivot */
	float *_pos;	/** Cartesian coordinates in cm-frame (pivotized) */

	unsigned* _type; 	/** atom type */
	std::vector<unsigned> _mappedTypes; /* for receptor grid mapping */
	float* _charge;	/** charge of the atoms/particle */

	unsigned _numModes; /** number of modes */
	float* _modes; /** normal mode deformation vectors */
};

}


#endif /* PROTEIN_H_ */

