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

#include <iostream>
#include <iomanip>

#include "as/Protein.h"

/* Constructor */
as::Protein::Protein() :
	_tag(),
	_numAtoms(0), _pivot(0,0,0),
	_pos(nullptr),
	_type(nullptr), _charge(nullptr),
	_numModes(0), _modes(nullptr) {};


/* TODO: Deprecated!!! */
as::Protein::Protein(ProteinDesc desc) :
	_tag(desc.id),
	_numAtoms(desc.numAtoms),
	_pivot(0,0,0),
	_type(desc.type),
	_mappedTypes(_numAtoms),
	_charge(desc.charge),
	_numModes(desc.numModes)

{
	if (desc.pos != NULL) {
		_pos = desc.pos;
	} else {
		std::cerr << "Error: Protein description: invalid position pointer"
				<< std::endl;
		exit(EXIT_FAILURE);
	}

	if (_numModes > 0 ) {
		if (desc.modes == nullptr) {
			std::cerr << "Error: Protein description: invalid mode pointer"
					<< std::endl;
			exit(EXIT_FAILURE);
		}
		_modes = desc.modes;
	} else { // _xModes is going to be explicitly deleted!
		_modes = nullptr;
	}

	for (unsigned i = 0; i < _numAtoms; ++i) {
		_type[i] = _type[i] == 0 ? 31 : _type[i]-1;
	}

}

/* Destructor */
as::Protein::~Protein() {
	delete[] _pos;
	delete[] _charge;
	delete[] _type;
	delete[] _modes;
}


/****************************
 * public member functions
 ****************************/

float* as::Protein::getOrCreatePosPtr() {
	if (_pos == nullptr) {
		if (_numAtoms == 0) {
			std::cerr << "Error: getOrCreatePosPtr(): the number of atoms must be set before." << std::endl;
			exit(EXIT_FAILURE);
		}
		_pos = new float[3*_numAtoms];
	}
	return _pos;
}

unsigned* as::Protein::getOrCreateTypePtr() {
	if (_type == nullptr) {
		if (_numAtoms == 0) {
			std::cerr << "Error: getOrCreateTypePtr(): the number of atoms must be set before" << std::endl;
			exit(EXIT_FAILURE);
		}
		_type = new unsigned[_numAtoms];
	}
	return _type;
}

float* as::Protein::getOrCreateChargePtr() {
	if (_charge == nullptr) {
		if (_numAtoms == 0) {
			std::cerr << "Error: getOrCreateChargePtr(): the number of atoms must be set before" << std::endl;
			exit(EXIT_FAILURE);
		}
		_charge = new float[_numAtoms];
	}
	return _charge;
}

float* as::Protein::getOrCreateModePtr() {
	if (_modes == nullptr) {
		if (_numAtoms == 0) {
			std::cerr << "Error: getOrCreateModePtr(): the number of atoms must be set before" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (_numModes == 0) {
			std::cerr << "Error: getOrCreateModePtr(): the number of modes must be set before" << std::endl;
			exit(EXIT_FAILURE);
		}
		_modes = new float[3*_numAtoms*_numModes];
	}
	return _modes;
}

void as::Protein::pivotize(asUtils::Vec3f pivot) {
	/* undo preceding pivotizing if necessary */
	if (!(_pivot == asUtils::Vec3f(0,0,0))) {
		for (unsigned i = 0; i < _numAtoms; ++i) {
			xPos()[i] += _pivot[0];
			yPos()[i] += _pivot[1];
			zPos()[i] += _pivot[2];
		}
	}
	_pivot = pivot;
	if (!(_pivot == asUtils::Vec3f(0,0,0))) {
		for (unsigned i = 0; i < _numAtoms; ++i) {
			xPos()[i] -= _pivot[0];
			yPos()[i] -= _pivot[1];
			zPos()[i] -= _pivot[2];
		}
	}
}

void as::Protein::auto_pivotize() {
	/* undo preceding pivotizing if necessary */
	if (!(_pivot == asUtils::Vec3f(0,0,0))) {
		for (unsigned i = 0; i < _numAtoms; ++i) {
			xPos()[i] += _pivot[0];
			yPos()[i] += _pivot[1];
			zPos()[i] += _pivot[2];
		}
		_pivot = asUtils::Vec3f(0,0,0);
	}
	/* calculate pivot by center of mass coordinates */
	asUtils::Vec3f pivot(0,0,0);
	for (unsigned i = 0; i < _numAtoms; ++i) {
		pivot[0] += xPos()[i];
		pivot[1] += yPos()[i];
		pivot[2] += zPos()[i];
	}
	pivot /= static_cast<double>(_numAtoms);
	pivotize(pivot);
}

void as::Protein::print() {
	using namespace std;
	int precisionSetting = cout.precision( );
	ios::fmtflags flagSettings = cout.flags();
	cout.setf(ios::dec | ios::showpoint | ios::showpos);
	cout.precision(6);

	int w = 13;
//	outStream 	<< setw(w) << "DOF"
//				<< setw(w) << dof.pos.x << setw(w) << dof.pos.y << setw(w) << dof.pos.z
//				<< setw(w) << dof.ang.x << setw(w) << dof.ang.y << setw(w) << dof.ang.z;

	cout << setw(5) << "#" << setw(w) << "X" << setw(w) << "Y" << setw(w) << "Z" << setw(6) << "TYPE" << setw(w) << "CHARGE" << endl;
	for(unsigned i = 0; i < _numAtoms; ++i) {
		cout << setw(5) << i+1
			 << setw(w) << xPos()[i]
		     << setw(w) << yPos()[i]
		     << setw(w) << zPos()[i]
		     << setw(6) << type()[i]
		     << setw(6) << mappedTypes()[i]
			 << setw(w) << charge()[i]
			 << endl;
	}

	cout.precision(precisionSetting);
	cout.flags(flagSettings);
}






