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

#include "as/TypeMap.h"
#include "defaultTypeInitializerList.h"

using as::TypeMap;
using std::initializer_list;
using std::pair;

const as::TypeMap as::TypeMap::defaultTypeMap = DEFAULT_TYPEMAP_INITIALIZERLIST;

as::TypeMap::TypeMap(std::initializer_list<std::pair<keyType const, valueType>> initList):
	_map(initList) {}

void as::applyMapping(const TypeMap& map, unsigned nAtoms, TypeMap::keyType const * typesIn, TypeMap::keyType* typesOut) {
	for (unsigned i = 0; i < nAtoms; ++i) {
		typesOut[i] = map.getValue(typesIn[i]);
	}
}

void as::applyDefaultMapping(unsigned nAtoms, TypeMap::keyType const * typesIn, TypeMap::keyType* typesOut) {
	applyMapping(TypeMap::defaultTypeMap, nAtoms, typesIn, typesOut);
}
