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

#ifndef TYPEMAP_H_
#define TYPEMAP_H_

#include <map>
#include <vector>
#include <initializer_list>

namespace as {

class TypeMap {
public:
	using valueType = unsigned;
	using keyType = unsigned;

	TypeMap() = default;
	TypeMap(std::initializer_list<std::pair<keyType const, valueType>> initList);

	keyType getValue(valueType key) const {
		return _map.at(key);
	}

	void setKeyValuePair(keyType key, valueType value) {
		_map[key] = value;
	}

	const static TypeMap defaultTypeMap;

private:

	std::map<valueType, keyType> _map;
};

inline TypeMap createTypeMapFromVector(const std::vector<TypeMap::keyType>& vec) {
	TypeMap map;

	for (unsigned i = 0; i < vec.size(); ++i) {
		map.setKeyValuePair(vec[i], (vec[i] == 99 ? 0: i+1)); /* valid type begins at 1 */
	}
	map.setKeyValuePair(0, 0);
	return map;
}

void applyMapping(const as::TypeMap& map, unsigned numAtoms, as::TypeMap::keyType const * typesIn, as::TypeMap::keyType* typesOut);

void applyDefaultMapping(unsigned numAtoms, as::TypeMap::keyType const * typesIn, as::TypeMap::keyType* typesOut);


}



#endif /* TYPEMAP_H_ */
