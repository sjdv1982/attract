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

#ifndef ASDB_READ_FILE_H_
#define ASDB_READ_FILE_H_

#include <string>
#include <vector>

#include "as/Protein.h"
#include "as/GridUnion.h"
#include "as/ParamTable.h"
#include "as/asTypes.h"
#include "asUtils/Vec3.h"

namespace asDB {

/*
 ** @brief: Creates Protein object, reads pdb and assignes values.
 ** Supports old and new recduced format. Memory management needs to
 ** handled outside (e.g. by using a shared_ptr or manual deletion).
 */
as::Protein* createProteinFromPDB (std::string filename);

/*
 ** @brief: as above but user provides Protein pointer.
 */
void readProteinFromPDB(as::Protein*, std::string filename);

/*
 ** @brief: read the number of atoms of a protein from a pdb-file.
 */
unsigned readProteinSizeFromPDB(std::string filename);


/*
 ** @brief: Creates GridUnion from ATTRACT grid file (original
 ** format). Memory management needs to handled outside
 ** (e.g. by using a shared_ptr or manual deletion)
 */
as::GridUnion* createGridFromGridFile(std::string filename);

/*
 ** @brief: as above but user provides GridUnion pointer.
 */
void readGridFromGridFile(as::GridUnion*, std::string filename);

/*
 ** @brief: deprecated
 ** Reads binary file according to GridDesc.
 ** Memory management needs to handled outside
 ** (e.g. by using a shared_ptr or manual deletion).
 */
as::GridUnion* createGridUnionFromDumpFile (std::string filename);

/*
 ** @brief: reads an ATTRACT forcefield parameter file and creates an
 ** ATTRACT parameter table object. Memory management needs to handled outside
 ** (e.g. by using a shared_ptr or manual deletion)
 */
as::AttrParamTable* createParamTableFromFile(std::string filename);

/*
 ** @brief: as above but user provides AttrParamTable pointer.
 */
void readParamTableFromFile(as::AttrParamTable* table, std::string filename);


/*
 ** @brief: reads a (ATTRACT-) .dat containing DOFs. The number of elements read
 ** is returned by the parameter numEl. Normal modes are not yet supported.
 ** Memory needs to deleted by the caller. (e.g. by using a shared_ptr or manual deletion)
 ** ToDo: Read normal modes
 */
void readDOFFromFile(std::string filename, std::vector<std::vector<as::DOF>>& DOF_molecules);
void readDOFFromFile(std::string filename, std::vector<std::vector<as::DOF>>& DOF_molecules, const std::vector<bool> &has_ensemble);

/*
 ** @brief: reads the header of a (ATTRACT-) .dat containing DOFs.
 */
void readDOFHeader(std::string filename, std::vector<asUtils::Vec3f>& pivots,
		bool& auto_pivot, bool& centered_receptor, bool& centered_ligands);

std::vector<std::string> readFileNamesFromEnsembleList(std::string filename);

std::vector<unsigned> readGridAlphabetFromFile(std::string filename);

} // namespace asDB









#endif /* ASDB_READ_FILE_H_ */
