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

#ifndef DEVICEDATAFACTORY_H_
#define DEVICEDATAFACTORY_H_

#include <string>
#include "as/IntrplGrid.h"
#include "as/NLGrid.h"
#include "as/GridUnion.h"
#include "as/Protein.h"
#include "as/ParamTable.h"
#include "as/asTypes.h"

namespace as {

/*
 ** @brief: Performs necessary steps to initialize a gridUnion,
 ** a Protein or a ParamTable object for use on the GPU.
 **
 */
class DeviceDataFactory {
public:
	/***************
	* G E T T E R
	***************/

	/***************
	* S E T T E R
	***************/

	/****************************
	 * public member functions
	 ****************************/
	static cudaGridUnionDesc initDeviceGridUnion(const GridUnion* gridUnion, int deviceId);
	static cudaProteinDesc initDeviceProtein(const Protein* protein, int deviceId);
	static cudaParamTableDesc initDeviceParamTable (const AttrParamTable* table, int deviceId);

	static void disposeDeviceGridUnion(hostGridUnionResource desc, int deviceId);
	static void disposeDeviceProtein(hostProteinResource desc, int deviceId);
	static void disposeDeviceParamTable(hostParamTableResource desc, int deviceId);
	

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

	/* Constructor */
	DeviceDataFactory() {}; // object instantiation is private!

	static cudaIntrplGridDesc initIntrpl(const IntrplGrid* grid);
	static cudaNLGridDesc initNL(const NLGrid* grid);

	/****************************
	 * private member variables
	 ****************************/

};

} // namespace


#endif /* DEVICEDATAFACTORY_H__ */
