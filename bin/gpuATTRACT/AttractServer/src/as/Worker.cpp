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

#include "as/Worker.h"
#include "as/ServerManagement.h"

/* Constructor */


/* Destructor */



/****************************
 * public member functions
 ****************************/

/****************************
 * protected member functions
 ****************************/
as::Protein* as::Worker::getProtein(const unsigned& globId) {
	return _S_mngt.DataMngt()->getProtein(globId);
}

as::GridUnion* as::Worker::getGridUnion(const unsigned& globId) {
	return _S_mngt.DataMngt()->getGridUnion(globId);
}

as::AttrParamTable* as::Worker::getParamTable() {
	return _S_mngt.DataMngt()->getParamTable();
}

as::SimParam* as::Worker::getSimParam() {
	return _S_mngt.DataMngt()->getSimParam();
}

/****************************
 * private member functions
 ****************************/



