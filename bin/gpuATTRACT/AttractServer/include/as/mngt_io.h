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

#ifndef MNGT_IO_H_
#define MNGT_IO_H_

#include <string>

#include "as/io_config.h"
#include "as/DataManagement.h"

namespace as {

class mngt_io {
public:
	/* Constructor */
	// private

	/* Destructor */

	/***************
	* G E T T E R
	***************/

	/***************
	* S E T T E R
	***************/

	/****************************
	 * public member functions
	 ****************************/
	static std::string heading();
	static std::string clientInfo(DataManagement* mngt);
	static std::string dataInfo(const DataManagement* mngt);
	static std::string deviceInfo(const DataManagement* mngt);


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
	mngt_io() {}
	/****************************
	 * private member functions
	 ****************************/

	/****************************
	 * private member variables
	 ****************************/
	static const char fill  = '#';
	static const char vsepa = '|';
	static char const hsepa = '_';
	static const unsigned colWidth = COLWIDTH;
};

}
#endif /* MNGT_IO_H_ */
