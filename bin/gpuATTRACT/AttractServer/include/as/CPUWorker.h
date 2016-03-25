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

#ifndef CPUWORKER_H_
#define CPUWORKER_H_

#include "as/Worker.h"

namespace as {

class CPUWorker : public as::Worker {
public:
	/* Constructor */
	CPUWorker(ServerManagement& S_mngt,
			unsigned atomBufferSize);

	/* Destructor */
	virtual ~CPUWorker();

	/***************
	* G E T T E R
	***************/
	int id() const {
		return _id;
	}

	/***************
	* S E T T E R
	***************/
	void set(int id) {
		_id = id;
	}

	/****************************
	 * public member functions
	 ****************************/
	void run();

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
	inline int deviceLocalProteinID(const int& globalId);


	inline int deviceLocalGridID(const int& globalId);

	/****************************
	 * private member variables
	 ****************************/
	int _id;

	unsigned _atomBufferSize;
};

} // namespace


#endif /* CPUWORKER_H_ */
