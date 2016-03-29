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
#include <cuda_runtime.h>

#include "asClient/cudaArchCheck.h"

using namespace std;


int main (int argc, char *argv[]) {
	using namespace std;
	/* Check Compute Capability of devices */
	try {
		asClient::checkComputeCapability();
	} 
	catch (std::exception& e) {
		cerr << "Error: " << e.what() << endl;
		exit(EXIT_FAILURE);
	}
	float *dummy;
        cudaMalloc(&dummy, 1000); //just to force the creation of a CUDA context
}


