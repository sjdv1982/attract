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
#include <cassert>

#include "SolverFactoryImpl.h"
#include "VA13Solver.h"
#include "BFGSSolver.h"
#include "LBFGS_B_Solver.h"




std::map<std::string, ema::SolverFactoryImpl::SolverType> ema::SolverFactoryImpl::_solverTypeMap =
	{
			{"VA13", SolverType::VA13},
			{"BFGS", SolverType::BFGS},
			{"LBFGS-B", SolverType::LBFGS_B}
	};

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

std::unique_ptr<ema::SolverBase> ema::SolverFactoryImpl::createSolverByNameImpl(const std::string& name) {
	using namespace std;
	/* assert that map contains element of name */
	assert(_solverTypeMap.find(name) != _solverTypeMap.end());
	switch (_solverTypeMap[name] ) {
	case VA13:
		return make_unique<VA13Solver>();
		break;
	case BFGS:
		return make_unique<BFGSSolver>();
		break;
	case LBFGS_B:
		return make_unique<LBFGS_B_Solver>();
		break;
	default:
		cerr << "Error: " << "Unknown solver specification." << endl;
		exit(EXIT_FAILURE);
		break;
	}
}


