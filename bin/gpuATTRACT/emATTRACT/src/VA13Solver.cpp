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

#include "VA13Solver.h"

using std::cerr;
using std::endl;

ema::VA13Solver::Options ema::VA13Solver::settings;

extern "C" void minfor_(void* FortranSmuggler_ptr, int const& maxFunEval,
		double const* state);


void ema::VA13Solver::run(coro_t::caller_type& ca) {
	/* Create Smuggler */
	VA13Solver::FortranSmuggler smuggler(ca, state, objective);

	/* create and fill state array */
	double state_array[state.rows()];
	for (int i = 0; i < state.rows(); ++i) {
		state_array[i] = state(i);
	}

	minfor_(&smuggler, settings.maxFunEval, state_array);
}


// Call back function for fortran to access the class.  This is a free function
// that is not a member or friend of MyClass, so it can only access public
// member variables and functions of MyClass.  BaseClass::operator() is such a
// public member function.
extern "C" void energy_for_fortran_to_call_(void* FortranSmuggler_ptr, double state_ptr[], double* energy, double grad[])
{
   // Cast to BaseClass.  If the pointer isn't a pointer to an object
   // derived from BaseClass, then the world will end.
	ema::VA13Solver::FortranSmuggler* smuggler = static_cast<ema::VA13Solver::FortranSmuggler*>(FortranSmuggler_ptr);

	/* set the state */
	ema::Vector& state = smuggler->state_ref();
	for (int i = 0; i < state.rows(); ++i) {
		state(i) = state_ptr[i];
	}

	/* call coroutine to break execution here until energy and gradients are available */
	smuggler->call_coro();

	/* get the state */
	ema::ObjGrad& objGrad = smuggler->objective_ref();
	*energy = objGrad.obj;
	for (int i = 0; i < objGrad.grad.rows(); ++i) {
		grad[i] = objGrad.grad(i);
	}

}
