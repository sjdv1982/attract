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

#include <cassert>

#include "SolverBase.h"
#include "BFGSSolver.h"
#include "VA13Solver.h"


void ema::SolverBase::start() {
	assert(state.rows() > 0);
	coro =  new coro_t(std::bind(&SolverBase::run, this, std::placeholders::_1));


}

void ema::SolverBase::step() {
	assert(this->converged() == false);
	/* make sure that objective is already set !!! */
	(*coro)();

	if (stats && !converged()) {
		Statistic* stats = internal_getStats();
		++stats->numRequests;
	}
}

void ema::SolverBase::finalize() {
	if (coro) {
		delete coro;
		coro = nullptr;
	}
}

bool ema::SolverBase::stats = false;
