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

#ifndef SOLVERFACTORY_H_
#define SOLVERFACTORY_H_

#include <string>
#include <memory>

#include "SolverBase.h"

namespace ema {

class SolverFactory {
public:
	virtual ~SolverFactory(){};

	std::unique_ptr<SolverBase> createSolverByName (const std::string& name) {
		return createSolverByNameImpl(name);
	}

private:
	virtual std::unique_ptr<SolverBase> createSolverByNameImpl(const std::string& name) = 0;
};

}  // namespace ema



#endif /* SOLVERFACTORY_H_ */
