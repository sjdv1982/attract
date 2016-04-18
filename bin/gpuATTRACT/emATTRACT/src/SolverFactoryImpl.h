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

#ifndef SOLVERFACTORYIMPL_H_
#define SOLVERFACTORYIMPL_H_

#include <map>

#include "SolverFactory.h"

namespace ema {

class SolverFactoryImpl : public SolverFactory {
private:
	std::unique_ptr<SolverBase> createSolverByNameImpl(const std::string& name) override;


	enum SolverType {
		VA13,
		BFGS,
		LBFGS_B,
		unspecified
	};
	static std::map<std::string, SolverType> _solverTypeMap;
};

}  // namespace ema




#endif /* SOLVERFACTORYIMPL_H_ */
