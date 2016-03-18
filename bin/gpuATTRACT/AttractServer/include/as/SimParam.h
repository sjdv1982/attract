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

#ifndef SIMPARAM_H_
#define SIMPARAM_H_

#include "config.h"

namespace as {

class SimParam {
public:
	dielec_t dielec = variable;		/** type of dielectric constant */
	float epsilon = 15;				/** dielectric constant */
	float ffelec = FELEC/epsilon;	/** precomputed factor felec/epsilon */
//	bool  useSwi;					/** using switching potential */
//	float swiOn;					/** min. switching potential distance */
//	float swiOff;					/** max. switching potential distance */
	bool useRecGrad = false;		/** using Receptor gradients */
	bool usePot = true;				/** use Potential grid */
};

}  // namespace AServ





#endif /* SIMPARAM_H_ */
