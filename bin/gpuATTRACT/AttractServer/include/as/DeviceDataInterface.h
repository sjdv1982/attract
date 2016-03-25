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

#ifndef DEVICEDATAINTERFACE_H_
#define DEVICEDATAINTERFACE_H_

#include "as/asTypes.h"
#include "as/SimParam.h"

namespace as {
	void setDeviceGridUnion(const deviceGridUnionDesc &desc, unsigned deviceId, unsigned localDeviceID);
	void unsetDeviceGridUnion(unsigned deviceId, unsigned localDeviceID);
	void setDeviceProtein(const deviceProteinDesc &desc, unsigned deviceId, unsigned localDeviceID);
	void unsetDeviceProtein(unsigned deviceId, unsigned localDeviceID);
	void setDeviceParamTable(const deviceParamTableDesc& desc, unsigned deviceId);
	void unsetDeviceParamTable(unsigned deviceId);
	void setDeviceSimParam(const SimParam& desc, unsigned deviceId);
	void unsetDeviceSimParam(unsigned deviceId);
}


#endif /* DEVICEDATAINTERFACE_H_ */
