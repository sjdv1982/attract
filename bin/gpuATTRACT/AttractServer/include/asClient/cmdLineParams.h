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

#ifndef PARAMS_H_
#define PARAMS_H_

#include <string>
#include <cstring>
#include <sstream>

namespace asClient {

// parameter processing
template<typename T>
bool getParam(std::string param, T &var, int argc, char **argv)
{
    const char *c_param = param.c_str();
    for(int i=argc-1; i>=1; i--)
    {
        if (argv[i][0]!='-') continue;
//        cout << argv[i]+1 << endl;
        if (strcmp(argv[i]+1, c_param)==0)
        {
            if (!(i+1<argc))
            	return true;
            std::stringstream ss;
            ss << argv[i+1];
            ss >> var;
            return true;
        }

    }
    return false;
}

}


#endif /* PARAMS_H_ */
