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

#include <sstream>
#include <iomanip>

#include "as/mngt_io.h"
#include "as/io_helpers.h"


std::string as::mngt_io::heading() {
	using namespace std;
	stringstream stream;

	string heading = " Data Info ";
	stream << setfill(fill) << setw(60) << centered(heading) << endl;

	return stream.str();

}

std::string as::mngt_io::clientInfo(DataManagement* mngt) {
	using namespace std;
	stringstream stream;

	string heading = "Clients";
	stream << "# "<< heading << " #" << endl;


	unsigned numCol = 6;
	unsigned width = numCol * colWidth;

	stream << setfill(hsepa) << setw(width) << "" << endl;
	stream << setfill(' ');

	stream	<< setw(2*colWidth - 1 )<< "" << vsepa << setw(4*colWidth-1) << centered("IDs") <<vsepa << endl;

	stream 	<< setw(colWidth - 1) << centered("clientID") << vsepa
			<< setw(colWidth - 1) << centered("name") << vsepa
			<< setw(colWidth - 1) << centered("loc") << vsepa
			<< setw(colWidth - 1) << centered("glob") << vsepa
			<< setw(colWidth - 1) << centered("dev") << vsepa
			<< setw(colWidth - 1) << centered("devloc") << vsepa << endl;
	stream << setfill(hsepa) << setw(width) << "" << endl;
	stream << setfill(' ');

	unsigned numClients = mngt->_clientData.size();
//		unsigned numValidClients = 0;
//
//		/* Collect Data */
////		vector<unsigned> clID;
////		vector<vector<int>> clientGrids;
////		vector<vector<int>>	clientProtein
//

//	cout << mngt->_sharedProteins.size() << endl;
//	cout << mngt->_sharedProteins[0].get()->tag() << endl;
	for (unsigned i = 0; i < numClients; ++i) {
		if (!mngt->_clientData[i].isValid()) {continue;}
		stream << setw(1*colWidth - 1 ) << i << vsepa << setw(1*colWidth - 1 ) << centered("Proteins") << vsepa << endl;
		for (unsigned j = 0; j < mngt->_clientData[i].proteins().size(); ++j) {
			int globId = mngt->_clientData[i].globalProteinId(j);
			stream 	<< setw(colWidth - 0) << vsepa
					<< setw(colWidth - 1) << mngt->_sharedProteins[globId].get()->tag() << vsepa
					<< setw(colWidth - 1) << j << vsepa
					<< setw(colWidth - 1) << globId << vsepa;
			stringstream dev, devloc;
			if (globId >= static_cast<int>(mngt->_prot_deviceIDs.size())) {
				dev << "n.a.";
				devloc << "n.a.";
			} else {

				const std::set<unsigned> &devices = mngt->_prot_deviceIDs[globId];
				if (devices.size() == 0) {
					dev << "n.a.";
					devloc << "n.a.";
				} else {
					for (std::set<unsigned>::iterator it = devices.begin(); it != devices.end(); ++it) {
						dev << " " << *it;
						devloc << " " << mngt->_prot_deviceOCC[*it].getLocation(globId);
					}
				}
			}

			stream 	<< setw(colWidth - 1) << dev.str() << vsepa
					<< setw(colWidth - 1) << devloc.str() << vsepa
					<< endl;
		}
		stream << setw(1*colWidth - 1 ) << i << vsepa << setw(1*colWidth - 1 ) << centered("Grids") << vsepa << endl;
		for (unsigned j = 0; j < mngt->_clientData[i].grids().size(); ++j) {
			int globId = mngt->_clientData[i].globalGridId(j);
			stream 	<< setw(colWidth - 0) << vsepa
					<< setw(colWidth - 1) << mngt->_sharedGridUnions[globId].get()->tag() << vsepa
					<< setw(colWidth - 1) << j << vsepa
					<< setw(colWidth - 1) << globId << vsepa;
			stringstream dev, devloc;
			if (globId >= static_cast<int>(mngt->_grid_deviceIDs.size())) {
				dev << "n.a.";
				devloc << "n.a.";
			} else {

				const std::set<unsigned> &devices = mngt->_grid_deviceIDs[globId];
				if (devices.size() == 0) {
					dev << "n.a.";
					devloc << "n.a.";
				} else {
					for (std::set<unsigned>::iterator it = devices.begin(); it != devices.end(); ++it) {
						dev << " " << *it;
						devloc << " " << mngt->_grid_deviceOCC[*it].getLocation(globId);
					}
				}
			}

			stream 	<< setw(colWidth - 1) << dev.str() << vsepa
					<< setw(colWidth - 1) << devloc.str() << vsepa
					<< endl;


		}
		stream << setfill(hsepa) << setw(width) << "" << endl;
		stream << setfill(' ');
	}

	return stream.str();
}
//std::string AServ::mngt_io::dataInfo(DataManagement* mngt){
//
//}
//std::string AServ::mngt_io::deviceInfo(DataManagement* mngt){
//
//}


