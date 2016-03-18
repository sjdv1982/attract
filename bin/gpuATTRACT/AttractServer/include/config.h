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

#ifndef CONFIG_H_
#define CONFIG_H_

/* DataManagement defines */

#define DEVICE_MAXGRIDS 6		/** Max. supported grid per Device */
#define DEVICE_MAXPROTEINS 20	/** Max. supported proteins per Device */

/* mngt defines */

#define NUMITEMS 100		/** number worker items provided in advance. MAX value = 10000; */

/* core defines */

#define MAXMODES 0 			/** Max. supported harmonic modes (receptor + ligand). */
#define DEVBUFSIZE 1000000 	/** Device buffer size = max number of atoms processed simultaneously */
#define HOSTBUFSIZE 2000 	/** Device buffer size = max number of atoms processed simultaneously */
#define MINPROTSIZE 500		/** Min. number of atoms per protein. To determine the DOF buffer size */
#define ITEMSIZE 1000		/** Size of worker item = max. number of DOFs that are simultaneously processed on the device */

#define BLSZ_TRAFO 128
#define BLSZ_REDUCE 512
#define BLSZ_INTRPL 128
#define WARP_SIZE 32		/** Device warp size = number of threads per warp */

#define USEMODES 0 /* Use normal modes */

/*
 * The felec constant is the electrostatic energy, in kcal/mol,
 * between two electrons, at 1 A distance, in vacuum
 * The formula is:
 *  e**2 * NA * KC * Ang * kcal
 * where:
 *  e = charge of the electron, in Coulomb
 *  NA = Avogadro's number
 *  KC = Coulomb force constant
 *  Ang = the size of an Angstrom (10**-10 meter)
 *  kcal = the amount of Joules per kcal, 4184
 */
#define FELEC 332.053986

//#define KERNELTIMING


#endif /* CONFIG_H_ */
