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

#include "asUtils/Thread.h"

/* Constructor */
asUtils::Thread::Thread() {}

/* Destructor */
asUtils::Thread::~Thread() {}


/****************************
 * public member functions
 ****************************/
std::thread::id asUtils::Thread::id() {
	return _thread.get_id();
}
void asUtils::Thread::start() {
	_thread = std::thread(&Thread::run, this);
}
void asUtils::Thread::join() {
	_thread.join();
}
void asUtils::Thread::detach() {
	_thread.detach();
}
void asUtils::Thread::joinable() {
	_thread.joinable();
}

/****************************
 * protected member functions
 ****************************/

/****************************
 * private member functions
 ****************************/


