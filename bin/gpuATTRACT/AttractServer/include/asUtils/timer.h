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

#ifndef TIMER_H_
#define TIMER_H_


namespace asUtils {

inline void cudaInitTimer(cudaEvent_t start, cudaStream_t stream = 0) {
	cudaEventRecord(start, stream);
}

inline double cudaGetTimer(cudaEvent_t start, cudaEvent_t stop, cudaStream_t stream = 0) {
	cudaEventRecord(stop, stream);
	cudaEventSynchronize(stop);
	float elapsedTime;
	cudaEventElapsedTime(&elapsedTime, start, stop);
	return (double) elapsedTime / 1000.0f;
}

inline double getAndResetTimer(cudaEvent_t start, cudaEvent_t stop, cudaStream_t stream = 0) {
	cudaEventRecord(stop, stream);
	cudaEventSynchronize(stop);
	float elapsedTime;
	cudaEventElapsedTime(&elapsedTime, start, stop);
	start = stop;
	return (double) elapsedTime / 1000.0f;
}

inline double getTimerDifference(cudaEvent_t start, cudaEvent_t stop) {
	float elapsedTime;
	cudaEventElapsedTime(&elapsedTime, start, stop);
	return (double) elapsedTime;
}

inline void getTimerAndPrint(cudaEvent_t start, cudaEvent_t stop, std::string str, bool reset = false, cudaStream_t stream = 0) {
//	static double time_old = 0;
	double time_new;
	time_new = cudaGetTimer(start, stop, stream);
	static double time_first = time_new;
	if (reset) {
		time_first = time_new;
	}
	std::cerr  << str << ":\tElapsed time_new: " << time_new << " s\tFactor: " << time_first/time_new <<  std::endl;
//	time_old = time_new;
}


typedef struct timespec Timer;

inline void initTimer(Timer* timer) {
	clock_gettime(CLOCK_MONOTONIC, timer);
}

inline double getTimer(Timer* timer) {
	Timer endTimer;
	clock_gettime(CLOCK_MONOTONIC, &endTimer);
	return (endTimer.tv_sec - timer->tv_sec) + (
			endTimer.tv_nsec - timer->tv_nsec)*1e-9;
}

inline double getAndResetTimer(Timer* timer) {
	Timer endTimer;
	clock_gettime(CLOCK_MONOTONIC, &endTimer);
	double result = (endTimer.tv_sec - timer->tv_sec) +
			(endTimer.tv_nsec - timer->tv_nsec)*1e-9;
	*timer = endTimer;
	return result;
}

inline double getTimerDifference(Timer* startTimer, Timer* endTimer) {
	return (endTimer->tv_sec - startTimer->tv_sec) +
			(endTimer->tv_nsec - startTimer->tv_nsec)*1e-9;
}

inline void getTimerAndPrint(Timer* timer, std::string str) {
//	static double time_old = 0;
	double time_new;
	time_new = getTimer(timer);
	static double time_first = time_new;
	std::cerr << str << ":\tElapsed time_new: " << time_new << " s\tFactor: " << time_first/time_new <<  std::endl;
//	time_old = time_new;
}

}



#endif /* TIMER_H_ */
