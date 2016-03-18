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

#ifndef LBFGSSOLVER_H_
#define LBFGSSOLVER_H_

#include "SolverBase.h"

namespace ema {

struct LBFGS_B_Statistic : public Statistic {

	virtual Statistic* getCopy() const override {
		return static_cast<Statistic*> (new LBFGS_B_Statistic(*this));
	}

	virtual std::ostream& print(std::ostream& stream) const override {
		using namespace std;
		int precisionSetting = stream.precision( );
		ios::fmtflags flagSettings = stream.flags();
		stream.setf(ios::scientific);
		stream.precision(5);

		stream << numRequests << endl;

		stream.precision(precisionSetting);
		stream.flags(flagSettings);
		return stream;
	}
};

class LBFGS_B_Solver : public SolverBase {
public:
	LBFGS_B_Solver() : SolverBase() {}
	virtual ~LBFGS_B_Solver() {};

	LBFGS_B_Solver(const LBFGS_B_Solver& ) = delete;
	LBFGS_B_Solver& operator= (const LBFGS_B_Solver& ) = delete;

	LBFGS_B_Solver(LBFGS_B_Solver &&) = default;
	LBFGS_B_Solver& operator= (LBFGS_B_Solver&& ) = default;

	std::unique_ptr<Statistic> getStats() const override {
		return std::unique_ptr<Statistic> (statistic.getCopy());
	}

	struct Options {
		/* Solver Options */
		unsigned maxFunEval = 500;
	};
	static void setOptions(Options opt) {settings = opt;}


private:

	void run(coro_t::caller_type& ca) override;

	/* solver options */
	static Options settings;

	/* Statistics */
	LBFGS_B_Statistic statistic;

	virtual Statistic* internal_getStats() override {
		return static_cast<Statistic*>(&statistic);
	}

};

}



#endif /* LBFGSSOLVER_H_ */
