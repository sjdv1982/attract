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

#ifndef BFGSSOLVER_H_
#define BFGSSOLVER_H_

#include "SolverBase.h"

#include <memory>

namespace ema {

struct BFGSStatistic : public Statistic {

	unsigned numIter = 0;
	unsigned numLineSearchIter = 0;
	double gradNorm = 0.0;

	enum Convergence {
		maxIter = 0,
		gradTolerance = 1,
		finitePrec = 2,
		illconditioned = 3,
		unspecified
	};
	Convergence convergence = unspecified;

	virtual Statistic* getCopy() const override {
		return static_cast<Statistic*> (new BFGSStatistic(*this));
	}

	virtual std::ostream& print(std::ostream& stream) const override {
		using namespace std;
		int precisionSetting = stream.precision( );
		ios::fmtflags flagSettings = stream.flags();
		stream.setf(ios::scientific);
		stream.precision(5);

		stream << numRequests << "\t" << numIter << "\t" << numLineSearchIter << "\t" << gradNorm << "\t" << convergence;

		stream.precision(precisionSetting);
		stream.flags(flagSettings);
		return stream;
	}
};


class BFGSSolver : public SolverBase {
public:
	BFGSSolver() : SolverBase(), search_type(WolfeRule) {}
	virtual ~BFGSSolver() {};

	BFGSSolver(const BFGSSolver& ) = delete;
	BFGSSolver& operator= (const BFGSSolver& ) = delete;

	BFGSSolver(BFGSSolver &&) = default;
	BFGSSolver& operator= (BFGSSolver&& ) = default;

	/* types */
	typedef enum linesearch {
		WolfeRule
	} linesearch_t;

	std::unique_ptr<Statistic> getStats() const override {
		return std::unique_ptr<Statistic> (statistic.getCopy());
	}

	struct Options {
		/* Solver Options */
		double gradTol = 1e-5; 				/** gradient tolerance: exit if gradient are < gradTol */
		double dObjTol = 0.5e-5;			/** delta objective tol: exit if difference of subsequent
												objective evaluations is < dObjTol */
		double dxTol = 1.0e-14;				/** spatial tolerance: exit if the norm spatial variation is < dxTol */
		double illTol = 1e-16;				/** ill-conditionedness tolerance: exit if y.dot(s) < illTol */
		unsigned maxIter = 100;				/** max solver iterations */

		/* Wolfe Conditions/Line Search */
		double c1 = 1.0e-4;					/** constant for sufficient decrease condition */
		double c2 = 0.9;					/** constant for curvature condition */
		/* Bracketing */
		unsigned maxLineSearchIter = 100; 	/** max line search iterations */
		double alpha_max = 20; 				/** max. allowed alpha. not used at the moment due to difficulties */
		/* Selection/Zoom */
		unsigned maxZoomIter = 100;			/** max. zoom iterations */
		unsigned maxEqualObj = 4;			/** max. allowed equal objective evaluations */
		double bndTol = 0.1;				/** boundary tolerace in search interval */
	};

	static void setOptions(Options opt) {settings = opt;}


private:

	void run(coro_t::caller_type& ca) override;

	double linesearch_WolfeRule(coro_t::caller_type& ca, const Vector& x0, const ObjGrad& objGrad0, const Vector& p,
			const double& dObj,
			/*OUT:*/ Vector& x_next, ObjGrad& objGrad_next);

	double zoom(coro_t::caller_type& ca, const Vector& x0, const Vector& p,
		const ObjGrad& objGrad0, const double& phi0_dash, const double& min_alpha,
		double& alpha_lo, double& alpha_hi,
		ObjGrad& objGrad_lo, ObjGrad& objGrad_hi,
		double& phi_dash_lo, double& phi_dash_hi,
		/*OUT:*/ Vector& x_next, ObjGrad& objGrad_next);

	/* solver options */
	static Options settings;

	linesearch_t search_type;

	/* Statistics */


	BFGSStatistic statistic;

	virtual Statistic* internal_getStats() override {
		return static_cast<Statistic*>(&statistic);
	}

};

} //namespace

#endif /* BFGSSOLVER_H_ */
