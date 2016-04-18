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

#include <iostream>

#include "BFGSSolver.h"
#include <Eigen/Dense>


using std::cerr;
using std::endl;

ema::BFGSSolver::Options ema::BFGSSolver::settings;

template <typename T> inline
int sign(T val) {
	return (T(0) < val) - (val < T(0));
}

void ema::BFGSSolver::run(coro_t::caller_type& ca) {

	/* Algorithm 6.1, J. Nocedal, S. J. Wright, Numerical Optimization (2nd Edition), page 140 */
	Vector x_curr = getState();
	Vector x_next;
	const unsigned DIM = x_curr.rows();
	Matrix H = Matrix::Identity(DIM, DIM);
	ObjGrad objGrad_curr;
	ObjGrad objGrad_next;


	Vector p;
	Vector s;
	Vector y;

	unsigned iter = 0;
	double dObj = 0.0;


	OBJGRAD(x_curr, objGrad_curr);
	do {

		if (stats) {
			++statistic.numIter;
		}

		assert(std::isnan(objGrad_curr.obj) == false);
		p = -1.0 * H * objGrad_curr.grad;

		switch (search_type) {
		case WolfeRule:
			linesearch_WolfeRule(ca, x_curr, objGrad_curr, p, dObj, /*out*/x_next, /*out*/objGrad_next);
			break;
		default:
			cerr << "Error: linesearch type unspecified" << endl;
		}


		s = x_next - x_curr;
		y = objGrad_next.grad - objGrad_curr.grad;
		dObj = objGrad_next.obj - objGrad_curr.obj;


		double ys = y.dot(s);

		if (fabs(dObj) < settings.dObjTol) {
			if (stats) {
				statistic.convergence = BFGSStatistic::finitePrec;
			}
			break;
		}
		if (s.norm() < settings.dxTol ) {
			if (stats) {
				statistic.convergence = BFGSStatistic::finitePrec;
			}
			break;
		}

		if (ys <= 0.0) {
			/* damped BFGS-update taken from */
			/* Al-Baali, Mehiddin, Lucio Grandinetti, and Ornella Pisacane.
			 * "Damped techniques for the limited memory BFGS method for large-scale optimization."
			 * Journal of Optimization Theory and Applications 161.2 (2014): 688-699. */
			/* and */
			/* Procedure 18.2, J. Nocedal, S. J. Wright, Numerical Optimization (2nd Edition), page 537 */

			Matrix B = H.inverse();

			double theta;
			Vector Bs = B*s;
			double sBs = s.dot(Bs);
			if (ys >= 0.2 * sBs) {
				theta = 1.0; // skip updating
			} else {
				theta = 0.8 * sBs / (sBs - ys);
			}
			y = theta * y  +  (1.0-theta)*Bs;
			assert(!isnan(theta) || !isinf(theta));

			ys = y.dot(s);
		}

		assert(ys > 0.0);

		if (ys < settings.illTol) {
			if (stats) {
				statistic.convergence = BFGSStatistic::illconditioned;
			}
			break;
		}

		const double rho = 1.0 / ys;

		if (iter == 0) {
			H = ((y.dot(s)) / (y.dot(y)) * Matrix::Identity(DIM, DIM));
		}

		H = (Matrix::Identity(DIM, DIM) - rho*s*y.transpose()) * H * (Matrix::Identity(DIM, DIM) - rho*y*s.transpose()) + rho*s * s.transpose();

		x_curr = x_next;
		objGrad_curr = objGrad_next;
		iter++;

	} while ((objGrad_next.grad.lpNorm<Eigen::Infinity>() > settings.gradTol) && (iter < settings.maxIter));

	if (stats) {
		statistic.gradNorm = objGrad_curr.grad.lpNorm<Eigen::Infinity>();
		if (statistic.gradNorm <= settings.gradTol) {
			statistic.convergence = BFGSStatistic::gradTolerance;
		} else if (iter >= settings.maxIter) {
			statistic.convergence = BFGSStatistic::maxIter;
		}
	}

	return;
}

namespace ema {

inline double cubic_interpolate(const double& alpha_curr, const double& alpha_last,
		const ObjGrad& objGrad_curr, const ObjGrad& objGrad_last,
		const double& phi_dash_curr, const double& phi_dash_last )
{
	/* eq. 3.59 J. Nocedal, S. J. Wright, Numerical Optimization (2nd Edition), page 59 */

	double d1 = phi_dash_last + phi_dash_curr - 3.0 * (objGrad_last.obj - objGrad_curr.obj) / (alpha_last - alpha_curr);
	double d2 = sign(alpha_curr - alpha_last) * sqrt(d1*d1 - phi_dash_last*phi_dash_curr);
	double alpha_next = alpha_curr - (alpha_curr - alpha_last)*( (phi_dash_curr + d2 - d1)/(phi_dash_curr - phi_dash_last + 2*d2) );
	return alpha_next;
}


inline double ema::BFGSSolver::zoom(coro_t::caller_type& ca, const Vector& x0, const Vector& p,
		const ObjGrad& objGrad0, const double& phi0_dash, const double& min_alpha,
		double& alpha_lo, double& alpha_hi,
		ObjGrad& objGrad_lo, ObjGrad& objGrad_hi,
		double& phi_dash_lo, double& phi_dash_hi,
		/*OUT:*/ Vector& x_next, ObjGrad& objGrad_next)
{
	/* Algorithm 3.6, J. Nocedal, S. J. Wright, Numerical Optimization (2nd Edition), page 61 */

	ObjGrad objGrad;
	Vector x_candidate;
	double alpha = 0.0;

	unsigned count_equal_obj = 0;
	unsigned i;
	for (i = 0; i < settings.maxZoomIter; ++i) {

		alpha = cubic_interpolate(alpha_lo, alpha_hi, objGrad_lo, objGrad_hi, phi_dash_lo, phi_dash_hi);


		double alpha_min = std::min(alpha_lo, alpha_hi);
		double alpha_max = std::max(alpha_lo, alpha_hi);
		assert(alpha <= alpha_max);
		assert(alpha >= alpha_min);

		double tol = (alpha_max - alpha_min)*settings.bndTol;
		if (alpha < alpha_min + tol) {
			alpha = alpha_min + tol;
		} else if (alpha > alpha_max - tol) {
			alpha = alpha_max - tol;
		}

		x_candidate = x0 + alpha*p;

		OBJGRAD(x_candidate, objGrad);

		double phi_dash = p.dot(objGrad.grad);
		assert(alpha <= MAX(alpha_hi,alpha_lo));
		assert(alpha > MIN(alpha_hi,alpha_lo));

		if ((objGrad.obj > objGrad0.obj + settings.c1*alpha*phi0_dash) || (objGrad.obj >= objGrad_lo.obj)) {
			if (fabs(objGrad.obj - objGrad_lo.obj) < 1.0e-9) {
				if ((objGrad.grad-objGrad_lo.grad).lpNorm<Eigen::Infinity>() < 1.0e-5) {
					break;
				}
				++count_equal_obj;
				if (count_equal_obj >= settings.maxEqualObj) {
					break;
				}
			} else if (alpha < min_alpha) {
				break;
			}
			alpha_hi = alpha;
			objGrad_hi = objGrad;
			phi_dash_hi = phi_dash;
		} else {
			if (fabs(phi_dash) <= -settings.c2*phi0_dash) {
				assert(fabs(phi_dash) <= fabs(settings.c2*phi0_dash));
				x_next = x_candidate;
				objGrad_next = objGrad;
				return alpha;
			}
			if (phi_dash*(alpha_hi - alpha_lo) >= 0.0) {

				alpha_hi = alpha_lo;
				objGrad_hi = objGrad_lo;
				phi_dash_hi = phi_dash_lo;
			}
			alpha_lo = alpha;
			objGrad_lo = objGrad;
			phi_dash_lo = phi_dash;
		}

	}

	/* return current alpha in case we found no better one */
	x_next = x_candidate;
	objGrad_next = objGrad;
	return alpha;
}

} // namespace

double ema::BFGSSolver::linesearch_WolfeRule(coro_t::caller_type& ca, const Vector& x0, const ObjGrad& objGrad0, const Vector& p,
		const double& dObj,
		/*OUT:*/ Vector& x_next, ObjGrad& objGrad_next)
{
	/* Algorithm 3.5, J. Nocedal, S. J. Wright, Numerical Optimization (2nd Edition), page 60 */

	Vector x_candidate;

	ObjGrad objGrad_last = objGrad0;
	ObjGrad objGrad_curr;

	const double phi0_dash = p.dot(objGrad0.grad);
	assert(phi0_dash < 0.0);

	double phi_dash_curr;
	double phi_dash_last = phi0_dash;

	/* set the initial step-length of the line search */
	double alpha_curr = 1.0;
	double alpha_last = 0.0;

	const double p_norm = p.norm();
	const double p_max = p.lpNorm<Eigen::Infinity>();

	if (dObj <= 0.0) {
		alpha_curr =  std::min(alpha_curr, 1.0/p_max);
	} else {
		/* scale alpha if direction has large norm (taken from orig. ATTRACT minfor.f)*/
		alpha_curr = std::min(alpha_curr, (dObj+dObj)/ (-phi0_dash));
	}

	double min_alpha = 1e-5/p_max;



	for (unsigned i = 0; i < settings.maxLineSearchIter; ++i) {
		if (stats) {
			++statistic.numLineSearchIter;
		}

		x_candidate = x0 + alpha_curr*p;
		OBJGRAD(x_candidate, objGrad_curr);

		phi_dash_curr = p.dot(objGrad_curr.grad);

		if ((objGrad_curr.obj > objGrad0.obj + settings.c1*alpha_curr*phi0_dash) || (objGrad_curr.obj >= objGrad_last.obj && i > 0)) {
			return zoom(ca, x0, p, objGrad0, phi0_dash, min_alpha,
					alpha_last		, 	alpha_curr		,
					objGrad_last	, 	objGrad_curr	,
					phi_dash_last	, 	phi_dash_curr	,
					/*OUT:*/ x_next , 	objGrad_next	);
		}

		if (fabs(phi_dash_curr) <= -settings.c2*phi0_dash) {
			assert(fabs(phi_dash_curr) <= fabs(settings.c2*phi0_dash));
			x_next = x_candidate;
			objGrad_next = objGrad_curr;
			return alpha_curr;
		}

		if (phi_dash_curr >= 0.0) {
			return zoom(ca, x0, p, objGrad0, phi0_dash, min_alpha,
					alpha_curr			, 	alpha_last		,
					objGrad_curr		, 	objGrad_last	,
					phi_dash_curr		, 	phi_dash_last	,
					/*OUT:*/ x_next 	, 	objGrad_next	);
		}

		double tmp;
		double cubic_interp = cubic_interpolate(alpha_curr, alpha_last, objGrad_curr, objGrad_last, phi_dash_curr, phi_dash_last);

		if (cubic_interp <= settings.alpha_max/p_norm && cubic_interp > alpha_curr ) {
			tmp = cubic_interp;
		} else {
			tmp = 1.5*alpha_curr;
		}

		alpha_last = alpha_curr;
		objGrad_last = objGrad_curr;
		phi_dash_last = phi_dash_curr;

		alpha_curr = tmp;
	}

	x_next = x_candidate;
	objGrad_next = objGrad_curr;
	return alpha_curr;
}
