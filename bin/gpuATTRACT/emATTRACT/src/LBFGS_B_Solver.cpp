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

#include <cfloat>
#include <cstring>

#include "LBFGS_B_Solver.h"

#define lbfgsb_FUN setulb_

extern "C" void lbfgsb_FUN(int* n, int* m, double* x, double* l, double* u,
        int* nbd, double* f, double* g, double* factr, double* pgtol,
        double* wa, int* iwa, char* task, int* iprint, char* csave, int*
        lsave, int* isave, double* dsave);

ema::LBFGS_B_Solver::Options ema::LBFGS_B_Solver::settings;

/* Struct to reserve storage for working variables */
struct LBFGS_B_WorkingStruct {
	/* ToDo: some of that variables might be shared by all Solver instances.
	 * So we could put them to the Solver options */

	int 	n;
	int 	m;
	double* l;
	double* u;
	int* 	nbd;
	int 	niter;     /* number of iterations so far */
	int 	max_iter;  /* max number of iters, default = 99, also check fun <= fun_min before stopping */
	double 	fun_min;    /* record the minimal obj. fun. val. */
	double 	factr;
	double 	pgtol;
	int 	iprint;   /* see the comment in lbfgsb.f for usage of this field */
	double 	*wa;
	int* 	iwa;
	char 	task[60];
	char 	csave[60];
	int 	lsave[4];
	int 	isave[44];
	double 	dsave[29];

	/* create an LBFGS_B_WorkingStruct object
	 * n is the dimension of the variable
	 * m is the number of the corrections used in BFGS update
	 * l is the array of lower bounds
	 * u is the array of upper bounds
	 * nbd is the indicator array:
	   nbd(i)=0 if x(i) is unbounded,
	          1 if x(i) has only a lower bound,
	          2 if x(i) has both lower and upper bounds,
	          3 if x(i) has only an upper bound.
	 * if all of l, u, nbd are NULL, this is equivalent to unconstrained
	 * optimization
	 */
	LBFGS_B_WorkingStruct(int _n, int _m, double* _l, double* _u, int* _nbd):
		n(_n),
		m(_m),
		l(_l),
		u(_u),
		nbd(_nbd),
	    niter(0),
	    max_iter(99),
	    fun_min(DBL_MAX),
	    factr(1.0E+1),
	    pgtol(1.0E-5),
	    iprint(105), // -1: no output
	    wa(new double[(2*m+4)*n + 11*m*m + 8*m]),
	    iwa(new int[3*n])
	{
		for (int i = 0; i < 60; ++i)
			task[i] = ' ';
		strncpy(task, "START", 5);
	}

	~LBFGS_B_WorkingStruct() {
		delete[] wa;
		delete[] iwa;
	}

	LBFGS_B_WorkingStruct(LBFGS_B_WorkingStruct const&) = delete;
	LBFGS_B_WorkingStruct& operator= (LBFGS_B_WorkingStruct const&) = delete;
	LBFGS_B_WorkingStruct(LBFGS_B_WorkingStruct&&) = delete;
	LBFGS_B_WorkingStruct& operator= (LBFGS_B_WorkingStruct&&) = delete;

};

/* wrapper for fortran entry point */
int lbfgsb_run(LBFGS_B_WorkingStruct& opt, double* x, double* f, double* g) {

	if (*f < opt.fun_min) {
		opt.fun_min = *f;
	}

	lbfgsb_FUN(&opt.n, &opt.m, x, opt.l, opt.u, opt.nbd, f, g,
			&opt.factr, &opt.pgtol, opt.wa, opt.iwa, opt.task,
			&opt.iprint, opt.csave, opt.lsave, opt.isave, opt.dsave);

	if (opt.task[0] == 'F' && opt.task[1] == 'G') {
		opt.niter++;
		return 1;
	} else if (strncmp(opt.task, "NEW_X", 5) == 0) {
		opt.niter++;
		if (opt.niter >= opt.max_iter) {
			if (*f > opt.fun_min) {
				/* fprintf(stdout, "(lbfgsb_run) max iter reached but fun val is not minimum, wait for another iteration\n"); */
				if (opt.niter >= opt.max_iter + 10) {
				/* really too many iterations made. Stopping */
				strncpy(opt.task, "STOP,ERR: could not converge", 28);
				}
			} else {
				strncpy(opt.task, "STOP, max iterations reached", 28);
			}
		}
		return 1;
	} else if (strncmp(opt.task, "STOP", 4) == 0) {
		return 0;
	} else if (strncmp(opt.task, "CONV", 4) == 0) {
		return 0;
	} else if (strncmp(opt.task, "ABNO", 4) == 0) {
		return -1;
	} else if (strncmp(opt.task, "ERROR", 5) == 0) {
		return -1;
	} else {
		opt.task[59] = '\0';
//		fprintf(stderr, "(lbfgsb_run) unknown return value in task:%s\n", opt.task);
		return -1;
	}
}

void ema::LBFGS_B_Solver::run(coro_t::caller_type& energyAndGradients) {
	// dimension of problem
	int n = state.rows();

	double l[n];
	double u[n];
	int nbd[n];

	for (int i=0;i<n; i++)
	{
		l[i]=0;
		u[i]=0;
		nbd[i]=0;
	}  //unconstrained problem

	int m = 5;
	LBFGS_B_WorkingStruct m_opt(n, m, &l[0], &u[0], &nbd[0]);


	m_opt.iprint=-1;

	objective.obj = DBL_MAX;

	m_opt.max_iter = settings.maxFunEval;

	while (1) {
		double& f = objective.obj;
		double* x = state.data();
		double* g = objective.grad.data();
		int rc = lbfgsb_run(m_opt, &x[0], &f, &g[0]);
		if (rc == 0) {
			break;
		} else if (rc < 0) {
//			printf("lbfgsb stop with an error");
			break;
		} else if (rc == 1) {
			energyAndGradients();
		} else {
			assert(!"can not reach here");
		}
	}
	m_opt.task[59]='\0'; //add a null terminating character to task[]

//	std::cout << m_opt.task  << " |  " << m_opt.niter << " iterations\n";

}




