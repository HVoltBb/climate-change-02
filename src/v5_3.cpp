/* Author/maintainer: Can Zhou [eidotog@gmail.com]
 * Date: Dec 1 2021
 * Version: 0.6
 */

#include <TMB.hpp>
#include "growth.h"

// v4 added observation error to initial length measurements
// v5 removed l
// v5.2 added regional effect
// v5.3 added sinusoidal month effect
// v6 c++11 compliant

// ToDo List
// 1. Variable name clean-up (Unify random effect variable names, internal variables, doubles vs ints)
// 2. Move the ability to specify an unknown number of indices to R code
// 3. ...

template <class Type>
Type objective_function<Type>::operator()()
{
	using namespace density;
	using namespace Eigen;
	grow::g<Type> tuna;

	// Flags
	DATA_INTEGER(PAR1); // Baseline model
	DATA_INTEGER(PAR2); // Deprecated. Pending clean-up
	DATA_SCALAR(PAR3);	// Deprecated. Pending clean-up
	DATA_INTEGER(PAR4); // Daily to monthly resolution breakup point
	DATA_INTEGER(PAR5); // Location of the intrinsic effects
	DATA_INTEGER(PAR6); // Observation error

	// Data
	DATA_VECTOR(l1);
	DATA_IVECTOR(t1); // Month #
	DATA_VECTOR(d1);  // Day #/ # of days in that month
	DATA_IVECTOR(t2);
	DATA_VECTOR(d2);
	DATA_VECTOR(l2);
	DATA_IVECTOR(DUPE);
	DATA_IVECTOR(MNS);

	// Daily resolution
	DATA_IVECTOR(d1_d); // Day #
	DATA_IVECTOR(d2_d);
	DATA_IVECTOR(MNS_map);
	DATA_VECTOR(days);

	// Monthly Climate Cubic spline
	DATA_MATRIX(X_nao);
	DATA_MATRIX(X_ao);
	DATA_MATRIX(X_pna);
	DATA_SPARSE_MATRIX(S);
	DATA_IVECTOR(Sdim);
	DATA_MATRIX(prediction_design_matrix_nao);
	DATA_MATRIX(prediction_design_matrix_ao);
	DATA_MATRIX(prediction_design_matrix_pna);
	DATA_INTEGER(Pdim);
	DATA_MATRIX(p_dm_nao);
	DATA_MATRIX(p_dm_ao);
	DATA_MATRIX(p_dm_pna);

	// Reference GI
	DATA_INTEGER(REFLEN);
	DATA_INTEGER(REFDAYS);
	DATA_INTEGER(REFMON);
	DATA_INTEGER(REFREG);
	DATA_INTEGER(STRTMON);

	// Fixed effects
	PARAMETER(c);
	PARAMETER(k);

	// NAO effects
	PARAMETER_VECTOR(gc);
	PARAMETER_VECTOR(gk);
	PARAMETER(logNaoc);
	PARAMETER(logNaok);

	// Residual AO effects
	PARAMETER_VECTOR(gc_ao);
	PARAMETER_VECTOR(gk_ao);
	PARAMETER(logAoc);
	PARAMETER(logAok);

	// Residual PNA effects
	PARAMETER_VECTOR(gc_pna);
	PARAMETER_VECTOR(gk_pna);
	PARAMETER(logPnac);
	PARAMETER(logPnak);

	// Observation error
	PARAMETER(logSig1);
	PARAMETER_VECTOR(e_o);

	// Month effects
	PARAMETER(month_k);
	PARAMETER(month_mg);

	// Individuality
	PARAMETER_VECTOR(indl);
	PARAMETER(logvInd);

	// Prediction
	DATA_VECTOR(days_p);

	// Generalized logistic
	PARAMETER(lognu);

	// Regional effects
	PARAMETER(logReg);
	PARAMETER_VECTOR(reg_);
	DATA_IVECTOR(reg);

#ifdef _OPENMP
	parallel_accumulator<Type> nll(this);
#else
	Type nll = .0;
#endif

	// Guarded individuality
	if (CppAD::Variable(logvInd))
	{
		nll -= sum(dnorm(indl, Type(.0), exp(logvInd), true));
	}

	// Guarded month effect
	vector<Type> theta_l(12);
	if (CppAD::Variable(month_k))
	{
		Type k = exp(month_k) / (1 + exp(month_k)) * 2 * M_PI;
		Type A = exp(month_mg);
		for (int i = 0; i < 12; i++)
		{
			theta_l(i) = A * sin(2 * M_PI * i / 12.0 + k);
		}
	}
	else
	{
		theta_l = 0;
	}

	// Guarded NAO effect
	vector<Type> naoc = X_nao * gc;
	if (CppAD::Variable(logNaoc))
		nll -= 0.5 * Sdim[0] * logNaoc - 0.5 * exp(logNaoc) * GMRF(tmbutils::asSparseMatrix<Type>(S.block(0, 0, Sdim[0], Sdim[0]))).Quadform(gc);

	vector<Type> naok = X_nao * gk;
	if (CppAD::Variable(logNaok))
		nll -= 0.5 * Sdim[0] * logNaok - 0.5 * exp(logNaok) * GMRF(tmbutils::asSparseMatrix<Type>(S.block(0, 0, Sdim[0], Sdim[0]))).Quadform(gk);

	// Guarded residual AO effect
	vector<Type> aoc = X_ao * gc_ao;
	if (CppAD::Variable(logAoc))
		nll -= 0.5 * Sdim[1] * logAoc - 0.5 * exp(logAoc) * GMRF(tmbutils::asSparseMatrix<Type>(S.block(Sdim[0], Sdim[0], Sdim[1], Sdim[1]))).Quadform(gc_ao);

	vector<Type> aok = X_ao * gk_ao;
	if (CppAD::Variable(logAok))
		nll -= 0.5 * Sdim[1] * logAok - 0.5 * exp(logAok) * GMRF(tmbutils::asSparseMatrix<Type>(S.block(Sdim[0], Sdim[0], Sdim[1], Sdim[1]))).Quadform(gk_ao);

	// Guarded residual PNA effect
	vector<Type> pnac = X_pna * gc_pna;
	if (CppAD::Variable(logPnac))
		nll -= 0.5 * Sdim[2] * logPnac - 0.5 * exp(logPnac) * GMRF(tmbutils::asSparseMatrix<Type>(S.block(Sdim[0] + Sdim[1], Sdim[0] + Sdim[1], Sdim[2], Sdim[2]))).Quadform(gc_pna);

	vector<Type> pnak = X_pna * gk_pna;
	if (CppAD::Variable(logPnak))
		nll -= 0.5 * Sdim[2] * logPnak - 0.5 * exp(logPnak) * GMRF(tmbutils::asSparseMatrix<Type>(S.block(Sdim[0] + Sdim[1], Sdim[0] + Sdim[1], Sdim[2], Sdim[2]))).Quadform(gk_pna);

	// Observation errorr
	if (PAR6 == 1)
		nll -= sum(dnorm(e_o, Type(.0), exp(logSig1), true));

	// Regional effects
	if (CppAD::Variable(logReg))
		nll -= sum(dnorm(reg_, Type(.0), exp(logReg), true));

	// Local variables
	Type lt;
	Type baseline;
	Type effects; // intrinsic effects
	Type rg;
	int ds;
	int td;
	int tb;
	int id;
	int ms;

	// Likelihood
	for (int i = 0; i < l1.size(); i++)
	{
		ds = d2_d(i) - d1_d(i);
		lt = l1(i) + e_o(i);
		id = DUPE(i); // specimen id map
		ms = MNS(i);
		rg = reg_(reg(i));

		tb = t1(i);
		td = t2(i) - t1(i);
		if (td == 0)
		{
			// Released and recaptured in the same month
			effects = indl(id) + theta_l(ms) + rg;
			lt = tuna.growth_fn(lt, c + PAR5 * effects + naoc(tb) + aoc(tb) + pnac(tb), k + (1 - PAR5) * effects + naok(tb) + aok(tb) + pnak(tb), d2(i) - d1(i), PAR1, exp(lognu));
		}
		else
		{
			// Released and recaptured in different months
			effects = indl(id) + theta_l(ms) + rg;
			lt = tuna.growth_fn(lt, c + PAR5 * effects + naoc(tb) + aoc(tb) + pnac(tb), k + (1 - PAR5) * effects + naok(tb) + aok(tb) + pnak(tb), 1 - d1(i), PAR1, exp(lognu));
			for (int j = 1; j < td; j++)
			{
				effects = indl(id) + theta_l((ms + j) % 12) + rg;
				lt = tuna.growth_fn(lt, c + PAR5 * effects + naoc(tb + j) + aoc(tb + j) + pnac(tb + j), k + (1 - PAR5) * effects + naok(tb + j) + aok(tb + j) + pnak(tb + j), days(tb + j), PAR1, exp(lognu));
			}
			// Second partial month
			effects = indl(id) + theta_l((ms + td) % 12) + rg;
			lt = tuna.growth_fn(lt, c + PAR5 * effects + naoc(t2(i)) + aoc(t2(i)) + pnac(t2(i)), k + (1 - PAR5) * effects + naok(t2(i)) + aok(t2(i)) + pnak(t2(i)), d2(i), PAR1, exp(lognu));
		}
		nll -= dnorm(l2(i), lt, exp(logSig1), true);
	}

#ifdef _GRAPH
	// Prediction
	vector<Type> gam_c_n = prediction_design_matrix_nao * gc;
	vector<Type> gam_k_n = prediction_design_matrix_nao * gk;

	vector<Type> gam_c_a = prediction_design_matrix_ao * gc_ao;
	vector<Type> gam_k_a = prediction_design_matrix_ao * gk_ao;

	vector<Type> gam_c_p = prediction_design_matrix_pna * gc_pna;
	vector<Type> gam_k_p = prediction_design_matrix_pna * gk_pna;

	Type ref_len = REFLEN;
	int ref_tim = REFDAYS;

	effects = theta_l(REFMON) + reg_(REFREG);

	vector<Type> gam_l(Pdim);
	for (int i = 0; i < Pdim; i++)
	{
		gam_l(i) = tuna.growth_fn(ref_len, c + gam_c_n(i) + PAR5 * effects, k + gam_k_n(i) + (1 - PAR5) * effects, (Type)ref_tim / 30.0, PAR1, exp(lognu));
	}

	ADREPORT(gam_l);

	vector<Type> gam_la(Pdim);
	for (int i = 0; i < Pdim; i++)
	{
		gam_la(i) = tuna.growth_fn(ref_len, c + gam_c_a(i) + PAR5 * effects, k + gam_k_a(i) + (1 - PAR5) * effects, (Type)ref_tim / 30.0, PAR1, exp(lognu));
	}

	ADREPORT(gam_la);

	vector<Type> gam_lp(Pdim);
	for (int i = 0; i < Pdim; i++)
	{
		gam_lp(i) = tuna.growth_fn(ref_len, c + gam_c_p(i) + PAR5 * effects, k + gam_k_p(i) + (1 - PAR5) * effects, (Type)ref_tim / 30.0, PAR1, exp(lognu));
	}

	ADREPORT(gam_lp);
#endif
#ifdef _PREDICT
	vector<Type> gam_cn = p_dm_nao * gc;
	vector<Type> gam_ca = p_dm_ao * gc_ao;
	vector<Type> gam_cp = p_dm_pna * gc_pna;

	vector<Type> ref_p(days_p.size());
	Type tmp;
	for (int i = 0; i < days_p.size(); i++)
	{
		effects = theta_l((i + STRTMON) % 12) + reg_(REFREG);
		ref_p(i) = tuna.growth_fn(ref_len, c + gam_cn(i) + gam_ca(i) + gam_cp(i) + PAR5 * effects, k + (1 - PAR5) * effects, days_p(i), PAR1, exp(lognu));
	}

	ADREPORT(ref_p);
#endif
	return nll;
}
