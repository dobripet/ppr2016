#include "CubicSpline.h"
#include <amp.h>  
#include <amp_math.h>
#include <iostream>
//solve tridiagonal systems
#if defined(PARALLEL_AMP)
HRESULT Tridiagonal(int count, TGlucoseLevel **levels, std::vector<floattype> *z, std::vector<floattype> *x, std::vector<floattype> *y);
#endif
#if !defined(PARALLEL_AMP)
HRESULT Tridiagonal(int count, TGlucoseLevel **levels, std::vector<floattype> *z, std::vector<floattype> *h);
#endif
HRESULT CCubicSpline::Approximate(TApproximationParams *params) {
	size_t count;
	mEnumeratedLevels->GetLevelsCount(&count);
	TGlucoseLevel *levels;
	mEnumeratedLevels->GetLevels(&levels);
	//invalid levels
	if (count < 2) {
		return S_FALSE;
	}
#if !defined(PARALLEL_AMP)
	std::vector<floattype> z(count);
	std::vector<floattype> h(count);
	//get vectors
	Tridiagonal(count, &levels, &z, &h);
	//calc params and store them
	for (int i = 0; i < count - 1; i++) {
		floattype b = (-h[i] / 6)*z[i + 1] - (h[i] / 3)*z[i] + (1 / h[i])*(levels[i+1].level - levels[i].level);
		floattype d = (1 / (6 * h[i]))*(z[i + 1] - z[i]);
		cParams.insert(std::pair<floattype, cubic_params>(levels[i].datetime, cubic_params(levels[i].level, b, z[i] / 2, d)));
	}
#endif 
#if defined(PARALLEL_AMP)
	std::vector<floattype> levelsV(count);
	std::vector<floattype> datetimesV(count);
	std::vector<floattype> zV(count);
	cubic_params *par = (cubic_params*)malloc((count - 1) * sizeof(cubic_params));
	//get vectors
	Tridiagonal(count, &levels, &zV, &datetimesV, &levelsV);
	// Create C++ AMP objects.  
	concurrency::array_view<floattype, 1> y(count, levelsV);
	concurrency::array_view<floattype, 1> x(count, datetimesV);
	concurrency::array_view<floattype, 1> z(count, zV);
	concurrency::array_view<cubic_params, 1> parP(count - 1, par);
	concurrency::extent<1> ext(count - 1);
	parP.discard_data();
	concurrency::parallel_for_each(
		// Define the compute domain, which is the set of threads that are created.  
		parP.extent,
		// Define the code to run on each thread on the accelerator.  
		[=](concurrency::index<1> i) restrict(amp)
	{
		floattype h = x[i + 1] - x[i];

		parP[i].a = y[i];
		parP[i].b = (-h / 6)*z[i + 1] - (h / 3)*z[i] + (1 / h)*(y[i + 1] - y[i]);
		parP[i].c = z[i] / 2;
		parP[i].d = (1 / (6 * h))*(z[i + 1] - z[i]);	
	});
	parP.synchronize();
	//pick up results
	for (int i = 0; i < count - 1; i++) {
		cParams.insert(std::pair<floattype, cubic_params>(levels[i].datetime, parP[i]));
	}
	free(par);

#endif
	lastTime = levels[count - 1].datetime;
	return S_OK;
	}


	HRESULT  CCubicSpline::GetLevels(floattype desiredtime, floattype stepping, size_t count, floattype *levels, size_t *filled, size_t derivationorder) {
		*filled = 0;
		//Not aproximated yet
		if (cParams.empty()) {
			std::cerr << "CubicError: getLevels - Not aproximated yet!";
			return S_FALSE;
		}
		//Invalid time (after end or too much before start)
		if (lastTime < desiredtime || desiredtime + (stepping*count) < cParams.begin()->first) {
			std::cerr << "CubicError: getLevels - Invalid time (after or before segment)!";
			return S_FALSE;
		}
		//Move time to start
		floattype currenttime = desiredtime;
		while (currenttime < cParams.begin()->first) {
			currenttime += stepping;
			count--;
		}
		//Walk trough
		std::map<floattype, cubic_params>::iterator it;
		floattype res;
		if (derivationorder == 0) {
			while (currenttime <= lastTime && count > 0) {
				it = cParams.lower_bound(currenttime);
				//workaround lower bound;
				if (it != cParams.begin() && it->first != currenttime)
				{
					it--;
				}
				//workround for last value
				if (it == cParams.end()) {
					it--;
				}
				//y = a{i} + b{i}*(x-x{i})+ c{i}*(x-x{i})^2+d{i}*(x-x{i})^3
				res = it->second.a + it->second.b*(currenttime - it->first) + it->second.c*(currenttime - it->first)*(currenttime - it->first) + it->second.d*(currenttime - it->first)*(currenttime - it->first)*(currenttime - it->first);
				levels[(*filled)] = res;
				(*filled)++;
				currenttime += stepping;
				count--;
			}
		}
		else if (derivationorder == 1) {
			while (currenttime <= lastTime && count > 0) {
				it = cParams.lower_bound(currenttime);

				//workaround lower bound;
				if (it != cParams.begin() && it->first != currenttime)
				{
					it--;
				}
				//workround for last value
				if (it == cParams.end()) {
					it--;
				}
				//y = b{i} + 2*c{i}*(x-x{i}) + 3*d{i}*(x-x{i})^2
				res = it->second.b + 2 * it->second.c*(currenttime - it->first) + 3 * it->second.d*(currenttime - it->first)*(currenttime - it->first);
				levels[(*filled)] = res;
				(*filled)++;
				currenttime += stepping;
				count--;
			}
		}
		//wrong derivative order
		else {
			std::cerr << "CubicError: getLevels - Invalid derivative order!";
			return S_FALSE;
		}
		return S_OK;
	}

	CCubicSpline::CCubicSpline(IGlucoseLevels *levels) : CCommonApprox(levels) {
	}

	CCubicSpline::~CCubicSpline() {
	}

#if defined(PARALLEL_AMP)
	HRESULT Tridiagonal(int count, TGlucoseLevel **levels, std::vector<floattype> *z, std::vector<floattype> *x, std::vector<floattype> *y) {
		floattype *h = (floattype *)malloc((count - 1) * sizeof(floattype));
		floattype *b = (floattype *)malloc((count - 1) * sizeof(floattype));
		floattype *u = (floattype *)malloc((count - 1) * sizeof(floattype));
		floattype *v = (floattype *)malloc((count - 1) * sizeof(floattype));
		//prepare for 0,1 and count -1
		(*x)[count - 1] = (*levels)[count - 1].datetime;
		(*y)[count - 1] = (*levels)[count - 1].level;
		(*x)[0] = (*levels)[0].datetime;
		(*y)[0] = (*levels)[0].level;
		(*x)[1] = (*levels)[1].datetime;
		(*y)[1] = (*levels)[1].level;
		h[0] = (*levels)[1].datetime - (*levels)[0].datetime;
		b[0] = ((*levels)[1].level - (*levels)[0].level) / h[0];
		h[1] = (*levels)[2].datetime - (*levels)[1].datetime;
		b[1] = ((*levels)[2].level - (*levels)[1].level) / h[1];
		u[1] = 2.0*(h[0] + h[1]);
		v[1] = 6.0*(b[1] - b[0]);
		//filling x, y and calc h, b, u, v
		for (int i = 2; i < count - 1; i++) {
			(*x)[i] = (*levels)[i].datetime;
			(*y)[i] = (*levels)[i].level;
			h[i] = (*levels)[i + 1].datetime - (*levels)[i].datetime;
			b[i] = ((*levels)[i + 1].level - (*levels)[i].level) / h[i];
			u[i] = 2.0*(h[i] + h[i - 1]) - h[i - 1] * h[i - 1] / u[i - 1];
			v[i] = 6.0*(b[i] - b[i - 1]) - h[i - 1] * v[i - 1] / u[i - 1];
		}
		(*z)[count - 1] = 0;
		for (int i = count - 2; i > 0; i--) {
			(*z)[i] = (v[i] - h[i] * (*z)[i + 1]) / u[i];
		}
		(*z)[0] = 0;
		free(h);
		free(b);
		free(u);
		free(v);
		return S_OK;
	}
#endif
#if !defined(PARALLEL_AMP)
HRESULT Tridiagonal(int count, TGlucoseLevel **levels, std::vector<floattype> *z, std::vector<floattype> *h) {
	floattype *b = (floattype *)malloc((count - 1) * sizeof(floattype));
	floattype *u = (floattype *)malloc((count - 1) * sizeof(floattype));
	floattype *v = (floattype *)malloc((count - 1) * sizeof(floattype));
	//prepare for 0, 1
	(*h)[0] = (*levels)[1].datetime - (*levels)[0].datetime;
	b[0] = ((*levels)[1].level - (*levels)[0].level) / (*h)[0];
	(*h)[1] = (*levels)[2].datetime - (*levels)[1].datetime;
	b[1] = ((*levels)[2].level - (*levels)[1].level) / (*h)[1];
	u[1] = 2.0*((*h)[0] + (*h)[1]);
	v[1] = 6.0*(b[1] - b[0]);
	//calc h, b, u, v
	for (int i = 2; i < count - 1; i++) {
		(*h)[i] = (*levels)[i + 1].datetime - (*levels)[i].datetime;
		b[i] = ((*levels)[i + 1].level - (*levels)[i].level) / (*h)[i];
		u[i] = 2.0*((*h)[i] + (*h)[i - 1]) - (*h)[i - 1] * (*h)[i - 1] / u[i - 1];
		v[i] = 6.0*(b[i] - b[i - 1]) - (*h)[i - 1] * v[i - 1] / u[i - 1];
	}
	(*z)[count - 1] = 0;
	for (int i = count - 2; i > 0; i--) {
		(*z)[i] = (v[i] - (*h)[i] * (*z)[i + 1]) / u[i];
	}
	(*z)[0] = 0;
	free(b);
	free(u);
	free(v);
	return S_OK;
}
#endif