#include "CubicSpline.h"
#include <amp.h>  
#include <amp_math.h>
#include <iostream>


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
	/*invalid levels*/
	if (count < 2) {
		return S_FALSE;
	}
#if !defined(PARALLEL_AMP)
	/*//size_t n = count - 1;
	floattype *h = (floattype*)malloc(count * sizeof(floattype));
	floattype *alpha = (floattype*)malloc(count * sizeof(floattype));
	//x = datetime, a = level
	//hi = x{i+1} - xi
	h[0] = levels[1].datetime - levels[0].datetime;


	//solve tridiagonal linear system
	floattype *l = (floattype*)malloc(count * sizeof(floattype));
	floattype *z = (floattype*)malloc(count * sizeof(floattype));
	floattype *mu = (floattype*)malloc(count * sizeof(floattype));
	l[0] = 1;
	z[0] = 0;
	mu[0] = 0;
	for (size_t i = 1; i < count-1; i++) {
		h[i] = levels[i + 1].datetime - levels[i].datetime;
		alpha[i] = (3 / h[i])*(levels[i + 1].level - levels[i].level) - (3 / h[i - 1])*(levels[i].level - levels[i - 1].level);


		l[i] = 2 * (levels[i + 1].datetime - levels[i - 1].datetime) - h[i - 1] * mu[i - 1];
		mu[i] = h[i] / l[i];
		z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
	}

	floattype *c = (floattype*)malloc(count * sizeof(floattype));
	//floattype *b = (floattype*)malloc(count * sizeof(floattype));
	//floattype *d = (floattype*)malloc(count * sizeof(floattype));
	l[count-1] = 1;
	z[count-1] = 0;
	c[count-1] = 0;
	floattype b;
	floattype d;

	for (long i = (long)(count-2); i >= 0; i--) {
		c[i] = z[i] - mu[i] * c[i + 1];
		b = (levels[i + 1].level - levels[i].level) / h[i] - h[i] * (c[i + 1] + 2 * c[i]) / 3;
		d = (c[i + 1] - c[i]) / (3 * h[i]);
		cParams.insert(std::pair<floattype, cubic_params>(levels[i].datetime, cubic_params(levels[i].level, b, c[i], d)));
	}
	//add last timestamp to check bounds in getLevels
	lastTime = levels[count - 1].datetime;
	free(h);
	free(alpha);
	free(z);
	free(l);
	free(mu);
	free(c);*/


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
	lastTime = levels[count - 1].datetime;
	free(par);

#endif
		//printf("aproximuju %d %f %f\n", count, cParams.begin()->first, lastTime);
		return S_OK;
	}
		HRESULT  CCubicSpline::GetLevels(floattype desiredtime, floattype stepping, size_t count, floattype *levels, size_t *filled, size_t derivationorder) {
		*filled = 0;
		/*Not aproximated yet*/
		if (cParams.empty()) {
			printf("err empty cparams\n");
			return S_FALSE;
		}
		//printf("vracim hodnoty %f %f %f", desiredtime, cParams.begin()->first, cParams.rbegin()->first);
		//cParams.end()->first
		//printf("test bond %f\n", cParams.lower_bound(desiredtime)->first);
		//printf("co se deje?\n");
		//std::map<floattype, cubic_params>::iterator itlow;
		//itlow = cParams.lower_bound(desiredtime);
		/*Invalid time (after end or too much before start)*/
		if (lastTime < desiredtime || desiredtime + (stepping*count) < cParams.begin()->first) {
			printf("err time\n");
			return S_FALSE;
		}
		//printf("prdelni testik %d\n", (cParams.find(desiredtime)));
		//std::map<floattype, cubic_params>::iterator it = cParams.find(desiredtime);
		/*Move time to start*/
		floattype currenttime = desiredtime;
		while (currenttime < cParams.begin()->first) {
			currenttime += stepping;
			count--;
		}
		/*Walk trough*/
		std::map<floattype, cubic_params>::iterator it;
		floattype res;
		if (derivationorder == 0) {
			//printf("tu jsem1 %lu %f %f\n", *filled, currenttime, lastTime);
			while (currenttime <= lastTime && count > 0) {
				it = cParams.lower_bound(currenttime);
				//workaround lower bound;
				if (it != cParams.begin() && it->first != currenttime)
				{
					it--;
				}
				/*workround for last value*/
				if (it == cParams.end()) {
					it--;
				}
				printf("tu jsem %f %f %f \n", currenttime, it->first, currenttime - it->first);
				//y = a{i} + b{i}*(x-x{i})+ c{i}*(x-x{i})^2+d{i}*(x-x{i})^3
				res = it->second.a + it->second.b*(currenttime - it->first) + it->second.c*(currenttime - it->first)*(currenttime - it->first) + it->second.d*(currenttime - it->first)*(currenttime - it->first)*(currenttime - it->first);
				levels[(*filled)] = res;
				(*filled)++;
				currenttime += stepping;
				count--;
				//printf("tu jsem %f %f %f %f %f %f\n", res, it->second.a, it->second.b, it->second.c, it->second.d, currenttime - it->first);
				//printf("tu jsem2 %lu\n", *filled);
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
				/*workround for last value*/
				if (it == cParams.end()) {
					it--;
				}
				//y = b{i} + 2*c{i}*(x-x{i}) + 3*d{i}*(x-x{i})^2
				res = it->second.b + 2 * it->second.c*(currenttime - it->first) + 3 * it->second.d*(currenttime - it->first)*(currenttime - it->first);
				levels[(*filled)] = res;
				(*filled)++;
				currenttime += stepping;
				count--;
				//printf("tu jsem roblhl %lu\n", *filled);
			}
		}
		/*wrong derivative order*/
		else {
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