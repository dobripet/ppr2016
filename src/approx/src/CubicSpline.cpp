#include "CubicSpline.h"


HRESULT CCubicSpline::Approximate(TApproximationParams *params) {
	size_t count;
	mEnumeratedLevels->GetLevelsCount(&count);
	TGlucoseLevel *levels;
	mEnumeratedLevels->GetLevels(&levels);	
	/*invalid levels*/
	if (count < 2) {
		return S_FALSE;
	}
	//size_t n = count - 1;
	floattype *h = (floattype*)malloc(count * sizeof(floattype));
	floattype *alpha = (floattype*)malloc(count * sizeof(floattype));
	//x = datetime, a = level
	//hi = x{i+1} - xi
	h[0] = levels[1].datetime - levels[0].datetime;
	for (size_t i = 1; i < count-1; i++) {
		h[i] = levels[i + 1].datetime - levels[i].datetime;
		alpha[i] = (3 / h[i])*(levels[i + 1].level - levels[i].level) - (3 / h[i - 1])*(levels[i].level - levels[i - 1].level);
	}
	//solve tridiagonal linear system

	floattype *l = (floattype*)malloc(count * sizeof(floattype));
	floattype *z = (floattype*)malloc(count * sizeof(floattype));
	floattype *mu = (floattype*)malloc(count * sizeof(floattype));
	l[0] = 1;
	z[0] = 0;
	mu[0] = 0;

	for (size_t i = 1; i < count-1; i++) {
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
		cParams[levels[i].datetime] = cubic_params(levels[i].level, b, c[i], d);
	}
	/*add last timestamp to check bounds in getLevels*/
	lastTime = levels[count - 1].datetime;
	free(h);
	free(alpha);
	free(z);
	free(l);
	free(mu);
	free(c);

	printf("aproximuju %d %f %f\n", count, cParams.begin()->first, lastTime);
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
	if (lastTime < desiredtime || desiredtime+(stepping*count) < cParams.begin()->first) {
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
			/*workround for last value*/
			if (it == cParams.end()) {
				it--;
			}
			res = it->second.a + it->second.b*(currenttime - it->first) + it->second.c*(currenttime - it->first)*(currenttime - it->first) + it->second.d*(currenttime - it->first)*(currenttime - it->first)*(currenttime - it->first);
			levels[(*filled)] = res;
			(*filled)++;
			currenttime += stepping;
			count--;
			printf("tu jsem2 %lu\n", *filled);
		}
	}
	else if (derivationorder == 1){
		while (currenttime <= lastTime && count > 0) {
			it = cParams.lower_bound(currenttime);
			res = it->second.b + 2*it->second.c*(currenttime - it->first) + 3*it->second.d*(currenttime - it->first)*(currenttime - it->first);
			levels[(*filled)] = res;
			(*filled)++;
			currenttime += stepping;
			count--;
			printf("tu jsem roblhl %lu\n", *filled);
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
