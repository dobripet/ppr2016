#include "QuadraticSpline.h"


HRESULT CQuadraticSpline::Approximate(TApproximationParams *params) {
	size_t count;
	mEnumeratedLevels->GetLevelsCount(&count);
	TGlucoseLevel *levels;
	mEnumeratedLevels->GetLevels(&levels);	
	/*invalid levels*/
	if (count < 2) {
		return S_FALSE;
	}
	//size_t n = count - 1;
	floattype *z = (floattype*)malloc(count * sizeof(floattype));	
	//prepare z
	//t is datetime, y is level
	z[0] = 0;
	for (size_t i = 0; i < count - 1; i++) {
		//z{i+1} = -z{i} + 2*((y{i+1}-y{i})/(t{i+1}-t{i}))
		z[i + 1] = -z[i] + 2 * ((levels[i + 1].level - levels[i].level) / (levels[i + 1].datetime - levels[i].datetime));
		//a{i} is y{i}, b{i} is z{i} and c{i}=(z{i+1}-z{i})/(2*(t{i+1}-t{i}))
		qParams.insert(std::pair<floattype,quadratic_params>(levels[i].datetime, quadratic_params(levels[i].level, z[i], (z[i+1]-z[i])/(2* (levels[i + 1].datetime - levels[i].datetime)))));
	}
	/*add last timestamp to check bounds in getLevels*/
	lastTime = levels[count - 1].datetime;
	free(z);
	//printf("aproximuju %d %f %f\n", count, cParams.begin()->first, lastTime);
	return S_OK;
}
HRESULT  CQuadraticSpline::GetLevels(floattype desiredtime, floattype stepping, size_t count, floattype *levels, size_t *filled, size_t derivationorder) {
	*filled = 0;
	/*Not aproximated yet*/
	if (qParams.empty()) {
		printf("err empty qparams\n");
		return S_FALSE;
	}
	//printf("vracim hodnoty %f %f %f", desiredtime, cParams.begin()->first, cParams.rbegin()->first);
	//cParams.end()->first
	//printf("test bond %f\n", cParams.lower_bound(desiredtime)->first);
	//printf("co se deje?\n");
	//std::map<floattype, cubic_params>::iterator itlow;
	//itlow = cParams.lower_bound(desiredtime);
	/*Invalid time (after end or too much before start)*/
	if (lastTime < desiredtime || desiredtime+(stepping*count) < qParams.begin()->first) {
		printf("err time\n");
		return S_FALSE;
	}
	//printf("prdelni testik %d\n", (cParams.find(desiredtime)));
	//std::map<floattype, cubic_params>::iterator it = cParams.find(desiredtime);
	/*Move time to start*/
	floattype currenttime = desiredtime;
	while (currenttime < qParams.begin()->first) {
		currenttime += stepping;
		count--;
	}
	/*Walk trough*/
	std::map<floattype, quadratic_params>::iterator it;
	floattype res;
	if (derivationorder == 0) {
		//printf("tu jsem1 %lu %f %f\n", *filled, currenttime, lastTime);
		while (currenttime <= lastTime && count > 0) {
			it = qParams.lower_bound(currenttime);
			/*workround for last value*/
			if (it == qParams.end()) {
				it--;
			}
			res = it->second.a + it->second.b*(currenttime - it->first) + it->second.c*(currenttime - it->first)*(currenttime - it->first);
			levels[(*filled)] = res;
			(*filled)++;
			currenttime += stepping;
			count--;
			//printf("tu jsem2 %lu\n", *filled);
		}
	}
	else if (derivationorder == 1){
		while (currenttime <= lastTime && count > 0) {
			it = qParams.lower_bound(currenttime);
			/*workround for last value*/
			if (it == qParams.end()) {
				it--;
			}
			res = it->second.b + 2*it->second.c*(currenttime - it->first);
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

CQuadraticSpline::CQuadraticSpline(IGlucoseLevels *levels) : CCommonApprox(levels) {
}

CQuadraticSpline::~CQuadraticSpline() {
}
