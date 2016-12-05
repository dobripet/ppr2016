#include "QuadraticSpline.h"
#include <iostream>  

HRESULT CQuadraticSpline::Approximate(TApproximationParams *params) {
	size_t count;
	mEnumeratedLevels->GetLevelsCount(&count);
	TGlucoseLevel *levels;
	mEnumeratedLevels->GetLevels(&levels);	
	//invalid levels
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
	//add last timestamp to check bounds in getLevels
	lastTime = levels[count - 1].datetime;
	free(z);
	return S_OK;
}
HRESULT  CQuadraticSpline::GetLevels(floattype desiredtime, floattype stepping, size_t count, floattype *levels, size_t *filled, size_t derivationorder) {
	*filled = 0;
	//Not aproximated yet
	if (qParams.empty()) {
		std::cerr << "QuadraticError: getLevels - Not aproximated yet!";
		return S_FALSE;
	}
	//Invalid time (after end or too much before start)
	if (lastTime < desiredtime || desiredtime+(stepping*count) < qParams.begin()->first) {
		std::cerr << "QuadraticError: getLevels - Invalid time (after or before segment)!";
		return S_FALSE;
	}
	//Move time to start
	floattype currenttime = desiredtime;
	while (currenttime < qParams.begin()->first) {
		currenttime += stepping;
		count--;
	}
	//Walk trough
	std::map<floattype, quadratic_params>::iterator it;
	floattype res;
	if (derivationorder == 0) {
		while (currenttime <= lastTime && count > 0) {
			it = qParams.lower_bound(currenttime);

			//workaround lower bound;
			if (it != qParams.begin() && it->first != currenttime)
			{
				it--;
			}
			//workround for last value
			if (it == qParams.end()) {
				it--;
			}
			res = it->second.a + it->second.b*(currenttime - it->first) + it->second.c*(currenttime - it->first)*(currenttime - it->first);
			levels[(*filled)] = res;
			(*filled)++;
			currenttime += stepping;
			count--;
		}
	}
	else if (derivationorder == 1){
		while (currenttime <= lastTime && count > 0) {
			it = qParams.lower_bound(currenttime);
			//workaround lower bound;
			if (it != qParams.begin() && it->first != currenttime)
			{
				it--;
			}
			//workround for last value
			if (it == qParams.end()) {
				it--;
			}
			res = it->second.b + 2*it->second.c*(currenttime - it->first);
			levels[(*filled)] = res;
			(*filled)++;
			currenttime += stepping;
			count--;
		}
	}
	//wrong derivative order
	else {
		std::cerr << "QuadraticError: getLevels - Invalid derivative order!";
		return S_FALSE;
	}
	return S_OK;
}

CQuadraticSpline::CQuadraticSpline(IGlucoseLevels *levels) : CCommonApprox(levels) {
}

CQuadraticSpline::~CQuadraticSpline() {
}
