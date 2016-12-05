#include "CatmullRomSpline.h"
#include <amp.h>  
#include <amp_math.h>
#include <iostream>

//calc sqrt alpha distance
floattype getDT(TGlucoseLevel p0, TGlucoseLevel p1) {
	const floattype exponent = 0.25;//cetripetal => alhpa = 0.5 and sqrt => 0.5
	// (sqrt((x{i+1}-x{i})^2 + (y{i+1}-y{i})^2))^alhpa
	return std::pow(((p1.datetime - p0.datetime)*(p1.datetime - p0.datetime) + (p1.level - p0.level)*(p1.level - p0.level)), exponent);
};
//calc sqrt alpha distance
floattype getDT(floattype x0, floattype x1, floattype y0, floattype y1) restrict(amp){
	const floattype exponent = 0.25;//cetripetal => alhpa = 0.5 and sqrt => 0.5
	//t{i+1} = t{i} + (sqrt((x{i+1}-x{i})^2 + (y{i+1}-y{i})^2))^alhpa
	return concurrency::precise_math::pow(((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0)), exponent);
};
//calc tangents with rescale to parametrization [0,1]
floattype getTan(floattype y0, floattype y1, floattype y2, floattype dt0, floattype dt1) restrict(cpu, amp) {
	return ((y1 - y0) / dt0 - (y2 - y0) / (dt0 + dt1) + (y2 - y1) / dt1) *dt1;
}
HRESULT CCatmullRomSpline::Approximate(TApproximationParams *params) {

	size_t count;
	mEnumeratedLevels->GetLevelsCount(&count);
	TGlucoseLevel *levels;
	mEnumeratedLevels->GetLevels(&levels);
	//invalid levels
	if (count < 2) {
		return S_FALSE;
	}
#if !defined(PARALLEL_AMP)
	//inner values
	for (size_t i = 1; i < count - 2; i++) {
		floattype dt0 = getDT(levels[i - 1], levels[i]);
		floattype dt1 = getDT(levels[i], levels[i + 1]);
		floattype dt2 = getDT(levels[i + 1], levels[i + 2]); 
		//cals tangents
		floattype tan1 = getTan(levels[i - 1].level, levels[i].level, levels[i + 1].level, dt0, dt1);
		floattype tan2 = getTan(levels[i].level, levels[i + 1].level, levels[i + 2].level, dt1, dt2);
		//c{i} = 3*(p{i+1}-p{i})-tan{2}-2*tan{1}
		floattype c = 3 * (levels[i + 1].level - levels[i].level) - tan2 - 2 * tan1;
		//d{i} = -2*(p{i+1}-p{i})+tan{2}+tan{1}
		floattype d = -2 * (levels[i + 1].level - levels[i].level) + tan2 + tan1;
		// a = y{i}, b = tan{1}
		crParams.insert(std::pair<floattype, cr_params>(levels[i].datetime, cr_params(levels[i].level, tan1, c, d)));
	}
	//boundary, adding one identical point to head and tail
	floattype dt0 = 1.0;
	floattype dt1 = getDT(levels[0], levels[1]);
	floattype dt2 = getDT(levels[1], levels[2]);
	floattype tan1 = getTan(levels[0].level, levels[0].level, levels[1].level, dt0, dt1);
	floattype tan2 = getTan(levels[0].level, levels[1].level, levels[2].level, dt1, dt2);
	floattype c = 3 * (levels[1].level - levels[0].level) - tan2 - 2 * tan1;
	floattype d = -2 * (levels[1].level - levels[0].level) + tan2 + tan1;
	crParams.insert(std::pair<floattype, cr_params>(levels[0].datetime, cr_params(levels[0].level, tan1, c, d)));
	dt0 = getDT(levels[count - 3], levels[count - 2]);
	dt1 = getDT(levels[count - 2], levels[count - 1]);
	dt2 = 1;
	tan1 = getTan(levels[count - 3].level, levels[count - 2].level, levels[count - 1].level, dt0, dt1);
	tan2 = getTan(levels[count - 2].level, levels[count - 1].level, levels[count - 1].level, dt1, dt2);
	c = 3 * (levels[count - 1].level - levels[count - 2].level) - tan2 - 2 * tan1;
	d = -2 * (levels[count - 1].level - levels[count - 2].level) + tan2 + tan1;
	crParams.insert(std::pair<floattype, cr_params>(levels[count - 2].datetime, cr_params(levels[count-2].level, tan1, c, d)));
#endif 
#if defined(PARALLEL_AMP)
	std::vector<floattype> levelsV(count);
	std::vector<floattype> datetimesV(count);
	cr_params *par = (cr_params*)malloc((count-1) * sizeof(cr_params));
	//prepare vectors
	for (size_t i = 0; i < count; i++) {
		levelsV[i] = levels[i].level;
		datetimesV[i] = levels[i].datetime;
	}
	// Create C++ AMP objects.  
	concurrency::extent<1> ext(count - 3);
	concurrency::array_view<floattype, 1> y(count, levelsV);
	concurrency::array_view<floattype, 1> x(count, datetimesV);
	concurrency::array_view<cr_params, 1> parP(count - 1, par);
	parP.discard_data();
	concurrency::parallel_for_each(ext, [=](concurrency::index<1> idx) restrict(amp) {
		int i = idx[0] + 1;
		floattype dt0 = getDT(x[i-1], x[i], y[i-1], y[i]);
		floattype dt1 = getDT(x[i], x[i + 1], y[i], y[i + 1]);
		floattype dt2 = getDT(x[i+1], x[i + 2], y[i+1], y[i + 2]);
		floattype tan1 = getTan(y[i - 1], y[i], y[i + 1], dt0, dt1);
		floattype tan2 = getTan(y[i], y[i + 1], y[i + 2], dt1, dt2);
		// a = y{i}, b = tan{1}
		parP[i].a = y[i]; 
		parP[i].b = tan1;
		//c{i} = 3*(p{i+1}-p{i})-tan{2}-2*tan{1}
		parP[i].c = 3 * (y[i + 1] - y[i]) - tan2 - 2 * tan1;
		//d{i} = -2*(p{i+1}-p{i})+tan{2}+tan{1}
		parP[i].d = -2 * (y[i + 1] -y[i]) + tan2 + tan1;
		
	});
	parP.synchronize();
	//boundary, adding one identical point to head and tail
	floattype dt0 = 1.0;
	floattype dt1 = getDT(levels[0], levels[1]);
	floattype dt2 = getDT(levels[1], levels[2]);
	floattype tan1 = getTan(levels[0].level, levels[0].level, levels[1].level, dt0, dt1);
	floattype tan2 = getTan(levels[0].level, levels[1].level, levels[2].level, dt1, dt2);
	parP[0].a = levels[0].level;
	parP[0].b = tan1;
	parP[0].c = 3 * (levels[1].level - levels[0].level) - tan2 - 2 * tan1;
	parP[0].d = -2 * (levels[1].level - levels[0].level) + tan2 + tan1;
	dt0 = getDT(levels[count - 3], levels[count - 2]);
	dt1 = getDT(levels[count - 2], levels[count - 1]);
	dt2 = 1;
	tan1 = getTan(levels[count - 3].level, levels[count - 2].level, levels[count - 1].level, dt0, dt1);
	tan2 = getTan(levels[count - 2].level, levels[count - 1].level, levels[count - 1].level, dt1, dt2);
	parP[count - 2].a = levels[count - 2].level;
	parP[count - 2].b = tan1;
	parP[count - 2].c = 3 * (levels[count - 1].level - levels[count - 2].level) - tan2 - 2 * tan1;
	parP[count - 2].d = -2 * (levels[count - 1].level - levels[count - 2].level) + tan2 + tan1;
	//add to maps
	for (int i = 0; i < count - 1; i++) {
		crParams.insert(std::pair<floattype, cr_params>(levels[i].datetime, parP[i]));
	}
	free(par);
#endif
	//add last timestamp to check bounds in getLevels
	lastTime = levels[count - 1].datetime;
	return S_OK;
}
HRESULT  CCatmullRomSpline::GetLevels(floattype desiredtime, floattype stepping, size_t count, floattype *levels, size_t *filled, size_t derivationorder) {
	*filled = 0;
	//Not aproximated yet
	if (crParams.empty()) {
		std::cerr << "CatmullRomError: getLevels - Not aproximated yet!";
		return S_FALSE;
	}
	//Invalid time (after end or too much before start)
	if (lastTime < desiredtime || desiredtime + (stepping*count) < crParams.begin()->first) {
		std::cerr << "CatmullRomError: getLevels - Invalid time (after or before segment)!";
		return S_FALSE;
	}
	//Move time to start
	floattype currenttime = desiredtime;
	while (currenttime < crParams.begin()->first) {
		currenttime += stepping;
		count--;
	}
	//Walk trough
	std::map<floattype, cr_params>::iterator it;
	floattype res;
	floattype nextTime;
	if (derivationorder == 0) {
		while (currenttime <= lastTime && count > 0) {
			it = crParams.lower_bound(currenttime);
			//workaround lower bound
			if (it != crParams.begin() && it->first != currenttime) {
				it--;

			}
			//workround for last value
			if (it == crParams.end()) {
				it--;					
			}
			auto it2 = it;
			it2++;
			nextTime = it2->first;
			if (it2 == crParams.end()) {
				nextTime = lastTime;
			}
			//time normalization
			floattype time = (currenttime - it->first) / (nextTime - it->first);
			//y = a{i} + b{i}*(x-x{i})+ c{i}*(x-x{i})^2+d{i}*(x-x{i})^3
			res = it->second.a + it->second.b*time + it->second.c*time*time + it->second.d*time*time*time;
			levels[(*filled)] = res;
			(*filled)++;
			currenttime += stepping;
			count--;
		}
	}
	else if (derivationorder == 1) {
		while (currenttime <= lastTime && count > 0) {
			it = crParams.lower_bound(currenttime);
			//workaround lower bound
			if (it != crParams.begin() && it->first != currenttime) {
				it--;

			}
			//workround for last value
			if (it == crParams.end()) {
				it--;
			}
			auto it2 = it;
			it2++;
			nextTime = it2->first;
			if (it2 == crParams.end()) {
				nextTime = lastTime;
			}
			//time normalization
			floattype time = (currenttime - it->first) / (nextTime - it->first);
			//y = b{i} + 2*c{i}*(x-x{i}) + 3*d{i}*(x-x{i})^2
			res = it->second.b + 2 * it->second.c*time + 3 * it->second.d*time*time;
			levels[(*filled)] = res;
			(*filled)++;
			currenttime += stepping;
			count--;
		}
	}
	//wrong derivative order
	else {
		std::cerr << "CatmullRomError: getLevels - Invalid derivative order!";
		return S_FALSE;
	}
	return S_OK;
}

CCatmullRomSpline::CCatmullRomSpline(IGlucoseLevels *levels) : CCommonApprox(levels) {
}

CCatmullRomSpline::~CCatmullRomSpline() {
}