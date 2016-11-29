#include "AkimaSpline.h"


HRESULT CAkimaSpline::Approximate(TApproximationParams *params) {
	size_t count;
	mEnumeratedLevels->GetLevelsCount(&count);
	TGlucoseLevel *levels;
	mEnumeratedLevels->GetLevels(&levels);	
	/*invalid levels*/
	if (count < 2) {
		return S_FALSE;
	}
	//size_t n = count - 1;
	floattype *t = (floattype*)malloc(count * sizeof(floattype));
	floattype *m = (floattype*)malloc((count+1) * sizeof(floattype));
	//x is datetime, y is level
	//boundary values
	//m{i} = (y{i+1}-y{i})/(x{i+1}-x{i})
	m[count - 2] = (levels[count - 1].level - levels[count - 2].level) / (levels[count - 1].datetime - levels[count - 2].datetime);
	m[count - 3] = (levels[count - 2].level - levels[count - 3].level) / (levels[count - 2].datetime - levels[count - 3].datetime);	
	//m{i-1} = 2m{i-2} - m{i-3}
	m[count-1] = 2 * m[count - 2] - m[count - 3];
	//m{i} = 2m{i-1} - m{i-2}
	m[count] = 2 * m[count-1] - m[count - 2];	
	//t{i+1}
	//zero check
	if (m[count] == m[count - 1] && m[count - 2] == m[count - 3]) {
		t[count - 1] = 0.5 * (m[count - 1] + m[count - 2]);
	}
	else {
		t[count - 1] = (abs(m[count] - m[count - 1])*m[count - 2] + abs(m[count - 2] - m[count - 3])*m[count - 1]) / (abs(m[count] - m[count - 1]) + abs(m[count - 2] - m[count - 3]));
	}
	for (size_t i = count-2; i > 1; i--) {
		//zero check, not a function
		if (levels[i + 1].datetime == levels[i].datetime) {
			return S_FALSE;
		}
		floattype dx = levels[i + 1].datetime - levels[i].datetime;
		//m{i} = (y{i+1}-y{i})/(x{i+1}-x{i})
		m[i] = (levels[i + 1].level - levels[i].level) / dx;
		//zero check
		if (m[i+1] == m[i] && m[i-1] == m[i-2]) {
			t[i] = 0.5 * (m[i] + m[i-1]);
		}
		else {
			//t{i} = (abs(m{i + 1} - m{i})*m{i - 1} + abs(m{i - 1} - m{i - 2})*m{i}) / (abs(m{i + 1} - m{i}) + abs(m{i - 1} - m{i - 2}))
			t[i] = (abs(m[i+1] - m[i])*m[i-1] + abs(m[i-1] - m[i-2])*m[i]) / (abs(m[i + 1] - m[i]) + abs(m[i - 1] - m[i - 2]));
		}
		//a{i} is y{i}, b{i} is t{i}, c{i} is (3 * m{i} - 2 * t{i} - t{i+1}) / (x{i+1} -x{i}) and d{i} is (t{i} + t{i+1} - 2*m{i}) / (x{i+1} -x{i})^2
		aParams.insert(std::pair<floattype, akima_params>(levels[i].datetime, akima_params(levels[i].level, t[i],
			(3*m[i]-2*t[i]-t[i+1])/dx,( t[i] + t[i + 1]-2*m[i])/(dx*dx))));
	}
	//boundary values
	m[0] = (levels[1].level - levels[0].level) / (levels[1].datetime - levels[0].datetime);
	m[1] = (levels[2].level - levels[1].level) / (levels[2].datetime - levels[1].datetime);
	//m{-1} = 2m{0} - 2m{1}
	floattype mm1 = 2 * m[0] - m[1];
	//m{-2} = 2m{-1} - m{0}
	floattype mm2 = 2 * mm1 - m[0];
	//t{1}
	//zero check
	if (m[0] == mm1 && m[2] == m[1]) {
		t[1] = 0.5 * (m[1] + m[0]);
	}
	else {
		t[1] = (abs(m[2] - m[1])*m[0] + abs(m[0] - mm1)*m[1]) / (abs(m[2] - m[1]) + abs(m[0] - mm1));
	}
	aParams.insert(std::pair<floattype, akima_params>(levels[1].datetime, 
		akima_params(levels[1].level, t[1], (3 * m[1] - 2 * t[1] - t[2]) / (levels[2].datetime - levels[1].datetime), 
		(t[1] + t[2] - 2 * m[1]) / ((levels[2].datetime - levels[1].datetime)*(levels[2].datetime - levels[1].datetime)))));
	//t{0}
	//zero check
	if (mm1 == mm2 && m[1] == m[0]) {
		t[0] = 0.5 * (m[0] + mm1);
	}
	else {
		t[0] = (abs(m[1] - m[0])*mm1 + abs(mm1 - mm2)*m[0]) / (abs(m[1] - m[0]) + abs(mm1 - mm2));
	}
	aParams.insert(std::pair<floattype, akima_params>(levels[0].datetime,
		akima_params(levels[0].level, t[0], (3 * m[0] - 2 * t[0] - t[1]) / (levels[1].datetime - levels[0].datetime),
		(t[0] + t[1] - 2 * m[0]) / ((levels[1].datetime - levels[0].datetime)*(levels[1].datetime - levels[0].datetime)))));



	/*add last timestamp to check bounds in getLevels*/
	lastTime = levels[count - 1].datetime;
	free(m);
	free(t);
	//printf("aproximuju %d %f %f\n", count, cParams.begin()->first, lastTime);
	return S_OK;
}
HRESULT  CAkimaSpline::GetLevels(floattype desiredtime, floattype stepping, size_t count, floattype *levels, size_t *filled, size_t derivationorder) {
	*filled = 0;
	/*Not aproximated yet*/
	if (aParams.empty()) {
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
	if (lastTime < desiredtime || desiredtime + (stepping*count) < aParams.begin()->first) {
		printf("err time\n");
		return S_FALSE;
	}
	//printf("prdelni testik %d\n", (cParams.find(desiredtime)));
	//std::map<floattype, cubic_params>::iterator it = cParams.find(desiredtime);
	/*Move time to start*/
	floattype currenttime = desiredtime;
	while (currenttime < aParams.begin()->first) {
		currenttime += stepping;
		count--;
	}
	/*Walk trough*/
	std::map<floattype, akima_params>::iterator it;
	floattype res;
	if (derivationorder == 0) {
		//printf("tu jsem1 %lu %f %f\n", *filled, currenttime, lastTime);
		while (currenttime <= lastTime && count > 0) {
			it = aParams.lower_bound(currenttime);
			/*workround for last value*/
			if (it == aParams.end()) {
				it--;
			}
			//y = a{i} + b{i}*(x-x{i})+ c{i}*(x-x{i})^2+d{i}*(x-x{i})^3
			res = it->second.a + it->second.b*(currenttime - it->first) + it->second.c*(currenttime - it->first)*(currenttime - it->first) + it->second.d*(currenttime - it->first)*(currenttime - it->first)*(currenttime - it->first);
			levels[(*filled)] = res;
			(*filled)++;
			currenttime += stepping;
			count--;
			//printf("tu jsem2 %lu\n", *filled);
		}
	}
	else if (derivationorder == 1) {
		while (currenttime <= lastTime && count > 0) {
			it = aParams.lower_bound(currenttime);
			/*workround for last value*/
			if (it == aParams.end()) {
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

CAkimaSpline::CAkimaSpline(IGlucoseLevels *levels) : CCommonApprox(levels) {
}

CAkimaSpline::~CAkimaSpline() {
}
