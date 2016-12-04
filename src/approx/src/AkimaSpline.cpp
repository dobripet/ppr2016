#include "AkimaSpline.h"
#include <amp.h>  
#include <amp_math.h>
#include <iostream>
int a();

floattype get_p3(floattype &xx, floattype &yy, floattype &ti, floattype &ti1) restrict(amp) {
	return (ti + ti1 - 2 * yy / xx) / concurrency::precise_math::pow(xx, 2);
}
floattype get_p2(floattype &xx, floattype &yy, floattype &ti, floattype &ti1) restrict(cpu, amp) {
	return (3 * yy / xx - 2 * ti - ti1) / xx;
}

floattype get_t(floattype m0, floattype m1, floattype m2, floattype m3) restrict(amp) {
	floattype numerator, denominator;
	numerator = concurrency::precise_math::fabs(m3 - m2) * m1 + concurrency::precise_math::fabs(m1 - m0) * m2;
	denominator = concurrency::precise_math::fabs(m3 - m2) + concurrency::precise_math::fabs(m1 - m0);
	return (denominator == 0) ? 0.5 * (m2 + m1) : numerator / denominator;
}

HRESULT CAkimaSpline::Approximate(TApproximationParams *params) {

	size_t count;
	mEnumeratedLevels->GetLevelsCount(&count);
	TGlucoseLevel *levels;
	mEnumeratedLevels->GetLevels(&levels);
	//invalid levels
	if (count < 2) {
		return S_FALSE;
	}

	//x is datetime, y is level
#if !defined(PARALLEL_AMP)

	floattype *t = (floattype*)malloc(count * sizeof(floattype));
	floattype *m = (floattype*)malloc((count+1) * sizeof(floattype));
	//tail boundary values
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
		//slopes
		//m{i} = (y{i+1}-y{i})/(x{i+1}-x{i})
		m[i - 2] = (levels[i - 1].level - levels[i - 2].level) / (levels[i - 1].datetime - levels[i - 2].datetime);
		//std::cout << "wtf " << m[0] << "\n";
		//m[i] = (levels[i + 1].level - levels[i].level) / dx;
		//zero check
		if (m[i+1] == m[i] && m[i-1] == m[i-2]) {
			t[i] = 0.5 * (m[i] + m[i-1]);
		}
		else {
			//t{i} = (abs(m{i + 1} - m{i})*m{i - 1} + abs(m{i - 1} - m{i - 2})*m{i}) / (abs(m{i + 1} - m{i}) + abs(m{i - 1} - m{i - 2}))
			t[i] = (abs(m[i + 1] - m[i])*m[i - 1] + abs(m[i - 1] - m[i - 2])*m[i]) / (abs(m[i + 1] - m[i]) + abs(m[i - 1] - m[i - 2]));;
		}
		floattype dx = levels[i + 1].datetime - levels[i].datetime;
		//a{i} is y{i}, b{i} is t{i}, c{i} is (3 * m{i} - 2 * t{i} - t{i+1}) / (x{i+1} -x{i}) and d{i} is (t{i} + t{i+1} - 2*m{i}) / (x{i+1} -x{i})^2
		aParams.insert(std::pair<floattype, akima_params>(levels[i].datetime, akima_params(levels[i].level, t[i],
			(3 * m[i] - 2 * t[i] - t[i + 1]) / dx, (t[i] + t[i + 1] - 2 * m[i]) / (dx*dx))));
	}
	// head boundary values
	//m[0] = (levels[1].level - levels[0].level) / (levels[1].datetime - levels[0].datetime);
	//m[1] = (levels[2].level - levels[1].level) / (levels[2].datetime - levels[1].datetime);
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
		//std::cout << "nechapu tohle " << t[0] << " " << m[0] << " " << mm1 << " " << mm2<< " " << m[1]<< " " << m[2] << " " << m[3] <<std::endl;
		//std::cout << "Wtf asdasd " << t[0] << "abs " << abs(m[1] - m[0])*mm1 << "mka " << mm1 << "dv " << (levels[1].level - levels[0].level) / (levels[1].datetime - levels[0].datetime)
			//<< std::endl;
	}
	aParams.insert(std::pair<floattype, akima_params>(levels[0].datetime,
		akima_params(levels[0].level, t[0], (3 * m[0] - 2 * t[0] - t[1]) / (levels[1].datetime - levels[0].datetime),
		(t[0] + t[1] - 2 * m[0]) / ((levels[1].datetime - levels[0].datetime)*(levels[1].datetime - levels[0].datetime)))));



	//add last timestamp to check bounds in getLevels
	lastTime = levels[count - 1].datetime;
	free(m);
	free(t);
	//printf("aproximuju %d %f %f\n", count, cParams.begin()->first, lastTime);
	//return S_OK;
#endif 
#ifdef PARALLEL_AMP
	std::vector<floattype> levelsV;
	std::vector<floattype> datetimesV;
	akima_params *par = (akima_params*)malloc((count-1) * sizeof(akima_params));
	//prepare vectors
	for (size_t i = 0; i < count; i++) {
		levelsV.push_back(levels[i].level);
		datetimesV.push_back(levels[i].datetime);
	}
	// Create C++ AMP objects.  
	concurrency::array_view<floattype, 1> y(count, levelsV);
	concurrency::array_view<floattype, 1> x(count, datetimesV);
	concurrency::array_view<akima_params, 1> parP(count - 1, par);
	concurrency::extent<1> ext(count - 3);
	parP.discard_data();
	//inside params
	/*
	concurrency::parallel_for_each(
		// Define the compute domain, which is the set of threads that are created.  
		parP.extent,
		// Define the code to run on each thread on the accelerator.  
		[=](concurrency::index<1> idx) restrict(amp)
	{
		//m{i} = (y{i+1}-y{i})/(x{i+1}-x{i})
		floattype t, t1;
		const int i = idx[0] + 2;
		floattype dy = y[i + 1] - y[i];
		floattype dx = x[i + 1] - x[i];
		floattype m0 = (y[i - 1] - y[i - 2]) / (x[i - 1] - x[i - 2]);
		floattype m1 = (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
		floattype m2 = dy / dx;
		floattype m3 = (y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]);
		floattype m4 = (y[i + 3] - y[i + 2]) / (x[i + 3] - x[i + 2]);
		//t{i} = (abs(m{i + 1} - m{i})*m{i - 1} + abs(m{i - 1} - m{i - 2})*m{i}) / (abs(m{i + 1} - m{i}) + abs(m{i - 1} - m{i - 2}))
		if (m3 == m2 && m1 == m0) {
			t = 0.5 * (m2 + m1);
		}
		else {
			t = (concurrency::fast_math::fabs(m3 - m2)*m1 + concurrency::fast_math::fabs(m1 - m0)*m2) /
				(concurrency::fast_math::fabs(m3 - m2) + concurrency::fast_math::fabs(m1 - m0));
		}
		if (m4 == m3 && m2 == m1) {
			t1 = 0.5 * (m3 + m2);
		}
		else {
			t1 = (concurrency::fast_math::fabs(m4 - m3)*m2 + concurrency::fast_math::fabs(m2 - m1)*m3) /
				(concurrency::fast_math::fabs(m4 - m3) + concurrency::fast_math::fabs(m2 - m1)); 
		}
		//a{i} is y{i}, b{i} is t{i}, c{i} is (3 * m{i} - 2 * t{i} - t{i+1}) / (x{i+1} -x{i}) and d{i} is (t{i} + t{i+1} - 2*m{i}) / (x{i+1} -x{i})^2
		parP[i].a = y[i]; 
		parP[i].b = t;
		parP[i].c = (3 * m2 - 2 * t - t1) / dx;
		parP[i].d = (t + t1 - 2 * m2) / (dx*dx);	
	}
	);*/
	std::vector<floattype> vec_times(count), vec_lvls(count);
	for (size_t i = 0; i < count; i++) {
		vec_times[i] = levels[i].datetime;
		vec_lvls[i] = levels[i].level;
	}

//	concurrency::extent<1> ext(count - 5);
	const concurrency::array_view<const floattype, 1> times(count, vec_times), lvls(count, vec_lvls);
	//concurrency::array_view<floattype, 1> p1s(size, p1), p2s(size, p2), p3s(size, p3);
	//p1s.discard_data();
	//p2s.discard_data();
	//p3s.discard_data();
	concurrency::parallel_for_each(ext, [=](concurrency::index<1> idx) restrict(amp) {
		floattype yy, xx, ti, ti1, m0, m1, m2, m3, m4;
		const int i = idx[0] + 2; // offset
		m0 = (lvls[i - 1] - lvls[i - 2]) / (times[i - 1] - times[i - 2]);
		m1 = (lvls[i] - lvls[i - 1]) / (times[i] - times[i - 1]);
		m2 = (lvls[i + 1] - lvls[i]) / (times[i + 1] - times[i]);
		m3 = (lvls[i + 2] - lvls[i + 1]) / (times[i + 2] - times[i + 1]);
		m4 = (lvls[i + 3] - lvls[i + 2]) / (times[i + 3] - times[i + 2]);
		ti = get_t(m0, m1, m2, m3);
		ti1 = get_t(m1, m2, m3, m4);

		yy = lvls[i + 1] - lvls[i];
		xx = times[i + 1] - times[i];
		parP[i].a = y[i];
		parP[i].b = ti;
		parP[i].c = get_p2(xx, yy, ti, ti1);
		parP[i].d = get_p3(xx, yy, ti, ti1);
	});
	parP.synchronize();
	//tail boundary params
	floattype tb, tb1, tb2;
	//m{i} = (y{i+1}-y{i})/(x{i+1}-x{i})
	floattype dy0 = levels[count - 2].level - levels[count - 3].level;
	floattype dx0 = levels[count - 2].datetime - levels[count - 3].datetime;
	floattype dy1 = levels[count - 1].level - levels[count - 2].level;
	floattype dx1 = levels[count - 1].datetime - levels[count - 2].datetime;
	floattype m3 = dy1 / dx1;
	floattype m2 = dy0 / dx0;
	floattype m1 = (levels[count - 3].level - levels[count - 4].level) / (levels[count - 3].datetime - levels[count - 4].datetime);
	floattype m0 = (levels[count - 4].level - levels[count - 5].level) / (levels[count - 4].datetime - levels[count - 5].datetime);
	//m{i-1} = 2m{i-2} - m{i-3}
	floattype m4 = 2 * m3 - m2;
	//m{i} = 2m{i-1} - m{i-2}
	floattype m5 = 2 * m4 - m3;
	//t{i+1}
	//zero check
	if (m3 == m2 && m1 == m0) {
		tb = 0.5 * (m2 + m1);
	}
	else {
		tb = (abs(m3 - m2)*m1 + abs(m1 - m0)*m2) / (abs(m3 - m2) + abs(m1 - m0));
	}
	if (m4 == m3 && m2 == m1) {
		tb1 = 0.5 * (m3 + m2);
	}
	else {
		tb1 = (abs(m4 - m3)*m2 + abs(m2 - m1)*m3) / (abs(m4 - m3) + abs(m2 - m1));
	}

	if (m5 == m4 && m3 == m2) {
		tb2 = 0.5 * (m4 + m3);
	}
	else {
		tb2 = (abs(m5 - m4)*m3 + abs(m3 - m2)*m4) / (abs(m5 - m4) + abs(m3 - m2));
	}
	parP[count-2] = akima_params(levels[count - 2].level, tb1, (3 * m3 - 2 * tb1 - tb2) / dx1,
		(tb1 + tb2 - 2 * m3) / (dx1*dx1));


	parP[count - 3] = akima_params(levels[count - 3].level, tb, (3 * m2 - 2 * tb - tb1) / dx0,
		(tb + tb1 - 2 * m2) / (dx0*dx0));
	//head boundary params

	//m{i} = (y{i+1}-y{i})/(x{i+1}-x{i})
	dy1 = levels[2].level - levels[1].level;
	dx1 = levels[2].datetime - levels[1].datetime;
	dy0 = levels[1].level - levels[0].level;
	dx0 = levels[1].datetime - levels[0].datetime;	
	//m{i-1} = 2m{i-2} - m{i-3}
	//m{i} = 2m{i-1} - m{i-2}
	m5 = (levels[4].level - levels[3].level) / (levels[4].datetime - levels[3].datetime);
	m4 = (levels[3].level - levels[2].level) / (levels[3].datetime - levels[2].datetime);
	m3 = dy1 / dx1;
	m2 = dy0 / dx0;
	m1 = 2 * m2 - m3;
	m0 = 2 * m1 - m2;
	//zero check
	if (m3 == m2 && m1 == m0) {
		tb = 0.5 * (m2 + m1);
	}
	else {
		tb = (abs(m3 - m2)*m1 + abs(m1 - m0)*m2) / (abs(m3 - m2) + abs(m1 - m0));
	}
	if (m4 == m3 && m2 == m1) {
		tb1 = 0.5 * (m3 + m2);
	}
	else {
		tb1 = (abs(m4 - m3)*m2 + abs(m2 - m1)*m3) / (abs(m4 - m3) + abs(m2 - m1));
	}

	if (m5 == m4 && m3 == m2) {
		tb2 = 0.5 * (m4 + m3);
	}
	else {
		tb2 = (abs(m5 - m4)*m3 + abs(m3 - m2)*m4) / (abs(m5 - m4) + abs(m3 - m2));
	}
	parP[1] = akima_params(levels[1].level, tb1, (3 * m3 - 2 * tb1 - tb2) / dx1,
		(tb1 + tb2 - 2 * m3) / (dx1*dx1));


	//std::cout << "nechapu tohle " << tb << " " << m2 << " " << m1 << " " << m0 << " " << m3 << " " << m4 << " " << m5 << std::endl;
	parP[0] = akima_params(levels[0].level, tb, (3 * m2 - 2 * tb - tb1) / dx0,
		(tb + tb1 - 2 * m2) / (dx0*dx0));
	//fill map
	for (int i = 0; i < count - 1; i++) {
		aParams.insert(std::pair<floattype, akima_params>(levels[i].datetime, parP[i])); 
		printf("jedu %f %f %f %f\n", parP[i].a, parP[i].b, parP[i].c, parP[i].d);// aParams[levels[i].datetime].b);// mC[i]);
	}
	free(par);
#endif
	//add last timestamp to check bounds in getLevels
	lastTime = levels[count - 1].datetime;
	return S_OK;
}
HRESULT  CAkimaSpline::GetLevels(floattype desiredtime, floattype stepping, size_t count, floattype *levels, size_t *filled, size_t derivationorder) {
	*filled = 0;
	/*Not aproximated yet*/
	if (aParams.empty()) {
		printf("err empty aparams\n");
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
			//workaround lower bound;
			if (it != aParams.begin() && it->first != currenttime)
			{
				it--;
			}
			/*workround for last value*/
			if (it == aParams.end()) {
				it--;
			}
			floattype time = currenttime - it->first;
			//y = a{i} + b{i}*(x-x{i})+ c{i}*(x-x{i})^2+d{i}*(x-x{i})^3
			res = it->second.a + it->second.b*time + it->second.c*time*time + it->second.d*time*time*time;
			levels[(*filled)] = res;
			(*filled)++;
			currenttime += stepping;
			count--;
			printf("tu jsem %f %f %f %f %f %f\n", res, it->second.a, it->second.b, it->second.c, it->second.d, time);
		}
	}
	else if (derivationorder == 1) {
		while (currenttime <= lastTime && count > 0) {
			it = aParams.lower_bound(currenttime);
			//workaround lower bound;
			if (it != aParams.begin() && it->first != currenttime)
			{
				it--;
			}
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
	//wrong derivative order
	else {
		return S_FALSE;
	}
	return S_OK;
}

CAkimaSpline::CAkimaSpline(IGlucoseLevels *levels) : CCommonApprox(levels) {
}

CAkimaSpline::~CAkimaSpline() {
}