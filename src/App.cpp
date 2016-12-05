#pragma once
#include "common\iface\ApproxIface.h"
#include "approx\src\GlucoseLevels.h"
#include "approx\src\CubicSpline.h"
#include "approx\src\QuadraticSpline.h"
#include "approx\src\AkimaSpline.h"
#include "approx\src\CatmullRomSpline.h"
#include "approx\src\Masker.h"
#include "approx\src\Statistics.h"
#include <iostream>
#include <sstream>
#include <sqlite3.h>
#include <vector>
#include <string>
#include <amp.h> 
#include <tbb/tbb.h>  
HRESULT processSegment(int method, CGlucoseLevels *segment);
HRESULT processMask(int method, int mask, CMasker *masker, std::vector<CCommonApprox*> *instances);
HRESULT processStats(CMasker *masker, CStatistics *stats, std::vector<CCommonApprox*> *instances, std::vector<std::string> *results);
HRESULT processStatsMask(int mask, CStatistics *stats, std::map<floattype, floattype> derivations, std::vector<CCommonApprox*> *instances, std::vector<std::string> *results);
//Loads segment levels to struct and vector
static int getSegmentLevels(void *glucoseLevels, int argc, char **argv, char **azColName) {
	std::vector<TGlucoseLevel> *gl = (std::vector<TGlucoseLevel>*)glucoseLevels;
	TGlucoseLevel val = TGlucoseLevel();
	val.datetime = atof(argv[0]) - 2415020.5;// Conversion to 1.1.1900
	val.level = atof(argv[1]);
	gl->push_back(val);	
	return 0;

}
//Loads segment ids to vector
static int getSegsIds(void *segmentIds, int argc, char **argv, char **azColName) {
	std::vector<std::string> *segs = (std::vector<std::string>*) segmentIds;
	segs->push_back(argv[0]);
	return 0;

}

int main(int argc, void *argv[]) {
	//time measurment
	CTimer *ty = new CTimer();
	ty->start();
	int method;
	if (argc != 3) {
		std::cerr << "Invalid parameters! Method flags Akima spline: 1, Cubic spline: 2 Catmull-Rom Spline 3 Quadratic spline: 4. Usage: ppr2016 PATH_TO_DB_FILE METHOD";
		return 1;
	}
	method = atoi((char*)argv[2]);
	if (method < 1 || method > 4) {
		std::cerr << "Invalid parameters! Method flags Akima spline: 1, Cubic spline: 2 Catmull-Rom Spline 3 Quadratic spline: 4. Usage: ppr2016 PATH_TO_DB_FILE METHOD";
		return 1;
	}
#ifdef PARALLEL_AMP
	//AMP init
	CTimer *t = new CTimer();
	t->start();
	std::wcout << "Default accelerator is now "
		<< concurrency::accelerator(concurrency::accelerator::default_accelerator).description << "\n";
	t->stop();
	delete t;
#endif
	//Databse handling
	sqlite3 *db;
	char *zErrMsg = 0;
	if (sqlite3_open((char*)argv[1], &db)) {
		std::cerr << "Can't open database: " << sqlite3_errmsg(db)<< "\n";
		return 1;
	}
	std::vector<std::string> segmentIds;
	sqlite3_exec(db, "SELECT id FROM timesegment", getSegsIds, &segmentIds, &zErrMsg);
	//print method
	switch (method) {
	case 1: {
		std::cout << "Akima spline\n";
		break;
	}
	case 2: {
		std::cout << "Cubic spline\n";
		break;
	}
	case 3: {
		std::cout << "Cetripetal Catmull-Rom spline\n";
		break;
	}
	case 4: {
		std::cout << "Quadratic spline\n";
		break;
	}
	}
	for (int i = 0; i < segmentIds.size(); i++) {
		std::vector<TGlucoseLevel> levels = std::vector<TGlucoseLevel>();
		std::string query = "SELECT  julianday(measuredat), ist FROM measuredvalue WHERE segmentid = " + segmentIds[i] + " AND ist IS NOT NULL";
		sqlite3_exec(db, query.c_str(), getSegmentLevels, &levels, &zErrMsg);
		//process segment with levels
		if (levels.size() > 0) {
			CGlucoseLevels *segment = new CGlucoseLevels();
			segment->AddRef();
			segment->SetLevelsCount(levels.size());
			TGlucoseLevel *ptr;
			segment->GetLevels(&ptr);
			memcpy(ptr, levels.data(), levels.size() * sizeof(TGlucoseLevel));
			//print segment id
			std::cout << "\tsegment_id: " << segmentIds[i].c_str() << "\n";
			//process segment values for method
			processSegment(method, segment);
			segment->Release();
		}
	}
	sqlite3_close(db);

	//stop timer
	ty->stop();
	delete ty;

	return 0;
}

HRESULT processSegment(int method, CGlucoseLevels *segment) {
	CMasker *masker = new CMasker(segment);
	CStatistics *stats = new CStatistics(segment);
	std::vector<CCommonApprox*> instances(255);
	//process masks
#if !defined(PARALLEL_TBB)
	for (int mask = 255; mask > 0; mask--) {
		processMask(method, mask, masker, &instances);
	}
#endif
#if defined(PARALLEL_TBB)	
	tbb::parallel_for(1, 256, [&](int mask) {
		processMask(method, mask, masker, &instances);
	});
#endif
	//prcess stats and print results
	std::vector<std::string> results(255);
	processStats(masker, stats, &instances, &results);
	for (int i = 255; i > 0; i--) {
		std::cout << results[i - 1];
	}
	delete masker;
	delete stats;
	return S_OK;
}


//process mask for given method
HRESULT processMask(int method, int mask, CMasker *masker, std::vector<CCommonApprox*> *instances) {
	CGlucoseLevels *gl;
	masker->GetLevels(mask, (IGlucoseLevels **)&gl);
	switch (method) {
	case 1: {
		CAkimaSpline *akima = new CAkimaSpline(gl);
		akima->Approximate(nullptr);
		(*instances)[mask - 1] = akima;
		break;
	}
	case 2: {
		CCubicSpline *cubic = new CCubicSpline(gl);
		cubic->Approximate(nullptr); 
		(*instances)[mask - 1] = cubic;
		break;
	}
	case 3: {
		CCatmullRomSpline *cat = new CCatmullRomSpline(gl);
		cat->Approximate(nullptr);
		(*instances)[mask - 1] = cat;
		break;
	}
	case 4: {
		CQuadraticSpline *quadra = new CQuadraticSpline(gl);
		quadra->Approximate(nullptr);
		(*instances)[mask - 1] = quadra;
		break;
	}
	}
	gl->Release();
	return S_OK;
}
//get all statistics for segment
HRESULT processStats(CMasker *masker, CStatistics *stats, std::vector<CCommonApprox*> *instances, std::vector<std::string> *results) {
	//get derivations for mask 255
	CGlucoseLevels *gl;
	masker->GetLevels(255, (IGlucoseLevels **)&gl);
	std::map<floattype, floattype> derivations;
	size_t size;
	gl->GetLevelsCount(&size);
	TGlucoseLevel *ptr;
	gl->GetLevels(&ptr);
	for (size_t i = 0; i < size; i++) {
		floattype value;
		size_t filled;
		if((*instances)[254]->GetLevels(ptr[i].datetime, 0, 1, &value, &filled, 1) == S_OK && filled == 1) {
			derivations.insert(std::pair<floattype, floattype>(ptr[i].datetime, value));
		}
		else {
			return S_FALSE;
		}
	}
	gl->Release();
	//process stat for masks
#if !defined(PARALLEL_TBB)
	for (int mask = 255; mask > 0; mask--) {
		processStatsMask(mask, stats, derivations, instances, results);
	}
#endif

#if defined(PARALLEL_TBB)
	tbb::parallel_for(1, 256, [&](int mask) {
		processStatsMask(mask, stats, derivations, instances, results);
	});
#endif
	return S_OK;
}
//process stats for individual mask and result add to result vector
HRESULT processStatsMask(int mask, CStatistics *stats, std::map<floattype, floattype> derivations, std::vector<CCommonApprox*> *instances, std::vector<std::string> *results) {
	std::stringstream ss;
	ss.clear();
	char buf[16];
	sprintf(buf, "\t\tmask: %#X\n", mask); 
	ss << buf;
	std::string res;
	stats->GetStats((*instances)[mask - 1], derivations, mask, &res);
	delete (*instances)[mask - 1];
	ss << res;
	(*results)[mask - 1] = ss.str();
	return S_OK;
}