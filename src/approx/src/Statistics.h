#pragma once

#include "CommonApprox.h"
#include "GlucoseLevels.h"
#include <map>

class CStatistics{
protected:
	//levels for 255 mask
	IGlucoseLevels *levels;
	//calc quartils
	HRESULT getQuartils(floattype *min, floattype *firstQ, floattype *median, floattype *thirdQ, floattype *max, std::vector<floattype> error);
	//calc statistics for given mask
	HRESULT getStatsForMask(floattype *absMean, floattype *absMin, floattype *absFirstQ, floattype *absMedian, floattype *absThirdQ, floattype *absMax,
		floattype *absStandardDeviation, floattype *relMean, floattype *relMin, floattype *relFirstQ, floattype *relMedian, floattype *relThirdQ,
		floattype *relMax, floattype *relStandardDeviation, IApproximatedGlucoseLevels *method, unsigned int mask, int derivation, std::map<floattype, floattype> derivations);
	//prints given statistics
	HRESULT printStats(floattype absMean, floattype absMin, floattype absFirstQ, floattype absMedian, floattype absThirdQ, floattype absMax,
		floattype absStandardDeviation, floattype relMean, floattype relMin, floattype relFirstQ, floattype relMedian, floattype relThirdQ,
		floattype relMax, floattype relStandardDeviation, unsigned int mask);
public:
	CStatistics(IGlucoseLevels *levels);
	~CStatistics();
	//calls getLevels method on given method instance, mask and derivations(255mask for given segment), return statistic filled to result
	HRESULT IfaceCalling  CStatistics::GetStats(IApproximatedGlucoseLevels *method, std::map<floattype, floattype> derivations, unsigned int mask, std::string *result);
};
//timer class
class CTimer {
protected:
	ULARGE_INTEGER startTime, stopTime;
public:
	//starts timer
	HRESULT IfaceCalling  CTimer::start();
	//stops timer and prints measured time to stdout in ms
	HRESULT IfaceCalling  CTimer::stop();

};