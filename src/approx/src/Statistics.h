#pragma once

#include "CommonApprox.h"
#include "GlucoseLevels.h"
#include <map>

class CStatistics{
protected:
	IGlucoseLevels *levels;
	HRESULT getQuartils(floattype *min, floattype *firstQ, floattype *median, floattype *thirdQ, floattype *max, std::vector<floattype> error); 
	HRESULT getStatsForMask(floattype *absMean, floattype *absMin, floattype *absFirstQ, floattype *absMedian, floattype *absThirdQ, floattype *absMax,
		floattype *absStandardDeviation, floattype *relMean, floattype *relMin, floattype *relFirstQ, floattype *relMedian, floattype *relThirdQ,
		floattype *relMax, floattype *relStandardDeviation, IApproximatedGlucoseLevels *method, unsigned int mask, int derivation, std::map<floattype, floattype> derivations);
	HRESULT printStats(floattype absMean, floattype absMin, floattype absFirstQ, floattype absMedian, floattype absThirdQ, floattype absMax,
		floattype absStandardDeviation, floattype relMean, floattype relMin, floattype relFirstQ, floattype relMedian, floattype relThirdQ,
		floattype relMax, floattype relStandardDeviation, unsigned int mask);
public:
	CStatistics(IGlucoseLevels *levels);
	~CStatistics();

	HRESULT IfaceCalling  CStatistics::GetStats(IApproximatedGlucoseLevels *method, std::map<floattype, floattype> derivations, unsigned int mask);
	//HRESULT  IfaceCalling GetLevels(int mask, IGlucoseLevels**levels);
};
