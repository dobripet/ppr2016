#include "Statistics.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm> 

#include <stdio.h>
HRESULT CStatistics::printStats(floattype mean, floattype absMin, floattype absFirstQ, floattype absMedian, floattype absThirdQ, floattype absMax,
	floattype absStandardDeviation, floattype relMean, floattype relMin, floattype relFirstQ, floattype relMedian, floattype relThirdQ,
	floattype relMax, floattype relStandardDeviation, unsigned int mask) {
	return S_OK;
}
HRESULT CStatistics::getStatsForMask(floattype *absMean, floattype *absMin, floattype *absFirstQ, floattype *absMedian, floattype *absThirdQ, floattype *absMax,
	floattype *absStandardDeviation, floattype *relMean, floattype *relMin, floattype *relFirstQ, floattype *relMedian, floattype *relThirdQ,
	floattype *relMax, floattype *relStandardDeviation, IApproximatedGlucoseLevels *method, unsigned int mask, int derivation, std::map<floattype, floattype> derivations) {
	floattype absSum = 0.;
	floattype absSumsq = 0.;
	floattype relSum = 0.;
	floattype relSumsq = 0.;
	std::vector<floattype> absError;
	std::vector<floattype> relError;
	size_t size;
	levels->GetLevelsCount(&size);
	TGlucoseLevel *ptr;
	levels->GetLevels(&ptr);
	std::ofstream file;
	/*
	file.open("C:/Users/Petr/Dropbox/Skola/PPR/semestralka/ppr2016/data/output.csv", std::ios::trunc);
	file.precision(5);
	file << "date" << ", " << "calc" << ", " << "real" << "\n";
	*/
	for (size_t i = 0; i < size; i++) {
		int mod = i % 8;
		//get value for mask
		if ((mask >> (7 - mod)) & 1) {
			floattype value;
			size_t filled;
			if (derivation == 0) {
				if (method->GetLevels(ptr[i].datetime, 0, 1, &value, &filled, 0) == S_OK && filled == 1) {
					floattype err = abs(value - ptr[i].level);
					floattype relerr = err / ptr[i].level;
					absError.push_back(err);
					relError.push_back(relerr);
					//file << (ptr[i].datetime - ptr[0].datetime) << ", " << value << ", " << ptr[i].level << "\n";
					absSum += err;
					absSumsq += err * err;
					relSum += relerr;
					relSumsq += relerr * relerr;
				}
			}
			//first order derivation
			else if (derivation == 1) {
				if (method->GetLevels(ptr[i].datetime, 0, 1, &value, &filled, 1) == S_OK && filled == 1) {
					floattype err = abs(value - derivations[ptr[i].datetime]);
					absError.push_back(err);
					file << (ptr[i].datetime - ptr[0].datetime) << ", " << value << ", " << derivations[ptr[i].datetime] << "\n";
					absSum += err;
					absSumsq += err * err;
				}
			}

		}

	}
	//too few values for statistics
	if (size < 5) {
		return S_FALSE;
	}
	//absolute
	getQuartils(absMin, absFirstQ, absMedian, absThirdQ, absMax, absError);
	*absMean = absSum / size;
	*absStandardDeviation = sqrt((absSumsq - (((*absMean)*(*absMean)) / size)) / size);//sqrt from variance
	if (derivation == 0) {
		//relative
		getQuartils(relMin, relFirstQ, relMedian, relThirdQ, relMax, relError);
		*relMean = relSum / size;
		*relStandardDeviation = sqrt((relSumsq - (((*relMean)*(*relMean)) / size)) / size);//sqrt from variance
	}
	//file.close();
	return S_OK;
}
HRESULT CStatistics::getQuartils(floattype *min, floattype *firstQ, floattype *median, floattype *thirdQ, floattype *max, std::vector<floattype> error) {
	size_t half = error.size() / 2;
	size_t quarter = error.size() / 4;
	size_t threeQuarters = half + quarter;
	*min = 0.;
	*firstQ = 0.;
	*median = 0.;
	*thirdQ = 0.;
	*max = 0.;
	if (error.size() > 0) {		
		if (error.size() & 1) {
			//odd size
			//median
			std::nth_element(error.begin(), error.begin() + half, error.end());
			*median = error[half];
			//first quartil
			std::nth_element(error.begin(), error.begin() + quarter, error.begin() + half);
			*firstQ = error[quarter];
			//third quartil
			std::nth_element(error.begin() + half, error.begin() + threeQuarters, error.end());
			*thirdQ = error[threeQuarters];

		}
		else {
			//even size
			//median
			std::nth_element(error.begin(), error.begin() + half, error.end());
			std::nth_element(error.begin(), error.begin() + half - 1, error.begin() + half);
			*median = (error[half] + error[half - 1]) / 2;
			//first quartil
			std::nth_element(error.begin(), error.begin() + quarter, error.begin() + half);
			std::nth_element(error.begin(), error.begin() + quarter - 1, error.begin() + quarter);
			*firstQ = (error[quarter] + error[quarter - 1]) / 2;
			//third quartil
			std::nth_element(error.begin() + half, error.begin() + threeQuarters, error.end());
			std::nth_element(error.begin() + half, error.begin() + threeQuarters - 1, error.begin() + threeQuarters);
			*thirdQ = (error[threeQuarters] + error[threeQuarters - 1]) / 2;
		}
		//min
		std::nth_element(error.begin(), error.begin(), error.begin() + quarter);
		*min = error[0];
		//max
		std::nth_element(error.begin() + threeQuarters, error.end() - 1, error.end());
		*max = error[error.size() - 1];
	}
	return S_OK;
}

CStatistics::CStatistics(IGlucoseLevels *levels) : levels(levels) {
	if (levels != NULL) levels->AddRef();
}
HRESULT CStatistics::GetStats(IApproximatedGlucoseLevels *method, std::map<floattype, floattype> derivations, unsigned int mask, std::string *result){
	if (levels == NULL || method == NULL) {
		return S_FALSE;
	}
	//for all concentrations
	floattype absMin = 0.;
	floattype absFirstQ = 0.;
	floattype absMedian = 0.;
	floattype absThirdQ = 0.;
	floattype absMax = 0.;
	floattype relMin = 0.;
	floattype relFirstQ = 0.;
	floattype relMedian = 0.;
	floattype relThirdQ = 0.;
	floattype relMax = 0.;
	floattype absStandardDeviation = 0.;
	floattype relStandardDeviation = 0.;
	floattype absMean = 0.;
	floattype relMean = 0.;
	//get stats for all points
	char buf[128];
	std::stringstream ss;
	ss << "\t\t\tall concentrations\n";
	getStatsForMask(&absMean, &absMin, &absFirstQ, &absMedian, &absThirdQ, &absMax, &absStandardDeviation, &relMean, &relMin,
		&relFirstQ, &relMedian, &relThirdQ, &relMax, &relStandardDeviation, method, 255, 0, derivations);
	sprintf(buf, "\t\t\t\tabsolute: %f, %f, %f, %f, %f, %f, %f\n", absMean, absMin, absFirstQ, absMedian, absThirdQ, absMax, absStandardDeviation);
	ss << buf;
	sprintf(buf, "\t\t\t\trelative: %f, %f, %f, %f, %f, %f, %f\n", relMean, relMin, relFirstQ, relMedian, relThirdQ, relMax, relStandardDeviation);
	ss << buf;
	/*//get stats for points in mask 0
	ss << "\t\t\tbit mask 0\n";
	getStatsForMask(&absMean, &absMin, &absFirstQ, &absMedian, &absThirdQ, &absMax, &absStandardDeviation, &relMean, &relMin,
		&relFirstQ, &relMedian, &relThirdQ, &relMax, &relStandardDeviation, method, ((~mask)&0xFF), 0, derivations);
	sprintf(buf, "\t\t\t\tabsolute: %f, %f, %f, %f, %f, %f, %f\n", absMean, absMin, absFirstQ, absMedian, absThirdQ, absMax, absStandardDeviation);
	ss << buf;
	sprintf(buf, "\t\t\t\trelative: %f, %f, %f, %f, %f, %f, %f\n", relMean, relMin, relFirstQ, relMedian, relThirdQ, relMax, relStandardDeviation);
	ss << buf;
	//get stats for points in mask 1
	ss << "\t\t\tbit mask 1\n";
	getStatsForMask(&absMean, &absMin, &absFirstQ, &absMedian, &absThirdQ, &absMax, &absStandardDeviation, &relMean, &relMin,
		&relFirstQ, &relMedian, &relThirdQ, &relMax, &relStandardDeviation, method, (mask), 0, derivations);
	sprintf(buf, "\t\t\t\tabsolute: %f, %f, %f, %f, %f, %f, %f\n", absMean, absMin, absFirstQ, absMedian, absThirdQ, absMax, absStandardDeviation);
	ss << buf;
	sprintf(buf, "\t\t\t\trelative: %f, %f, %f, %f, %f, %f, %f\n", relMean, relMin, relFirstQ, relMedian, relThirdQ, relMax, relStandardDeviation);
	ss << buf;

	*/
	/*ss << "\t\t\tfirst continuous derivative all concentrations\n";
	//get stats for all points derivation 1
	getStatsForMask(&absMean, &absMin, &absFirstQ, &absMedian, &absThirdQ, &absMax, &absStandardDeviation, &relMean, &relMin,
		&relFirstQ, &relMedian, &relThirdQ, &relMax, &relStandardDeviation, method, 255, 1, derivations);
	sprintf(buf, "\t\t\t\tabsolute: %f, %f, %f, %f, %f, %f, %f\n", absMean, absMin, absFirstQ, absMedian, absThirdQ, absMax, absStandardDeviation);
	ss << buf;*/
	/*ss << "\t\t\tfirst continuous derivative bit mask 0\n";
	//get stats for points in mask 0 derivation 1
	getStatsForMask(&absMean, &absMin, &absFirstQ, &absMedian, &absThirdQ, &absMax, &absStandardDeviation, &relMean, &relMin,
		&relFirstQ, &relMedian, &relThirdQ, &relMax, &relStandardDeviation, method, ((~mask) & 0xFF), 1, derivations);
	sprintf(buf, "\t\t\t\tabsolute: %f, %f, %f, %f, %f, %f, %f\n", absMean, absMin, absFirstQ, absMedian, absThirdQ, absMax, absStandardDeviation);
	ss << buf;
	ss <<"\t\t\tfirst continuous derivative bit mask 1\n";
	//get stats for points in mask 1 derivation 1
	getStatsForMask(&absMean, &absMin, &absFirstQ, &absMedian, &absThirdQ, &absMax, &absStandardDeviation, &relMean, &relMin,
		&relFirstQ, &relMedian, &relThirdQ, &relMax, &relStandardDeviation, method, (mask), 1, derivations);
	sprintf(buf, "\t\t\t\tabsolute: %f, %f, %f, %f, %f, %f, %f\n", absMean, absMin, absFirstQ, absMedian, absThirdQ, absMax, absStandardDeviation);
	ss << buf;*/
	(*result) = ss.str();
	return S_OK;
}
CStatistics::~CStatistics() {
	if (levels != NULL) levels->Release();
}
HRESULT CTimer::start() {
	FILETIME time;
	GetSystemTimeAsFileTime(&time);
	startTime.HighPart = time.dwHighDateTime;
	startTime.LowPart = time.dwLowDateTime;
	return S_OK;
}
HRESULT CTimer::stop() {
	FILETIME time;
	GetSystemTimeAsFileTime(&time);
	stopTime.HighPart = time.dwHighDateTime;
	stopTime.LowPart = time.dwLowDateTime;
	std::cout << "Timer measured: " << (floattype)((stopTime.QuadPart - startTime.QuadPart)) /10000 << "ms" << std::endl;
	return S_OK;
}