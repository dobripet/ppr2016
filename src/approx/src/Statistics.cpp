#include "Statistics.h"
#include "CubicSpline.h"
#include <iostream>
#include <fstream>

#include <stdio.h>
CStatistics::CStatistics(IGlucoseLevels *levels) : levels(levels) {
	if (levels != NULL) levels->AddRef();
}
HRESULT CStatistics::GetStats(IApproximatedGlucoseLevels *method){
	if (levels == NULL || method == NULL) {
		return S_FALSE;
	}
	std::vector<floattype> absError;
	std::vector<floattype> relError;
	size_t size;
	levels->GetLevelsCount(&size);
	TGlucoseLevel *ptr;
	levels->GetLevels(&ptr);
	std::ofstream file;
	file.open("C:/Users/Petr/Dropbox/Skola/PPR/semestralka/ppr2016/data/output.csv", std::ios::trunc);
	file.precision(5);
	file << "date" << ", " << "calc" << ", " << "real" << "\n";
	for (size_t i = 0; i < size; i++) {
		floattype value;
		size_t filled;
		if (method->GetLevels(ptr[i].datetime, 0, 1, &value, &filled, 0) == S_OK && filled == 1) {
			floattype err = abs(value - ptr[i].level);
			absError.push_back(err);
			relError.push_back(err/ptr[i].level);
			file << (ptr[i].datetime - 36500.) << ", " << value << ", " << ptr[i].level << "\n";
		}

	}
	file.close();
	return S_OK;
}
CStatistics::~CStatistics() {
	if (levels != NULL) levels->Release();
}