#pragma once

#include "CommonApprox.h"
#include "GlucoseLevels.h"

class CMasker{
protected:
	std::vector<IGlucoseLevels *> maskedLevels;
	std::vector<std::vector<floattype>> maskedArrLevels;
	std::vector<std::vector<floattype>> maskedArrDatetimes;
public:
	CMasker(IGlucoseLevels *levels);
	~CMasker();
	HRESULT  IfaceCalling GetLevels(int mask, IGlucoseLevels**levels);
};
