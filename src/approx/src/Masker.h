#pragma once

#include "CommonApprox.h"
#include "GlucoseLevels.h"

class CMasker{
protected:
	//vector of all mask levels
	std::vector<IGlucoseLevels *> maskedLevels;
public:
	CMasker(IGlucoseLevels *levels);
	~CMasker();
	//returns levels for given mask
	HRESULT  IfaceCalling GetLevels(int mask, IGlucoseLevels**levels);
};
