#pragma once

#include "CommonApprox.h"
#include "GlucoseLevels.h"

class CMasker{
protected:
	std::vector<IGlucoseLevels *> maskedLevels;
public:
	CMasker(IGlucoseLevels *levels);
	~CMasker();
	HRESULT  IfaceCalling GetLevels(int mask, IGlucoseLevels**levels);
};
