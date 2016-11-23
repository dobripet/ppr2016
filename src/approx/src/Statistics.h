#pragma once

#include "CommonApprox.h"
#include "GlucoseLevels.h"

class CStatistics{
protected:
	IGlucoseLevels *levels;
public:
	CStatistics(IGlucoseLevels *levels);
	~CStatistics();

	HRESULT IfaceCalling  CStatistics::GetStats(IApproximatedGlucoseLevels *method);
	//HRESULT  IfaceCalling GetLevels(int mask, IGlucoseLevels**levels);
};
