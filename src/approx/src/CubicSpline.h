#pragma once

#include "CommonApprox.h"
#include "GlucoseLevels.h"
#include "..\..\common\rtl\referencedImpl.h"

#pragma warning( push )
#pragma warning( disable : 4250 ) // C4250 - 'class1' : inherits 'class2::member' via dominance

class CCubicSpline : public CCommonApprox, public virtual CReferenced {
public:
	CCubicSpline(IGlucoseLevels *levels);
	virtual ~CCubicSpline();
	//dctor has to be virtual, even if it is empty, due to the inheritance by dominance	
	HRESULT Approximate(TApproximationParams * params);
	HRESULT GetLevels(floattype desiredtime, floattype stepping, size_t count, floattype * levels, size_t * filled, size_t derivationorder);
};

#pragma warning( pop )