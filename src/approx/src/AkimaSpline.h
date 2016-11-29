#pragma once

#include "CommonApprox.h"
#include "GlucoseLevels.h"
#include "..\..\common\rtl\referencedImpl.h"
#include <map>

#pragma warning( push )
#pragma warning( disable : 4250 ) // C4250 - 'class1' : inherits 'class2::member' via dominance

typedef struct akima_params{
	floattype a;
	floattype b;
	floattype c;
	floattype d;
	akima_params(floattype a, floattype b, floattype c, floattype d) : a(a), b(b), c(c), d(d) {};
	akima_params() {};
} akima_params;


class CAkimaSpline : public CCommonApprox, public virtual CReferenced {
protected:
	std::map<floattype, akima_params> aParams;
	floattype lastTime;
public:
	CAkimaSpline(IGlucoseLevels *levels);
	virtual ~CAkimaSpline();
	//dctor has to be virtual, even if it is empty, due to the inheritance by dominance	
	HRESULT Approximate(TApproximationParams * params);
	HRESULT GetLevels(floattype desiredtime, floattype stepping, size_t count, floattype * levels, size_t * filled, size_t derivationorder);
};

#pragma warning( pop )