#pragma once

#include "CommonApprox.h"
#include "GlucoseLevels.h"
#include "..\..\common\rtl\referencedImpl.h"
#include <map>

#pragma warning( push )
#pragma warning( disable : 4250 ) // C4250 - 'class1' : inherits 'class2::member' via dominance

//struct to contain catmull-rom params
typedef struct cr_params{
	floattype a;
	floattype b;
	floattype c;
	floattype d;
	cr_params(floattype a, floattype b, floattype c, floattype d) : a(a), b(b), c(c), d(d) {};
	cr_params() {};
} cr_params;


class CCatmullRomSpline : public CCommonApprox, public virtual CReferenced {
protected:
	//mapping datetime to catmull-rom params
	std::map<floattype, cr_params> crParams;
	//lastpoint datetime
	floattype lastTime;
public:
	CCatmullRomSpline(IGlucoseLevels *levels);
	virtual ~CCatmullRomSpline();
	//dctor has to be virtual, even if it is empty, due to the inheritance by dominance	
	HRESULT Approximate(TApproximationParams * params);
	HRESULT GetLevels(floattype desiredtime, floattype stepping, size_t count, floattype * levels, size_t * filled, size_t derivationorder);
};

#pragma warning( pop )