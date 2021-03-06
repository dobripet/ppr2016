#pragma once

#include "CommonApprox.h"
#include "GlucoseLevels.h"
#include "..\..\common\rtl\referencedImpl.h"
#include <map>

#pragma warning( push )
#pragma warning( disable : 4250 ) // C4250 - 'class1' : inherits 'class2::member' via dominance

//struct to contain cubic params
typedef struct cubic_params{
	floattype a;
	floattype b;
	floattype c;
	floattype d;
	cubic_params(floattype a, floattype b, floattype c, floattype d) : a(a), b(b), c(c), d(d) {};
	cubic_params() {};
} cubic_params;


class CCubicSpline : public CCommonApprox, public virtual CReferenced {
protected:
	//mapping datetime to cubic params
	std::map<floattype,cubic_params> cParams;
	//lastpoint datetime
	floattype lastTime;
public:
	CCubicSpline(IGlucoseLevels *levels);
	virtual ~CCubicSpline();
	//dctor has to be virtual, even if it is empty, due to the inheritance by dominance	
	HRESULT Approximate(TApproximationParams * params);
	HRESULT GetLevels(floattype desiredtime, floattype stepping, size_t count, floattype * levels, size_t * filled, size_t derivationorder);
};

#pragma warning( pop )