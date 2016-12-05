#pragma once

#include "CommonApprox.h"
#include "GlucoseLevels.h"
#include "..\..\common\rtl\referencedImpl.h"
#include <map>

#pragma warning( push )
#pragma warning( disable : 4250 ) // C4250 - 'class1' : inherits 'class2::member' via dominance

//struct to contain quadratic params
typedef struct quadratic_params{
	floattype a;
	floattype b;
	floattype c;
	quadratic_params(floattype a, floattype b, floattype c) : a(a), b(b), c(c) {};
	quadratic_params() {};
} quadratic_params;


class CQuadraticSpline : public CCommonApprox, public virtual CReferenced {
protected:
	//mapping datetime to quadratic params
	std::map<floattype, quadratic_params> qParams;
	//lastpoint datetime
	floattype lastTime;
public:
	CQuadraticSpline(IGlucoseLevels *levels);
	virtual ~CQuadraticSpline();
	//dctor has to be virtual, even if it is empty, due to the inheritance by dominance	
	HRESULT Approximate(TApproximationParams * params);
	HRESULT GetLevels(floattype desiredtime, floattype stepping, size_t count, floattype * levels, size_t * filled, size_t derivationorder);
};

#pragma warning( pop )
