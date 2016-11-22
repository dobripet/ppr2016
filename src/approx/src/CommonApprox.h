#pragma once

#include "../../common/iface/ApproxIface.h"
#include "..\..\common\rtl\hresult.h"
#include "..\..\common\rtl\referencedImpl.h"


#pragma warning( push )
#pragma warning( disable : 4250 ) // C4250 - 'class1' : inherits 'class2::member' via dominance

class CCommonApprox : public IApproximatedGlucoseLevels, public virtual CReferenced {
protected:
	IGlucoseLevels *mEnumeratedLevels;	
public:
	CCommonApprox(IGlucoseLevels *levels);
	virtual ~CCommonApprox();
	//dctor has to be virtual, even if it is empty, due to the inheritance by dominance	
};

#pragma warning( pop )