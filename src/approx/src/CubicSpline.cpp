#include "CubicSpline.h"


HRESULT CCubicSpline::Approximate(TApproximationParams *params) {
	printf("aproximuju %d", params->mask);
	return S_OK;
}
HRESULT  CCubicSpline::GetLevels(floattype desiredtime, floattype stepping, size_t count, floattype *levels, size_t *filled, size_t derivationorder) {
	printf("vracim hodnoty");
	return S_OK;
}

CCubicSpline::CCubicSpline(IGlucoseLevels *levels) : CCommonApprox(levels) {
}

CCubicSpline::~CCubicSpline() {
}
