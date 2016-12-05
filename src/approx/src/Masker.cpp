#include "Masker.h"
#include <iostream>

CMasker::CMasker(IGlucoseLevels *levels) {
	if (levels != NULL) {
		size_t size;
		levels->GetLevelsCount(&size);
		TGlucoseLevel *ptr;
		levels->GetLevels(&ptr);
		//Calculate masks, always add first and last element of segment
		for (int m = 1; m < 255; m++) {
			int y;
			std::vector<TGlucoseLevel> masked = std::vector<TGlucoseLevel>();
			//if (size > 7) {
			for (y = 0; y < size; y++) {
				int mod = y % 8;
				if ((m >> (7 - mod)) & 1 || y == 0 || y == (size - 1)) {
					masked.push_back(ptr[y]);
				}
			}
			//creating masked glucose levels
			CGlucoseLevels *mask = new CGlucoseLevels();
			mask->SetLevelsCount(masked.size());
			TGlucoseLevel *ptr2;
			mask->GetLevels(&ptr2);
			memcpy(ptr2, masked.data(), masked.size() * sizeof(TGlucoseLevel));
			mask->AddRef();
			maskedLevels.push_back(mask);	
		}
		//add 255 default mask
		levels->AddRef();
		maskedLevels.push_back(levels);
	}
}
HRESULT CMasker::GetLevels(int mask, IGlucoseLevels **levels){
	if ((mask < 1 &&  mask > 255) || maskedLevels.size() < 255) {
		std::cerr << "Error getting levels from masker. Invalid mask or constructor masking failed.";
		return S_FALSE;
	}
	size_t ss; 
	maskedLevels[mask - 1]->GetLevelsCount(&ss);
	*levels = maskedLevels[mask - 1];
	(*levels)->AddRef();
	return S_OK;
}

CMasker::~CMasker() {
	for (int i = 0; i < maskedLevels.size(); i++) {
		if (maskedLevels[i] != NULL) maskedLevels[i]->Release();
	}
}