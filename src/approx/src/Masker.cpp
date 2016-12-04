#include "Masker.h"

#include <stdio.h>
CMasker::CMasker(IGlucoseLevels *levels) {
	if (levels != NULL) {
		size_t size;
		levels->GetLevelsCount(&size);
		TGlucoseLevel *ptr;
		levels->GetLevels(&ptr);
		//printf("wtf  %lu \n", size);
		//Calculate masks, always add first and last element of segment
		for (int m = 1; m < 255; m++) {
			int y;
			std::vector<TGlucoseLevel> masked = std::vector<TGlucoseLevel>();
			//if (size > 7) {
			for (y = 0; y < size; y++) {
				int mod = y % 8;
				//printf("wtf %d %d %d\n", m >> (7 - mod), (m >> (7 - mod)) & 1, m);
				if ((m >> (7 - mod)) & 1 || y == 0 || y == (size - 1)) {
					masked.push_back(ptr[y]);
				}
		//creating masked glucose levels
			}
				CGlucoseLevels *mask = new CGlucoseLevels();
				mask->SetLevelsCount(masked.size());
				TGlucoseLevel *ptr2;
				mask->GetLevels(&ptr2);
				memcpy(ptr2, masked.data(), masked.size() * sizeof(TGlucoseLevel));
				//printf("cotoje %d %lu %f %f\n", m, masked.size(), masked[0].datetime, masked[masked.size() - 1].datetime);
				mask->AddRef();
				maskedLevels.push_back(mask);
				//delete(ptr2);			
		}
		//add 255 default mask
		levels->AddRef();
		maskedLevels.push_back(levels);
	}
	//free(ptr);
}
HRESULT CMasker::GetLevels(int mask, IGlucoseLevels **levels){
	if ((mask < 1 &&  mask > 255) || maskedLevels.size() < 255) {
		fprintf(stderr, "Error getting levels from masker. Invalid mask or constructor masking failed.");
		return S_FALSE;
	}
	size_t ss; 
	maskedLevels[mask - 1]->GetLevelsCount(&ss);
	//printf("ajaja %d %d\n", mask, ss);
	//maskedLevels[mask - 1]->AddRef();
	*levels = maskedLevels[mask - 1];
	(*levels)->AddRef();
	return S_OK;
}

CMasker::~CMasker() {
	for (int i = 0; i < maskedLevels.size(); i++) {
		if (maskedLevels[i] != NULL) maskedLevels[i]->Release();
	}
}