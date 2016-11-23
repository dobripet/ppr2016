#include "Masker.h"

#include <stdio.h>
CMasker::CMasker(IGlucoseLevels *levels) {
	if (levels != NULL) {
		size_t size;
		levels->GetLevelsCount(&size);
		TGlucoseLevel *ptr;
		levels->GetLevels(&ptr);
		printf("wtf  %lu \n", size);
		/*Calculate masks, always add first and last element of segment*/
		for (int i = 1; i < 255; i++) {
			int y;
			std::vector<TGlucoseLevel> masked = std::vector<TGlucoseLevel>();
			if (size > 7) {
				for (y = 0; y < size; y += 8) {
					//1. bit of mask
					if (i & 128 || y == 0) {
						masked.push_back(ptr[y]);
					}
					//2. bit of mask
					if (i & 64) {
						masked.push_back(ptr[y + 1]);
					}
					//3. bit of mask
					if (i & 32) {
						masked.push_back(ptr[y + 2]);
					}
					//4. bit of mask
					if (i & 16) {
						masked.push_back(ptr[y + 3]);
					}
					//5. bit of mask
					if (i & 8) {
						masked.push_back(ptr[y + 4]);
					}
					//6. bit of mask
					if (i & 4) {
						masked.push_back(ptr[y + 5]);
					}
					//7. bit of mask
					if (i & 2) {
						masked.push_back(ptr[y + 6]);
					}
					//8. bit of mask
					if (i & 1 || y == (size - 8)) {
						masked.push_back(ptr[y + 7]);
					}
					if (y > size - 8) {
						break;
					}
				}
			}
			/*tail*/
			//1. bit of mask
			if (i & 128 && y < size || y == 0 || y == (size - 1)) {
				masked.push_back(ptr[y]);
			}
			//2. bit of mask
			if (i & 64 && y + 1 < size || y == (size - 1)) {
				masked.push_back(ptr[y + 1]);
			}
			//3. bit of mask
			if (i & 32 && y + 2 < size || y == (size - 1)) {
				masked.push_back(ptr[y + 2]);
			}
			//4. bit of mask
			if (i & 16 && y + 3 < size || y == (size - 1)) {
				masked.push_back(ptr[y + 3]);
			}
			//5. bit of mask
			if (i & 8 && y + 4 < size || y == (size - 1)) {
				masked.push_back(ptr[y + 4]);
			}
			//6. bit of mask
			if (i & 4 && y + 5 < size || y == (size - 1)) {
				masked.push_back(ptr[y + 5]);
			}
			//7. bit of mask
			if (i & 2 && y + 6 < size || y == (size - 1)) {
				masked.push_back(ptr[y + 6]);
			}
			//8. bit of mask
			if (i & 1 && y + 7 < size || y == (size - 1)) {
				masked.push_back(ptr[y + 7]);
			}

			/*creating masked glucose levels*/
			CGlucoseLevels *mask = new CGlucoseLevels();
			mask->SetLevelsCount(masked.size());
			TGlucoseLevel *ptr2;
			mask->GetLevels(&ptr2);
			memcpy(ptr2, masked.data(), masked.size() * sizeof(TGlucoseLevel));
			printf("cotoje  %lu %f %f\n", masked.size(), masked[0].datetime, masked[masked.size() - 1].datetime);
			mask->AddRef();
			maskedLevels.push_back(mask);
			//delete(ptr2);
		}
		/*add 255 default mask*/
		levels->AddRef();
		maskedLevels.push_back(levels);
	}
	//free(ptr);
}
HRESULT CMasker::GetLevels(int mask, IGlucoseLevels **levels){
	if ((mask < 1 &&  mask > 255) || maskedLevels.size() < 255) {
		return S_FALSE;
	}
	size_t ss; 
	maskedLevels[mask - 1]->GetLevelsCount(&ss);
	//printf("ajaja %d\n", ss);
	//maskedLevels[mask - 1]->AddRef();
	*levels = maskedLevels[mask - 1];
	//(*levels)->AddRef();
	return S_OK;
}
CMasker::~CMasker() {
	for (int i = 0; i < maskedLevels.size(); i++) {
		if (maskedLevels[i] != NULL) maskedLevels[i]->Release();
	}
}