#pragma once
#include "common\iface\ApproxIface.h"
#include "approx\src\GlucoseLevels.h"
#include "approx\src\CubicSpline.h"
#include "approx\src\QuadraticSpline.h"
#include "approx\src\AkimaSpline.h"
#include "approx\src\Masker.h"
#include "approx\src\Statistics.h"
#include <stdio.h>
#include <sqlite3.h>
#include <vector>
#include <string>

HRESULT processSegment(CGlucoseLevels *segment);

/*Loads segment levels to struct and vector*/
static int getSegmentLevels(void *glucoseLevels, int argc, char **argv, char **azColName) {
	std::vector<TGlucoseLevel> *gl = (std::vector<TGlucoseLevel>*)glucoseLevels;
	TGlucoseLevel val = TGlucoseLevel();
	val.datetime = atof(argv[0]) - 2415020.5;// Conversion to 1.1.1900
	val.level = atof(argv[1]);
	//printf("hodnota: %f %s a cas %f %s\n", val.level, argv[1], val.datetime, argv[0]);
	gl->push_back(val);	
	return 0;

}
/*Loads segment ids to vector*/
static int getSegsIds(void *segmentIds, int argc, char **argv, char **azColName) {
	std::vector<std::string> *segs = (std::vector<std::string>*) segmentIds;
	segs->push_back(argv[0]);
	//printf("%s\n", argv[0]);
	return 0;

}
/*heler sql*/
static int getinfo(void *NotUsed, int argc, char **argv, char **azColName) {
	int i;
    for (i = 0; i<argc; i++) {
		     printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
		
	}
	printf("\n");
	return 0;
}
int wrapper();

int main(int argc, void *arrgv[]) {
	wrapper();

	return 0;
}


int wrapper() {
	/*Databse handling*/
	sqlite3 *db;
	char *zErrMsg = 0;
	int rc;
	/*TODO dynamic loading from file*/
	rc = sqlite3_open("C:/Users/Petr/Dropbox/Skola/PPR/semestralka/ppr2016/data/direcnet.sqlite", &db);
	if (rc) {
		fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));


		return(0);
	}
	else {
		//fprintf(stderr, "Opened database successfully\n");
	}

	std::vector<std::string> segmentIds;
	sqlite3_exec(db, "SELECT id FROM timesegment", getSegsIds, &segmentIds, &zErrMsg);
	//std::vector<CGlucoseLevels *> segments;
	for (int i = 0; i < segmentIds.size(); i++) {
		std::vector<TGlucoseLevel> levels = std::vector<TGlucoseLevel>();
		//vytvoris si isntanci CGlucoseLevel.Pres SetLevelCount reknes, kolik chces naalokkovat mista pro data.Pak si vytvoris pointer TGlucoseLevel * a ten prdnes do GetLevels.Ta metoda ti preda naalokovane misto na ten ukazatel a pak uz to jen jen preklopis tu pamet s datama na ten ukazatel pres memcpy
		std::string query = "SELECT  julianday(measuredat), ist FROM measuredvalue WHERE segmentid = " + segmentIds[i] + " AND ist IS NOT NULL";
		sqlite3_exec(db, query.c_str(), getSegmentLevels, &levels, &zErrMsg);
		if (levels.size() > 0) {
			CGlucoseLevels *segment = new CGlucoseLevels();
			segment->AddRef();
			segment->SetLevelsCount(levels.size());
			TGlucoseLevel *ptr;
			segment->GetLevels(&ptr);
			for (size_t i = 0; i < levels.size(); i++) {
				//printf("test %f\n", levels[i].level);
			}
			memcpy(ptr, levels.data(), levels.size() * sizeof(TGlucoseLevel));
			size_t size;
			segment->GetLevelsCount(&size);
			printf("cajk %d %d\n", levels.size(), size);
			//segments.push_back(segment);



			printf("segment_id: %s\n", segmentIds[i].c_str());
			processSegment(segment);

			segment->Release();

			printf("cajk \n");
			//printf("pocet %d %s \n", levels->size(), segmentIds[i].c_str());
		}
	}
/*
	printf("cajk \n");
	size_t size;
	segments[segments.size() - 1]->GetLevelsCount(&size);

	printf("cajk %d\n", size);
	TGlucoseLevel *ptr;
	segments[segments.size() - 1]->GetLevels(&ptr);
	for (size_t i = 0; i < segments.size(); i++) {
		size_t size;
		segments[i]->GetLevelsCount(&size);
		//printf("test %d\n", size);
	}

	CMasker *masker = new CMasker(segments[0]);
	size_t ss;
	CGlucoseLevels *test = new CGlucoseLevels();
	masker->GetLevels(1,(IGlucoseLevels **) &test);
	test->GetLevelsCount(&ss);
	printf("ssacek %d\n", ss);
	CCubicSpline *cubic = new CCubicSpline(test);*/
	/*for (size_t i = 255; i > 0; i--) {
		TApproximationParams cubicParams = {
			0,
			{	//mask
				i
			}
		};
		cubic->Approximate(&cubicParams);
	}*/
	/*TApproximationParams cubicParams = {
		0,
		{	//mask
			1
		}
	};
	cubic->Approximate(&cubicParams);*/
	//TApproximationParams cubicParams = TApproximationParams();
	//cubic->Approximate(&cubicParams);
	/*size_t filled;
	floattype *levels = (floattype*)malloc(500 * sizeof(floattype));
	cubic->GetLevels(0, 0.01, 50, levels, &filled, 0);
	free(levels);
	for (size_t i = 0; i < segments.size(); i++) {
		segments[i]->Release();
	}
	delete cubic;*/
	//sqlite3_exec(db, "SELECT * FROM measuredvalue WHERE segmentid = 2", getinfo, 0, &zErrMsg);
	/*TGlucoseLevel *glucoseLevels;
	sqlite3_exec(db, "SELECT * FROM ", fillLevels, glucoseLevels, &zErrMsg);*/
	sqlite3_close(db);

	return 0;
}


HRESULT processSegment(CGlucoseLevels *segment) {
	CMasker *masker = new CMasker(segment);
	CStatistics *stats = new CStatistics(segment);
	//size_t ss;
	/*process all mask for all methods*/
	CGlucoseLevels *gl;// = new CGlucoseLevels();
	//gl->AddRef();
	std::map<floattype, floattype> derivations;
	for (int mask = 255; mask > 0; mask--) {
		//test->GetLevelsCount(&ss);
		//printf("ssacek %d\n", ss);
		masker->GetLevels(mask, (IGlucoseLevels **)&gl);
		CCubicSpline *cubic = new CCubicSpline(gl);
		CQuadraticSpline *quadra = new CQuadraticSpline(gl);
		CAkimaSpline *akima = new CAkimaSpline(gl);
		/*dont need params*/
		cubic->Approximate(nullptr);
		quadra->Approximate(nullptr);
		akima->Approximate(nullptr);
		printf("\tmask: %#X\n", mask);
		if (mask == 255) {
			size_t size;
			gl->GetLevelsCount(&size);
			TGlucoseLevel *ptr;
			gl->GetLevels(&ptr);
			//printf("vel: %d\n", size);
			for (size_t i = 0; i < size; i++) {
				floattype value;
				size_t filled;
				if (cubic->GetLevels(ptr[i].datetime, 0, 1, &value, &filled, 1) == S_OK && filled == 1) {
					//printf("mask: %d %f %f\n", i, ptr[i].datetime, value);
					derivations.insert(std::pair<floattype, floattype>(ptr[i].datetime, value));					
				}
				else {
				}
			}
		}	


		//TODO getlevels zavola statisticky nastroj presne pro casy z masky 255
		floattype *levels = (floattype*)malloc(10 * sizeof(floattype));
		size_t filled;
		//printf("awtffasfaf\n");
		//cubic->GetLevels(36525, 0.1, 1, levels, &filled, 0);
		//printf("filled je %lu a hodnota %f\n", filled, levels[0]);
		//size_t size;
		//gl->GetLevelsCount(&size);
		//printf("mask: %lu\n", size);
		//stats->GetStats(cubic, derivations, mask);
		free(levels);
		delete cubic;
		delete quadra;
		delete akima;
		//delete gl;

	}
	//gl = nullptr;
	//delete gl;
	//printf("mrdka %d\n", gl->Release());
	delete masker;
	delete stats;

	return S_OK;
}