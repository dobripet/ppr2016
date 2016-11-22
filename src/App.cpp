#pragma once
#include "common\iface\ApproxIface.h"
#include "approx\src\GlucoseLevels.h"
#include "approx\src\CubicSpline.h"
#include <stdio.h>
#include <sqlite3.h>
#include <vector>
#include <string>
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
int main(int argc, void *arrgv[]) {
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
	std::vector<CGlucoseLevels *> segments;
	for (int i = 0; i < segmentIds.size(); i++) {
		std::vector<TGlucoseLevel> levels = std::vector<TGlucoseLevel>();
		//vytvoris si isntanci CGlucoseLevel.Pres SetLevelCount reknes, kolik chces naalokkovat mista pro data.Pak si vytvoris pointer TGlucoseLevel * a ten prdnes do GetLevels.Ta metoda ti preda naalokovane misto na ten ukazatel a pak uz to jen jen preklopis tu pamet s datama na ten ukazatel pres memcpy
		std::string query = "SELECT  julianday(measuredat), ist FROM measuredvalue WHERE segmentid = " + segmentIds[i] + " AND ist IS NOT NULL";
		sqlite3_exec(db, query.c_str(), getSegmentLevels, &levels, &zErrMsg);
		if (levels.size() > 0) {
			CGlucoseLevels *segment = new CGlucoseLevels();
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
			segments.push_back(segment);
			
			//printf("pocet %d %s \n", levels->size(), segmentIds[i].c_str());
		}
	}

	printf("cajk \n");
	size_t size;
	segments[segments.size() - 1]->GetLevelsCount(&size);

	printf("cajk %d\n", size);
	TGlucoseLevel *ptr;
	segments[segments.size()-1]->GetLevels(&ptr);
	for (size_t i = 0; i < segments.size(); i++) {
		size_t size; 
		segments[i]->GetLevelsCount(&size);
		//printf("test %d\n", size);
	}

	CCubicSpline *cubic = new CCubicSpline(segments[0]);
	for (int i = 255; i > 0; i--) {
		TApproximationParams cubicParams = {
			0,
			{	//mask
				i
			}
		};
		cubic->Approximate(&cubicParams);
	}
	//TApproximationParams cubicParams = TApproximationParams();
	//cubic->Approximate(&cubicParams);
	size_t filled;
	floattype *levels =(floattype*) malloc(500 * sizeof(floattype));
	cubic->GetLevels(0, 0.01, 50, levels, &filled, 0);
	free(levels);
	delete cubic;
	//sqlite3_exec(db, "SELECT * FROM measuredvalue WHERE segmentid = 2", getinfo, 0, &zErrMsg);
	/*TGlucoseLevel *glucoseLevels;
	sqlite3_exec(db, "SELECT * FROM ", fillLevels, glucoseLevels, &zErrMsg);*/
	sqlite3_close(db);


	return 0;
}