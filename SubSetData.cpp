//
// Created by haemish on 2020/04/11.
//
#include <DICe_ImageIO.h>
#include <DICe_Schema.h>
#include <DICe_Triangulation.h>

#include <fstream>
#include "opencv2/opencv.hpp"

#include "SubSetData.h"

using namespace DICe::field_enums;
using namespace DICe;
using namespace std;
using namespace cv;

static list<SubSetData> subSets;

list<SubSetData> *getSubSets(){
    return &subSets;
}

void drawSubsets(Mat *passedFrame, Teuchos::RCP<DICe::Schema> *passedSchema){
    if (subSets.size() > 0) {
        for (SubSetData &theSet : (subSets)) {
            Rect r = Rect(theSet.X_Coord - (theSet.Subset_Size / 2),
                          theSet.Y_Coord - (theSet.Subset_Size / 2),
                          theSet.Subset_Size, theSet.Subset_Size);
            std::vector<Mat> channels(3);
            cv:split(*passedFrame,channels);
            Mat extractedRoi;
            extractedRoi = channels.at(2)(r);
            extractedRoi += char(((*passedSchema)->local_field_value(theSet.Subset_Idx, MODEL_DISPLACEMENT_Z_FS))*-5);
            cv::merge(channels, *passedFrame);
            rectangle(*passedFrame, r, Scalar(255, 255, 0), 1, 8, 0);
        }
    }
}

void SubSetData::SetStartCorner(int x, int y){
    X_Start_Coord = x;
    Y_Start_Coord = y;
}

void SubSetData::SetEndCorner(int x, int y){
    X_End_Coord = x;
    Y_End_Coord = y;
}

void SubSetData::GenerateSubsetFile(){
    int x_accum,y_accum;
    std::stringstream subsetFileName;
    subsetFileName << "subsets.txt";
    std::ofstream ofs(subsetFileName.str(), std::ofstream::out);
    x_accum = X_Start_Coord;
    y_accum = Y_Start_Coord;
    ofs << "BEGIN SUBSET_COORDINATES" << endl;
    while (y_accum < Y_End_Coord){
        x_accum = X_Start_Coord;
        while (x_accum < X_End_Coord){
            ofs <<"    "<< x_accum << " " << y_accum << endl;
            x_accum += 31;
        }
        y_accum += 31;
    }
    ofs << "END" << endl;
    ofs.close();
}