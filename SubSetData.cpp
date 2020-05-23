//
// Created by haemish on 2020/04/11.
//
#include <DICe_ImageIO.h>
#include <DICe_Schema.h>
#include <DICe_Triangulation.h>


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
            rectangle(*passedFrame, r, Scalar(255, 0, 0), 1, 8, 0);
        }
    }
}
