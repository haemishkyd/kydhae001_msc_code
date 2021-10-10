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

vector<int> rgb(double ratio)
{
    vector<int> rgb_array;
    //we want to normalize ratio so that it fits in to 6 regions
    //where each region is 256 units long
    int normalized = int(ratio * 256 * 5);

    //find the distance to the start of the closest region
    int x = normalized % 256;

    int red = 0, grn = 0, blu = 0;
    // switch (normalized / 256)
    // {
    // case 0:
    //     red = 255;
    //     grn = x;
    //     blu = 0;
    //     break; //red
    // case 1:
    //     red = 255 - x;
    //     grn = 255;
    //     blu = 0;
    //     break; //yellow
    // case 2:
    //     red = 0;
    //     grn = 255;
    //     blu = x;
    //     break; //green
    // case 3:
    //     red = 0;
    //     grn = 255 - x;
    //     blu = 255;
    //     break; //cyan
    // case 4:
    //     red = x;
    //     grn = 0;
    //     blu = 255;
    //     break; //blue
    // case 5:
    //     red = 255;
    //     grn = 0;
    //     blu = 255 - x;
    //     break; //magenta
    // }
    switch(normalized/256)
    {
    case 0:
        red = 0;
        grn = x;
        blu = 255;
        break; //blue
    case 1:
        red = 0;
        grn = 255;
        blu = 255 - x;
        break; //cyan
    case 2:
        red = x;
        grn = 255;
        blu = 0;
        break; //green
    case 3:
        red = 255;
        grn = 255 - x;
        blu = 0;
        break; //yellow
    case 4:
        red = 255;
        grn = 0;
        blu = x;
        break; //red
    }

    rgb_array.push_back(red);
    rgb_array.push_back(grn);
    rgb_array.push_back(blu);
    return rgb_array;
}

void drawSubsets(Mat *passedFrame, Teuchos::RCP<DICe::Schema> *passedSchema, int which_axis_to_draw){
    double disp_value;
    vector<int> temp_col;
    string axis_description;
    
    if (subSets.size() > 0)
    {
        switch (which_axis_to_draw)
        {
        case 0:
            axis_description = "Z Displacement";            
            break;
        case 1:
            axis_description = "Y Displacement";
            break;
        case 2:
            axis_description = "X Displacement";
            break;
        }
        Scalar textColour = Scalar(255, 255, 0);
        putText(*passedFrame, axis_description, Point(10, 80), FONT_HERSHEY_SIMPLEX, 2, textColour);
        for (SubSetData &theSet : (subSets)) {
            Rect r = Rect(theSet.X_Coord - (theSet.Subset_Size / 2),
                          theSet.Y_Coord - (theSet.Subset_Size / 2),
                          theSet.Subset_Size, theSet.Subset_Size);
            

            switch (which_axis_to_draw)
            {
            case 0:                
                disp_value = (double)((*passedSchema)->local_field_value(theSet.Subset_Idx, MODEL_DISPLACEMENT_Z_FS));
                break;
            case 1:
                disp_value = (double)((*passedSchema)->local_field_value(theSet.Subset_Idx, MODEL_DISPLACEMENT_Y_FS));
                break;
            case 2:
                disp_value = (double)((*passedSchema)->local_field_value(theSet.Subset_Idx, MODEL_DISPLACEMENT_X_FS));
                break;
            default:
                disp_value = (double)((*passedSchema)->local_field_value(theSet.Subset_Idx, MODEL_DISPLACEMENT_Z_FS));
                break;
            }
            temp_col = rgb((double)(disp_value / 15.0));
            if (theSet.Subset_Idx = 0){
                cout << "Subset " << theSet.Subset_Idx << ":" << disp_value << endl;
                cout << "B:" << temp_col[2] << " G:" << temp_col[1] << " R:" << temp_col[0] << endl;
            }
                
            rectangle(*passedFrame, r, Scalar(temp_col[2], temp_col[1], temp_col[0]), -1, LINE_8);
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