//
// Created by haemish on 2020/04/11.
//
#include <DICe.h>
#include <DICe_Parser.h>
#include <DICe_Image.h>
#include <DICe_ImageIO.h>
#include <DICe_Schema.h>
#include <DICe_Triangulation.h>


#include <fstream>

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include "opencv2/opencv.hpp"

#include "SubSetData.h"

using namespace std;

static list<SubSetData> subSets;

list<SubSetData> *getSubSets(){
    return &subSets;
}
