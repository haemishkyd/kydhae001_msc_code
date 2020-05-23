//
// Created by haemish on 2020/04/11.
//

using namespace std;
using namespace cv;

class SubSetData {
public:
    int Subset_Idx;
    int X_Coord;
    int Y_Coord;
    int Subset_Size;

    double displacement_x;
    double displacement_y;
    double displacement_z;

    SubSetData(int x, int y, int s) {     // Constructor
        X_Coord = x;
        Y_Coord = y;
        Subset_Size = s;
    }


};

list<SubSetData> *getSubSets();
void drawSubsets(Mat *passedFrame, Teuchos::RCP<DICe::Schema> *passedSchema);

