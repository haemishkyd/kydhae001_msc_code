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

    int X_Start_Coord;
    int Y_Start_Coord;
    int X_End_Coord;
    int Y_End_Coord;

    int X_Start_Draw_Coord;
    int Y_Start_Draw_Coord;
    int X_End_Draw_Coord;
    int Y_End_Draw_Coord;

    double displacement_x;
    double displacement_y;
    double displacement_z;

    SubSetData(int x, int y, int s) {     // Constructor
        X_Coord = x;
        Y_Coord = y;
        Subset_Size = s;
    }

    void SetStartCorner(int x, int y);
    void SetEndCorner(int x, int y);
    void GenerateSubsetFile();
};

list<SubSetData> *getSubSets();
void drawSubsets(Mat *passedFrame, Teuchos::RCP<DICe::Schema> *passedSchema);

