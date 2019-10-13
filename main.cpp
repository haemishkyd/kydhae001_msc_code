//
// An example used to illustrate calling DICe routines from an external code
//

#include <iostream>
#include <string>

// DICe includes for classes used below:
#include <DICe.h>
#include <DICe_Schema.h>

#include "opencv2/opencv.hpp"

// See the custom_app example in the tutorial for more information about the code below

using namespace DICe::field_enums;
using namespace std;
using namespace cv;

static DICe::Schema dic_init(int argc, char *argv[]);
static void run_2D_dic(DICe::Schema schema);

float Brightness;
float FrameWidth;
float FrameHeight;
class SubsetData{
public:
    int X_Coord;
    int Y_Coord;
    SubsetData(int x,int y) {     // Constructor
        X_Coord = x;
        Y_Coord = y;
    }
};

list <SubsetData> subSets;

static DICe::Schema dic_init(int argc, char *argv[])
{
    //
    // STEP 0: initialize threading, etc if it is enabled
    //
    DICe::initialize(argc, argv);
    //
    // STEP 1:
    //
    // Create a DICe::Schema that holds all the correlation parameters, image names, solution results, and
    // provides a number of helpful methods (See class DIC::Schema in the documentation for the full list of methods)
    //
    // This step is only needed once, at the beginning of the analysis (not for multiple images or frames in a video).
    // If later, the parameters should change for the analysis, use the set_params(file_name)
    // or set_params(parameterlist) methods to change the parameters
    //
    // Using the parameters file method to create a schema (See params.xml for the correlation parameters and input.xml for the subset locations, etc.)
    DICe::Schema schema("input.xml", "params.xml");
    //
    // There are two different ways to set the parameters in the constructor, using an xml file (used here, see above)
    // or by creating a Teuchos::ParameterList manually and setting the parameters.
    // To use a Teuchos::ParameterList, #include<Teuchos_ParameterList.hpp> and link to library teuchosparameterlist (from Trilinos)
    // If the second method is used, set parameters as follows:
    //
    //     Teuchos::RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList());
    //     params->set("enable_rotation",false);
    //     params->set("enable_normal_strain",true);
    //     params->set("interpolation_method", DICe::BILINEAR);
    //     ... and so on for the rest of the desired parameters
    //
    // The schema constructor would then be
    //
    // DICe::Schema("ref.tif","def.tif",params);
    //
    // STEP 2:
    //
    // set the reference and deformed images
    //
    schema.set_ref_image("ref.tif");
    schema.set_def_image("def.tif");
    //
    // There are also set methods that take DICe::Images as input arguments or pointers to arrays of intensity values.
    // See DICe::Schema.h for these methods
    schema.print_fields();
    std::cout << "HK:" << schema.local_num_subsets() << std::endl;
    for(int subset_idx = 0;subset_idx<schema.local_num_subsets();subset_idx++){
        stringstream sx;
        stringstream sy;
        sx << schema.mesh()->get_field(schema.mesh()->get_field_spec("COORDINATE_X"))->local_value(subset_idx);
        sy << schema.mesh()->get_field(schema.mesh()->get_field_spec("COORDINATE_Y"))->local_value(subset_idx);
        SubsetData newOne(stoi(sx.str()),
                          stoi(sy.str()));
        subSets.push_front(newOne);
        std::cout << "HK: X:" << newOne.X_Coord << std::endl;
        std::cout << "HK: Y:" << newOne.Y_Coord << std::endl;
    }

    return schema;
}

static void run_2D_dic(DICe::Schema schema)
{
    //
    // STEP 3:
    //
    // Run the analysis
    //
    schema.execute_correlation();
    // post process the strain, etc.
    schema.execute_post_processors();
    //
    // STEP 4:
    //
    // Write the output
    schema.write_output("", "custom_app_output");

    std::cout << "End  masters DICe program\n";

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // The code below this point either provides extra illustrations of schema methods or is used
    // to test that the custom_app example is building and executing properly in the DICe regression test system

    //
    // Examples of other operations
    //
    // Change the deformed image
    //schema.set_def_image("def2.tif");
    // after this, execute_correlation can be called again on this image without having to re-init the schema
    //
    // Direct access to field values in the schema
    // schema.field_value( global_subset_id, field_name)
    std::cout << "The DISPLACEMENT_X field value for subset 0 is "
              << schema.local_field_value(0, SUBSET_DISPLACEMENT_X_FS) << std::endl;
    // The field_value() method can be used to set the value as well,
    // for example if you wanted to move subset 0 to a new x-location, the syntax would be
    schema.local_field_value(0, SUBSET_COORDINATES_X_FS) = 150;

    DICe::finalize();
}

int main(int argc, char *argv[]) {

    std::cout << "Begin masters DICe program\n";

   int image_cap_sm = 0;

    DICe::Schema schema = dic_init(argc, argv);

    VideoCapture cap(0); // open the default camera

    Brightness = cap.get(CV_CAP_PROP_BRIGHTNESS);
    FrameWidth = cap.get(CV_CAP_PROP_FRAME_WIDTH);
    FrameHeight = cap.get(CV_CAP_PROP_FRAME_HEIGHT);

    cout << "====================================" << endl << endl;
    cout << "Default Brightness -------> " << Brightness << endl;
    cout << "Default Width      -------> " << FrameWidth << endl;
    cout << "Default Height     -------> " << FrameHeight << endl;
    cout << "====================================" << endl;


    if (!cap.isOpened()) {  // check if we succeeded
        std::cout << "No camera can be found\n";
        return -1;
    } else {
        cout << "Camera is open\n";
    }

    Mat edges;
    namedWindow("edges", 1);
    for (;;) {
        Mat frame;
        cap >> frame; // get a new frame from camera
        for (const SubsetData & theSet : subSets)
        {
            //FIXME: We need to get the size from the file
            Rect r=Rect(theSet.X_Coord-14,theSet.Y_Coord-14,28,28);
            rectangle(frame,r,Scalar(255,0,0),1,8,0);
        }

        imshow("edges", frame);
        char c = waitKey(5);

        if (c == 'c') {
            switch (image_cap_sm) {
                case 0:
                    image_cap_sm++;
                    imwrite("ref.tif",frame);
                    break;
                case 1:
                    image_cap_sm = 0;
                    imwrite("def.tif",frame);
                    run_2D_dic(schema);
                    break;
                default:
                    image_cap_sm = 0;
                    break;
            }
        }

        if (c == 'q') {
            break;
        }
    }
    // the camera will be deinitialized automatically in VideoCapture destructor
    cout << "Program End\n";
    return 0;
}

