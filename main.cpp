//
// An example used to illustrate calling DICe routines from an external code
//

#include <iostream>
#include <string>
#include <chrono>

// DICe includes for classes used below:
#include <DICe.h>
#include <DICe_Schema.h>

#include "opencv2/opencv.hpp"

// See the custom_app example in the tutorial for more information about the code below

using namespace DICe::field_enums;
using namespace DICe;
using namespace std;
using namespace cv;

float Brightness;
float FrameWidth;
float FrameHeight;

Mat refImage;
Mat defImage;

class SubsetData {
public:
    int Subset_Idx;
    int X_Coord;
    int Y_Coord;
    int Subset_Size;

    double displacement_x;
    double displacement_y;
    double displacement_z;

    SubsetData(int x, int y, int s) {     // Constructor
        X_Coord = x;
        Y_Coord = y;
        Subset_Size = s;
    }
};

list<SubsetData> subSets;

static DICe::Schema dic_init(int argc, char *argv[]);

static void run_2D_dic(DICe::Schema schema);

bool updateListElement(SubsetData *update);

static void load_2D_data(DICe::Schema *schema);


static DICe::Schema dic_init(int argc, char *argv[]) {
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

    return schema;
}

static void load_2D_data(DICe::Schema *schema) {
    //
    // STEP 2:
    //
    // set the reference and deformed images
    //
//        schema->set_ref_image("ref.tif");
//        schema->set_def_image("def.tif");
    /* We are going to attempt to load the files using data in memory instead of writing the files to the drive */

    double refRCP[(int) (FrameHeight * FrameWidth)];
    double defRCP[(int) (FrameHeight * FrameWidth)];

    for (int img_idx = 0; img_idx < FrameHeight * FrameWidth; img_idx++) {
        refRCP[img_idx] = refImage.data[img_idx * 3] * 0.11 + refImage.data[img_idx * 3 + 1] * 0.59 +
                          refImage.data[img_idx * 3 + 2] * 0.30;
        defRCP[img_idx] = defImage.data[img_idx * 3] * 0.11 + defImage.data[img_idx * 3 + 1] * 0.59 +
                          defImage.data[img_idx * 3 + 2] * 0.30;
    }
    {
        int_t l_width = FrameWidth;
        int_t l_height = FrameHeight;

        Teuchos::ArrayRCP<double> fRefRCP(refRCP, 0, FrameHeight * FrameWidth, false);
        Teuchos::ArrayRCP<double> fDefRCP(defRCP, 0, FrameHeight * FrameWidth, false);

        schema->set_ref_image(l_width, l_height, fRefRCP);
        schema->set_def_image(l_width, l_height, fDefRCP, 0);
    }

    //
    // There are also set methods that take DICe::Images as input arguments or pointers to arrays of intensity values.
    // See DICe::Schema.h for these methods
    schema->print_fields();
    std::cout << "HK:" << schema->local_num_subsets() << std::endl;
    subSets.clear();
    for (int subset_idx = 0; subset_idx < schema->local_num_subsets(); subset_idx++) {
        stringstream sx;
        stringstream sy;
        stringstream ss;
        sx << schema->mesh()->get_field(schema->mesh()->get_field_spec("COORDINATE_X"))->local_value(subset_idx);
        sy << schema->mesh()->get_field(schema->mesh()->get_field_spec("COORDINATE_Y"))->local_value(subset_idx);
        SubsetData newOne(stoi(sx.str()),
                          stoi(sy.str()),
                          schema->subset_dim());
        newOne.Subset_Idx = subset_idx;
        subSets.push_front(newOne);
        std::cout << "HK: X:      " << newOne.X_Coord << std::endl;
        std::cout << "HK: Y:      " << newOne.Y_Coord << std::endl;
        std::cout << "HK: Subset: " << newOne.Subset_Size << std::endl;
    }
}

static void run_2D_dic(DICe::Schema schema) {
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
    for (int subset_idx = 0; subset_idx < schema.local_num_subsets(); subset_idx++) {
        SubsetData update_Data(0, 0, 0);
        update_Data.Subset_Idx = subset_idx;
        update_Data.displacement_x = schema.local_field_value(subset_idx, SUBSET_DISPLACEMENT_X_FS);
        update_Data.displacement_y = schema.local_field_value(subset_idx, SUBSET_DISPLACEMENT_Y_FS);
        update_Data.displacement_z = 0;

        std::cout << "The DISPLACEMENT_X field value for subset " << subset_idx << " is "
                  << update_Data.displacement_x << std::endl;
        std::cout << "The DISPLACEMENT_Y field value for subset " << subset_idx << " is "
                  << update_Data.displacement_y << std::endl;
        updateListElement(&update_Data);
    }
    // The field_value() method can be used to set the value as well,
    // for example if you wanted to move subset 0 to a new x-location, the syntax would be
    // schema.local_field_value(0, SUBSET_COORDINATES_X_FS) = 150;

    DICe::finalize();
}

bool updateListElement(SubsetData *update) {
    std::list<SubsetData>::iterator iObject;
    for (iObject = subSets.begin(); iObject != subSets.end(); ++iObject) {
        SubsetData temp = *iObject;
        if (temp.Subset_Idx == update->Subset_Idx) {
            update->Y_Coord = iObject->Y_Coord;
            update->X_Coord = iObject->X_Coord;
            update->Subset_Size = iObject->Subset_Size;
            *iObject = *update;
            return true;
        }
    }
    return false;
}

int main(int argc, char *argv[]) {

    int image_cap_sm = 0;
    int thickness = 1;
    double Total_X_Displacement[3] = {0.0, 0.0, 0.0};
    double Total_Y_Displacement[3] = {0.0, 0.0, 0.0};
    std::chrono::time_point<std::chrono::high_resolution_clock> timer_start, timer_end, main_timer;
    DICe::Schema schema;

    std::cout << "Begin masters DICe program\n";

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
    main_timer = std::chrono::high_resolution_clock::now();
    for (;;) {
        int ti = 0;
        Mat frame;
        Mat show_frame;
        cap >> frame; // get a new frame from camera
        show_frame = frame.clone();
        for (SubsetData &theSet : subSets) {
            Rect r = Rect(theSet.X_Coord - (theSet.Subset_Size / 2), theSet.Y_Coord - (theSet.Subset_Size / 2),
                          theSet.Subset_Size, theSet.Subset_Size);
            rectangle(show_frame, r, Scalar(255, 0, 0), 1, 8, 0);

            arrowedLine(show_frame, Point(theSet.X_Coord, theSet.Y_Coord),
                        Point(theSet.X_Coord + (theSet.displacement_x * 10),
                              theSet.Y_Coord + (theSet.displacement_y * 10)), Scalar(0, 255, 0), 1, 8, 0, 0.1);
            Total_X_Displacement[ti] += theSet.displacement_x;
            Total_Y_Displacement[ti] += theSet.displacement_y;
            thickness = (int) (sqrt(pow(Total_X_Displacement[ti], 2) + pow(Total_Y_Displacement[ti], 2)) / 10);
            if (thickness < 1) {
                thickness = 1;
            }
            arrowedLine(show_frame, Point(theSet.X_Coord, theSet.Y_Coord),
                        Point(theSet.X_Coord + (Total_X_Displacement[ti]),
                              theSet.Y_Coord + (Total_Y_Displacement[ti])), Scalar(0, 0, 255), thickness, 8, 0, 0.1);
            ti++;
        }

        imshow("edges", show_frame);
        char c = waitKey(5);

        if ((c == 'c') || (std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - main_timer).count() > 200)) {
            switch (image_cap_sm) {
                case 0:
                    timer_start = std::chrono::high_resolution_clock::now();
                    image_cap_sm = 1;
//                    imwrite("ref.tif", frame);
                    refImage = frame.clone();
                    std::cout << "Obtain first image\n";
                    break;
                case 1:
                    image_cap_sm = 0;
//                    imwrite("def.tif", frame);
                    defImage = frame.clone();
                    std::cout << "Obtain second image\n";
                    schema = dic_init(argc, argv);
                    load_2D_data(&schema);
                    run_2D_dic(schema);
                    timer_end = std::chrono::high_resolution_clock::now();
                    std::cout << "Completed in: "
                              << std::chrono::duration_cast<std::chrono::milliseconds>(timer_end - timer_start).count()
                              << "ms\n";
                    break;
                default:
                    image_cap_sm = 0;
                    break;
            }
            main_timer = std::chrono::high_resolution_clock::now();
        }

        if (c == 'q') {
            break;
        }
        if (c == 'r') {
            for (int reset_idx = 0; reset_idx < 3; reset_idx++) {
                Total_X_Displacement[reset_idx] = 0;
                Total_Y_Displacement[reset_idx] = 0;
            }
        }
    }
    // the camera will be deinitialized automatically in VideoCapture destructor
    cout << "Program End\n";
    return 0;
}

