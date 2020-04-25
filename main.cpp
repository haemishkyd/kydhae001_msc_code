// @HEADER
// ************************************************************************
//
//               Digital Image Correlation Engine (DICe)
//                 Copyright 2015 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact: Dan Turner (dzturne@sandia.gov)
//
// ************************************************************************
// @HEADER

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

#if DICE_MPI
#  include <mpi.h>
#endif

using namespace DICe::field_enums;
using namespace DICe;
using namespace cv;
using namespace std;

bool read_input_data_files();
void run_cross_correlation();
bool run_correlation_and_triangulation(int image_it);
void main_stereo_3d_correlation();
void information_extraction();
void create_image_in_memory(double *refRCP, double *defRCP);

typedef struct{
    int num_frames;
    std::string file_prefix;
    std::string stereo_file_prefix;
    bool is_stereo;
    int proc_size;
    int proc_rank;
    int FrameWidth;
    int FrameHeight;
}MainDataStructType;

Teuchos::RCP<DICe::Schema> schema;
Teuchos::RCP<DICe::Schema> stereo_schema;
Teuchos::RCP<DICe::Triangulation> triangulation;
Teuchos::RCP<Teuchos::ParameterList> input_params;
Teuchos::RCP<Teuchos::ParameterList> correlation_params;
std::vector<std::string> image_files;
std::vector<std::string> stereo_image_files;
Teuchos::RCP<std::ostream> outStream;
MainDataStructType MainDataStruct;
/**
 * Start a whole bunch of metrics to see how long stuff takes!!
 */
Teuchos::RCP<Teuchos::Time> total_time = Teuchos::TimeMonitor::getNewCounter("## Total Time ##");
Teuchos::RCP<Teuchos::Time> cross_time = Teuchos::TimeMonitor::getNewCounter("Cross-correlation");
Teuchos::RCP<Teuchos::Time> corr_time = Teuchos::TimeMonitor::getNewCounter("Correlation");
Teuchos::RCP<Teuchos::Time> write_time = Teuchos::TimeMonitor::getNewCounter("Write Output");
Mat frame1, frame2, data(500, 1200, CV_8UC3, Scalar(0, 0, 0));;

bool read_input_data_files() {
    /**
     * Get all of the input parameters from the input files.
     */
    input_params = Teuchos::rcp(new Teuchos::ParameterList());
    Teuchos::Ptr<Teuchos::ParameterList> inputParamsPtr(input_params.get());
    Teuchos::updateParametersFromXmlFile("input.xml", inputParamsPtr);
    TEUCHOS_TEST_FOR_EXCEPTION(input_params == Teuchos::null, std::runtime_error, "");

    *outStream << "Input Parameters: " << std::endl;
    input_params->print(*outStream);
    *outStream << "\n--- Input read successfully ---\n" << std::endl;

    /**
     * Get all of the correlation parameters from the input files.
     */
    bool is_error_est_run = false;
    if (input_params->isParameter(DICe::correlation_parameters_file)) {
        const std::string paramsFileName = input_params->get<std::string>(DICe::correlation_parameters_file);
        correlation_params = DICe::read_correlation_params(paramsFileName);
        *outStream << "User specified correlation Parameters: " << std::endl;
        correlation_params->print(*outStream);
        is_error_est_run = correlation_params->get<bool>(DICe::estimate_resolution_error, false);
        if (is_error_est_run) {
            // force the computing of the image laplacian for the reference image:
            correlation_params->set(DICe::compute_laplacian_image, true);
        }
        *outStream << "\n--- Correlation parameters read successfully ---\n" << std::endl;
    } else {
        *outStream << "Correlation parameters not specified by user" << std::endl;
    }

    return is_error_est_run;
}

void run_cross_correlation() {
    /* We know this is a stereo analysis so we just assume all is correct */
    Teuchos::TimeMonitor cross_time_monitor(*cross_time);
    TEUCHOS_TEST_FOR_EXCEPTION(schema->analysis_type() == GLOBAL_DIC, std::runtime_error,
                               "Error, global stereo not enabled yet");
    *outStream << "Processing cross correlation between left and right images" << std::endl;
    schema->initialize_cross_correlation(triangulation,
                                              input_params); // images don't need to be loaded by here they get loaded in this routine based on the input params
    schema->update_extents(true);

    /* Create the image in memory */

    double refRCP[(int) (MainDataStruct.FrameHeight * MainDataStruct.FrameWidth)];
    double defRCP[(int) (MainDataStruct.FrameHeight * MainDataStruct.FrameWidth)];

    create_image_in_memory(refRCP,defRCP);
    Teuchos::ArrayRCP<double> fRefRCP(refRCP, 0, MainDataStruct.FrameHeight * MainDataStruct.FrameWidth, false);
    Teuchos::ArrayRCP<double> fDefRCP(defRCP, 0, MainDataStruct.FrameHeight * MainDataStruct.FrameWidth, false);

    schema->set_ref_image(MainDataStruct.FrameWidth, MainDataStruct.FrameHeight, fRefRCP);
    schema->set_def_image(MainDataStruct.FrameWidth, MainDataStruct.FrameHeight, fDefRCP, 0);

    schema->set_ref_image(image_files[0]);
    schema->set_def_image(stereo_image_files[0]);
    if (schema->use_nonlinear_projection()) {
        schema->project_right_image_into_left_frame(triangulation, false);
    }
    schema->execute_cross_correlation();
    schema->save_cross_correlation_fields();
    stereo_schema = Teuchos::rcp(new DICe::Schema(input_params, correlation_params, schema));
    stereo_schema->update_extents();
    stereo_schema->set_ref_image(MainDataStruct.FrameWidth, MainDataStruct.FrameHeight, fDefRCP);
//    stereo_schema->set_ref_image(stereo_image_files[0]);
    assert(stereo_schema != Teuchos::null);
    //if(stereo_schema->use_nonlinear_projection())
    //  stereo_schema->project_right_image_into_left_frame(triangulation,true);
    stereo_schema->set_frame_range(0, 2);

    // go ahead and set up the model coordinates field
    schema->execute_triangulation(triangulation, schema);
}

void create_image_in_memory(double *refRCP, double *defRCP){
    for (int img_idx = 0; img_idx < MainDataStruct.FrameHeight * MainDataStruct.FrameWidth; img_idx++) {
        refRCP[img_idx] = frame1.data[img_idx * 3] * 0.11 + frame1.data[img_idx * 3 + 1] * 0.59 +
                frame1.data[img_idx * 3 + 2] * 0.30;
        defRCP[img_idx] = frame2.data[img_idx * 3] * 0.11 + frame2.data[img_idx * 3 + 1] * 0.59 +
                frame2.data[img_idx * 3 + 2] * 0.30;
    }
}

bool run_correlation_and_triangulation(int image_it) {
    bool failed_step = false;
    bool is_stereo = true; //This is definitely stereo

    std::string file_prefix = input_params->get<std::string>(DICe::output_prefix, "DICe_solution");
    std::string stereo_file_prefix = input_params->get<std::string>(DICe::output_prefix, "DICe_solution");
    stereo_file_prefix += "_stereo";
    const bool separate_header_file = input_params->get<bool>(DICe::create_separate_run_info_file, false);

    *outStream << "Processing frame: " << image_it << ", " << image_files[image_it]
               << std::endl;
    if (schema->use_incremental_formulation() && image_it > 1) {
        schema->set_ref_image(schema->def_img());
    }
    schema->update_extents();
    schema->set_def_image(image_files[image_it]);
    if (is_stereo) {
        if (stereo_schema->use_incremental_formulation() && image_it > 1) {
            stereo_schema->set_ref_image(stereo_schema->def_img());
        }
        stereo_schema->update_extents();
        stereo_schema->set_def_image(stereo_image_files[image_it]);
        //if(stereo_schema->use_nonlinear_projection())
        //  stereo_schema->project_right_image_into_left_frame(triangulation,false);
    }
    { // start the timer
        Teuchos::TimeMonitor corr_time_monitor(*corr_time);
        int_t corr_error = schema->execute_correlation();
        if (corr_error)
            failed_step = true;
        if (is_stereo) {
            corr_error = stereo_schema->execute_correlation();
            if (corr_error)
                failed_step = true;
        }
        schema->execute_triangulation(triangulation, stereo_schema);
        schema->execute_post_processors();
    }
    // write the output
    const bool no_text_output = input_params->get<bool>(DICe::no_text_output_files, false);
    {
        Teuchos::TimeMonitor write_time_monitor(*write_time);
        schema->write_output(output_folder, file_prefix, separate_output_file_for_each_subset,
                                  separate_header_file, no_text_output);
        schema->post_execution_tasks();
        // print the timing data with or without verbose flag
        if (input_params->get<bool>(DICe::print_stats, false)) {
            schema->mesh()->print_field_stats();
        }
        //if(subset_info->conformal_area_defs!=Teuchos::null&&image_it==1){
        //  schema->write_control_points_image("RegionOfInterest");
        //}
        if (is_stereo) {
            if (input_params->get<bool>(DICe::output_stereo_files, false)) {
                stereo_schema->write_output(output_folder, stereo_file_prefix,
                                                 separate_output_file_for_each_subset, separate_header_file,
                                                 no_text_output);
            }
            stereo_schema->post_execution_tasks();
        }
    }
    return failed_step;
}


void information_extraction(){
    std::string output_folder;
    bool is_error_est_run;
    int_t proc_size = 1;
    int_t proc_rank = 0;

    std::cout << "Start of process." << std::endl;

    outStream = Teuchos::rcp(&std::cout, false);

    if (proc_rank == 0) DEBUG_MSG("Parsing command line options");
    bool force_exit = false;

    /******* Get the input parameters */
    is_error_est_run = read_input_data_files();

    /******* Decipher the image file names (note: zero entry is the reference image) */
    DICe::decipher_image_file_names(input_params, image_files, stereo_image_files);
    const bool is_stereo = stereo_image_files.size() > 0;

    /******* Create the list of images */
    const int_t num_frames = image_files.size() - 1;
    int_t first_frame_id = 0;
    int_t image_width = 0;
    int_t image_height = 0;

    TEUCHOS_TEST_FOR_EXCEPTION(num_frames <= 0, std::runtime_error, "");
    *outStream << "Reference image: " << image_files[0] << std::endl;
    for (int_t i = 1; i <= num_frames; ++i) {
        if (i == 10 && num_frames != 10) *outStream << "..." << std::endl;
        else if (i > 10 && i < num_frames) continue;
        else
            *outStream << "Deformed image: " << image_files[i] << std::endl;
    }
    *outStream << "\n--- List of images constructed successfuly ---\n" << std::endl;

    /******* Get the size of the images being used */
    utils::read_image_dimensions(image_files[0].c_str(), image_width, image_height);
    *outStream << "Image dimensions: " << image_width << " x " << image_height << std::endl;

    /******* Where are we going to put the output information */
    output_folder = input_params->get<std::string>(DICe::output_folder);
    const bool separate_output_file_for_each_subset = input_params->get<bool>(
            DICe::separate_output_file_for_each_subset, false);
    if (separate_output_file_for_each_subset) {
        *outStream << "Output will be written to separate output files for each subset" << std::endl;
    } else {
        *outStream << "Output will be written to one file per frame with all subsets included" << std::endl;
    }
    const bool separate_header_file = input_params->get<bool>(DICe::create_separate_run_info_file, false);
    if (separate_header_file) {
        *outStream
                << "Execution information will be written to a separate file (not placed in the output headers)"
                << std::endl;
    }

    /******* create schemas: */
    schema = Teuchos::rcp(new DICe::Schema(input_params, correlation_params));
    // let the schema know how many images there are in the sequence and the first frame id:
    schema->set_frame_range(first_frame_id, num_frames);

    /******* Set up the subsets */
    *outStream << "Number of global subsets: " << schema->global_num_subsets() << std::endl;
    list<SubSetData> *subSets = getSubSets();
    subSets->clear();
    for (int_t i = 0; i < schema->local_num_subsets(); ++i) {

        stringstream sx;
        stringstream sy;
        stringstream ss;
        sx << schema->mesh()->get_field(schema->mesh()->get_field_spec("COORDINATE_X"))->local_value(i);
        sy << schema->mesh()->get_field(schema->mesh()->get_field_spec("COORDINATE_Y"))->local_value(i);
        SubSetData newOne(stoi(sx.str()),
                          stoi(sy.str()),
                          schema->subset_dim());
        newOne.Subset_Idx = i;
        subSets->push_front(newOne);

        if (i == 10 && schema->local_num_subsets() != 11) *outStream << "..." << std::endl;
        else if (i > 10 && i < schema->local_num_subsets() - 1) continue;
        else
            *outStream << "Proc 0: subset global id: " << schema->subset_global_id(i) << " global coordinates ("
                       << schema->local_field_value(i, DICe::field_enums::SUBSET_COORDINATES_X_FS) <<
                       "," << schema->local_field_value(i, DICe::field_enums::SUBSET_COORDINATES_Y_FS) << ")"
                       << std::endl;
    }
    *outStream << std::endl;

    std::string file_prefix = input_params->get<std::string>(DICe::output_prefix, "DICe_solution");
    std::string stereo_file_prefix = input_params->get<std::string>(DICe::output_prefix, "DICe_solution");
    stereo_file_prefix += "_stereo";

    // for backwards compatibility allow the user to specify either a calibration_parameters_file or a camera_system_file
    // (camera_system_file is the new preferred way)
    TEUCHOS_TEST_FOR_EXCEPTION(is_stereo && (!input_params->isParameter(DICe::calibration_parameters_file) &&
                                             !input_params->isParameter(DICe::camera_system_file)),
                               std::runtime_error,
                               "Error, calibration_parameters_file or camera_system_file required for stereo");
    TEUCHOS_TEST_FOR_EXCEPTION(input_params->isParameter(DICe::calibration_parameters_file) &&
                               input_params->isParameter(DICe::camera_system_file),
                               std::runtime_error,
                               "Error, both calibration_parameters_file and camera_system_file cannot be specified");

    if (input_params->isParameter(DICe::calibration_parameters_file) ||
        input_params->isParameter(DICe::camera_system_file)) {
        if (proc_rank == 0)
            update_legacy_txt_cal_input(
                    input_params); // in case an old txt format cal input file is being used it needs to have width and height added to it

        const std::string cal_file_name = input_params->isParameter(DICe::calibration_parameters_file)
                                          ? input_params->get<std::string>(DICe::calibration_parameters_file) :
                                          input_params->get<std::string>(DICe::camera_system_file);
        triangulation = Teuchos::rcp(new DICe::Triangulation(cal_file_name));
        *outStream << "\n--- Calibration parameters read successfully ---\n" << std::endl;
    } else {
        *outStream << "Calibration parameters not specified by user" << std::endl;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(is_stereo && triangulation == Teuchos::null, std::runtime_error,
                               "Error, triangulation should be instantiated at this point");
    MainDataStruct.num_frames = num_frames;
    MainDataStruct.is_stereo = is_stereo;
    MainDataStruct.file_prefix = file_prefix;
    MainDataStruct.stereo_file_prefix = stereo_file_prefix;
    MainDataStruct.proc_rank = proc_rank;
    MainDataStruct.proc_size = proc_size;
}


void main_stereo_3d_correlation() {

    // iterate through the images and perform the correlation:
    bool failed_step;

    for (int_t image_it = 1; image_it <= MainDataStruct.num_frames; ++image_it) {
        failed_step = run_correlation_and_triangulation(image_it);
    } // image loop

    schema->write_stats(output_folder, MainDataStruct.file_prefix);
    if (MainDataStruct.is_stereo)
        stereo_schema->write_stats(output_folder, MainDataStruct.stereo_file_prefix);

    if (failed_step)
        *outStream << "\n--- Failed Step Occurred ---\n" << std::endl;
    else
        *outStream << "\n--- Successful Completion ---\n" << std::endl;

    // output timing

    // print the timing data with or without verbose flag
    if (input_params->get<bool>(DICe::print_timing, false)) {
        Teuchos::TimeMonitor::summarize(*outStream, false, true, false/*zero timers*/);
    }
    //  write the time output to file:
    std::stringstream timeFileName;
    timeFileName << output_folder << "timing." << MainDataStruct.proc_size << "." << MainDataStruct.proc_rank << ".txt";
    std::ofstream ofs(timeFileName.str(), std::ofstream::out);
    Teuchos::TimeMonitor::summarize(ofs, false, true, false/*zero timers*/);
    ofs.close();
    if (MainDataStruct.proc_rank != 0) // only keep the process zero copy of the timing results
        std::remove(timeFileName.str().c_str());

    DICe::finalize();
}

int main(int argc, char *argv[]) {
    int return_val;
    float Brightness;
    int system_state = 0;

    VideoCapture cap2(0); // open the default camera
    VideoCapture cap1(2); // open the default camera

    Brightness = cap1.get(CV_CAP_PROP_BRIGHTNESS);
    MainDataStruct.FrameWidth = cap1.get(CV_CAP_PROP_FRAME_WIDTH);
    MainDataStruct.FrameHeight = cap1.get(CV_CAP_PROP_FRAME_HEIGHT);

    cout << "====================================" << endl << endl;
    cout << "Default Brightness -------> " << Brightness << endl;
    cout << "Default Width      -------> " << MainDataStruct.FrameWidth << endl;
    cout << "Default Height     -------> " << MainDataStruct.FrameHeight << endl;
    cout << "====================================" << endl;

    if (!cap1.isOpened()) {  // check if we succeeded
        std::cout << "First camera cannot be found\n";
        return -1;
    } else {
        cout << "Camera 1 is open\n";
    }
    if (!cap2.isOpened()) {  // check if we succeeded
        std::cout << "Second camera cannot be found\n";
        return -1;
    } else {
        cout << "Camera 2 is open\n";
    }

    namedWindow("Left", WINDOW_AUTOSIZE);
    namedWindow("Right", WINDOW_AUTOSIZE);
    for (;;) {
        switch (system_state) {
            case 0:
                cap2 >> frame1; // get a new frame from camera
                cap1 >> frame2;

                imshow("Left", frame1);
                imshow("Right", frame2);
                break;
            case 1:
                cap2 >> frame1; // get a new frame from camera
                cap1 >> frame2;
                DICe::initialize(argc, argv);
                information_extraction();
                run_cross_correlation();
                system_state = 2;
                break;
            case 2:
                cap2 >> frame1; // get a new frame from camera
                cap1 >> frame2;

                imwrite("Img_0001_0.jpeg", frame1);
                imwrite("Img_0001_1.jpeg", frame2);

                list<SubSetData> *subSets = getSubSets();
                if (subSets->size() > 0) {
                    for (SubSetData &theSet : (*subSets)) {
                        Rect r = Rect(theSet.X_Coord - (theSet.Subset_Size / 2),
                                      theSet.Y_Coord - (theSet.Subset_Size / 2),
                                      theSet.Subset_Size, theSet.Subset_Size);
                        rectangle(frame1, r, Scalar(255, 0, 0), 1, 8, 0);
                    }
                }

                imshow("Left", frame1);
                imshow("Right", frame2);
                imshow("Data", data);

                main_stereo_3d_correlation();
                data.release();
                data = Mat(500, 1200, CV_8UC3, Scalar(0, 0, 0));
                putText(data, "Subset 1", Point(0, 30), FONT_HERSHEY_SIMPLEX, 1, Scalar(128));
                putText(data, "Subset 2", Point(400, 30), FONT_HERSHEY_SIMPLEX, 1, Scalar(128));
                putText(data, "Subset 3", Point(800, 30), FONT_HERSHEY_SIMPLEX, 1, Scalar(128));
                for (int subset_idx = 0; subset_idx < schema->local_num_subsets(); subset_idx++) {

                    stringstream sx;
                    stringstream sy;
                    stringstream sz;
                    sx << schema->local_field_value(subset_idx, MODEL_DISPLACEMENT_X_FS);
                    sy << schema->local_field_value(subset_idx, MODEL_DISPLACEMENT_Y_FS);
                    sz << schema->local_field_value(subset_idx, MODEL_DISPLACEMENT_Z_FS);
                    putText(data, "X:", Point(subset_idx * 400, 80), FONT_HERSHEY_SIMPLEX, 1, Scalar(128));
                    putText(data, sx.str(), Point(subset_idx * 400 + 100, 80), FONT_HERSHEY_SIMPLEX, 1, Scalar(128));
                    putText(data, "Y:", Point(subset_idx * 400, 130), FONT_HERSHEY_SIMPLEX, 1, Scalar(128));
                    putText(data, sy.str(), Point(subset_idx * 400 + 100, 130), FONT_HERSHEY_SIMPLEX, 1, Scalar(128));
                    putText(data, "Z:", Point(subset_idx * 400, 180), FONT_HERSHEY_SIMPLEX, 1, Scalar(128));
                    putText(data, sz.str(), Point(subset_idx * 400 + 100, 180), FONT_HERSHEY_SIMPLEX, 1, Scalar(128));
                }
                break;
        }
        char c = waitKey(5);

        if (c == 'q') {
            break;
        }
        if (c == 'i') {
            imwrite("Img_0000_0.jpeg", frame1);
            imwrite("Img_0000_1.jpeg", frame2);
            system_state = 1;
        }
    }


    return return_val;
}
