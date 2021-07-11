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

#include <pthread.h>
#include <DICe.h>
#include <DICe_Parser.h>
#include <DICe_Image.h>
#include <DICe_ImageIO.h>
#include <DICe_Schema.h>
#include <DICe_Triangulation.h>


#include <fstream>

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <SerialPort.h>

#include "opencv2/opencv.hpp"

#include "SubSetData.h"
#include "Semaphore.h"
#include "ScriptRun.h"
#include "calibrate.h"

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

void write_timing_metrics();

void read_script(string script_name);

void dataMouseCallBack(int event, int x, int y, int flags, void* userdata);

void projectionMouseCallBack(int event, int x, int y, int flags, void* userdata);

void outputImageInformation();

string gstreamer_pipeline (int capture_width, int capture_height, int display_width, int display_height, int framerate, int flip_method);

typedef struct {
    int num_frames;
    std::string file_prefix;
    std::string stereo_file_prefix;
    bool is_stereo;
    int proc_size;
    int proc_rank;
    int FrameWidth;
    int FrameHeight;
    bool WriteThreadRunning;
    vector<scalar_t> Right_X;
    vector<scalar_t> Right_Y;
    vector<scalar_t> Left_X;
    vector<scalar_t> Left_Y;
    scalar_t User_Right_X;
    scalar_t User_Right_Y;
    scalar_t User_Left_X;
    scalar_t User_Left_Y;
    vector<string> script_stack;
    ScriptRun *myScript;
    bool use_arducam;
} MainDataStructType;

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

Teuchos::RCP<Teuchos::Time> cross_time = Teuchos::TimeMonitor::getNewCounter("Cross-Correlation and Triangulation");
Teuchos::RCP<Teuchos::Time> state_2_timer = Teuchos::TimeMonitor::getNewCounter("Entire Process Cycle Time");
Teuchos::RCP<Teuchos::Time> total_corr_loop = Teuchos::TimeMonitor::getNewCounter("Total Correlation Time");
Teuchos::RCP<Teuchos::Time> corr_time = Teuchos::TimeMonitor::getNewCounter("Correlation");
Teuchos::RCP<Teuchos::Time> write_time = Teuchos::TimeMonitor::getNewCounter("Write Output");
stringstream cross_and_trian_time_str;
stringstream state_2_total_time_str;

Rect z_button, x_button, y_button;
bool x_button_clicked, y_button_clicked, z_button_clicked;

Mat frame1, frame2, data(500, 1200, CV_8UC3, Scalar(0, 0, 0));;

Semaphore WriteComplete(0);
Semaphore ReadComplete(0);

VideoCapture cap2; // open the default camera
VideoCapture cap1; // open the default camera


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

    schema->set_ref_image(image_files[0]);
    schema->set_def_image(stereo_image_files[0]);
    if (schema->use_nonlinear_projection()) {
        schema->project_right_image_into_left_frame(triangulation, false);
    }
    schema->execute_cross_correlation();
    schema->save_cross_correlation_fields();
    stereo_schema = Teuchos::rcp(new DICe::Schema(input_params, correlation_params, schema));
    stereo_schema->update_extents();
    stereo_schema->set_ref_image(stereo_image_files[0]);
    assert(stereo_schema != Teuchos::null);
    //if(stereo_schema->use_nonlinear_projection())
    //  stereo_schema->project_right_image_into_left_frame(triangulation,true);
    stereo_schema->set_frame_range(0, 2);

    // go ahead and set up the model coordinates field
    schema->execute_triangulation(triangulation, schema);
}

bool run_correlation_and_triangulation(int image_it) {
    bool failed_step = false;
    bool is_stereo = true; //This is definitely stereo

    Teuchos::TimeMonitor corr_time_monitor(*corr_time);

    std::string file_prefix = input_params->get<std::string>(DICe::output_prefix, "DICe_solution");
    std::string stereo_file_prefix = input_params->get<std::string>(DICe::output_prefix, "DICe_solution");
    stereo_file_prefix += "_stereo";
    const bool separate_header_file = input_params->get<bool>(DICe::create_separate_run_info_file, false);

    *outStream << "Processing frame: " << image_it << ", " << image_files[image_it]
               << std::endl;

    schema->update_extents();
    schema->set_def_image(image_files[image_it]);
    if (is_stereo) {
        stereo_schema->update_extents();
        stereo_schema->set_def_image(stereo_image_files[image_it]);
    }

    int_t corr_error = schema->execute_correlation();
    if (corr_error)
        failed_step = true;
    if (is_stereo) {
        corr_error = stereo_schema->execute_correlation();
        if (corr_error)
            failed_step = true;
    }
    ReadComplete.notify(0);
    schema->execute_triangulation(triangulation, stereo_schema);
    schema->execute_post_processors();

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


void information_extraction() {
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
    DICe::utils::read_image_dimensions(image_files[0].c_str(), image_width, image_height);
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
    Teuchos::TimeMonitor stereo_3d_corr(*total_corr_loop);


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

    DICe::finalize();
}

void dataMouseCallBack(int event, int x, int y, int flags, void* userdata){
    if (event == EVENT_LBUTTONDOWN)
    {
        if (x_button.contains(Point(x, y)))
        {
            x_button_clicked = true;
            y_button_clicked = false;
            z_button_clicked = false;
        }
        if (y_button.contains(Point(x, y)))
        {
            x_button_clicked = false;
            y_button_clicked = true;
            z_button_clicked = false;
        }
        if (z_button.contains(Point(x, y)))
        {
            x_button_clicked = false;
            y_button_clicked = false;
            z_button_clicked = true;
        }
    }
}

void projectionMouseCallBack(int event, int x, int y, int flags, void* userdata) {
    if (event == EVENT_LBUTTONDOWN) {
        srand(time(0));
        MainDataStruct.Right_X.clear();
        MainDataStruct.Right_Y.clear();
        MainDataStruct.Left_X.clear();
        MainDataStruct.Left_Y.clear();
        for (int num_it = 0; num_it<10; num_it++) {
            scalar_t left_x = (rand() % (640/3))+((640/2)-(640/3/2));
            scalar_t left_y = (rand() % (480/3))+((480/2)-(480/3/2));
            scalar_t right_x = 0;
            scalar_t right_y = 0;
            triangulation->project_left_to_right_sensor_coords(left_x, left_y, right_x, right_y);
            MainDataStruct.Right_X.push_back(right_x);
            MainDataStruct.Right_Y.push_back(right_y);
            MainDataStruct.Left_X.push_back(left_x);
            MainDataStruct.Left_Y.push_back(left_y);
            cout << "Left: " << left_x << ":" << left_y << " Right: " << right_x << ":" << right_y << endl;
        }
    }
    if (event == EVENT_RBUTTONDOWN){
        scalar_t right_x = 0;
        scalar_t right_y = 0;
        triangulation->project_left_to_right_sensor_coords(x, y, right_x, right_y);
        MainDataStruct.User_Left_X = x;
        MainDataStruct.User_Left_Y = y;
        MainDataStruct.User_Right_X = right_x;
        MainDataStruct.User_Right_Y = right_y;
    }
}

void write_timing_metrics() {
    // output timing
    //  write the time output to file:
    std::stringstream timeFileName;
    timeFileName << "timing.txt";
    std::ofstream ofs(timeFileName.str(), std::ofstream::out | std::ofstream::app);

    ofs << cross_and_trian_time_str.str() << "," << state_2_timer.get()->totalElapsedTime() << ","
        << total_corr_loop.get()->totalElapsedTime() << "," << corr_time.get()->totalElapsedTime() << ","
        << write_time.get()->totalElapsedTime() << "," << int(1.0 / state_2_timer.get()->totalElapsedTime())<<endl;

    ofs.close();
}

void read_script(string script_name){
    ifstream myfile;
    std::string line;

    myfile.open (script_name);
    while (std::getline(myfile, line)){
        MainDataStruct.script_stack.push_back(line);
    }
}

void outputImageInformation(){

    drawSubsets(&frame1,&schema);

    Mat OutputFrame;
    Mat DisplayFrame;
    hconcat(frame1,frame2,OutputFrame);
    resize(OutputFrame,DisplayFrame,Size(1280,480));
    data.release();
    data = Mat(500, 1280, CV_8UC3, Scalar(0, 0, 0));
    Scalar textColour = Scalar(255,255,0);
    putText(data, "Subset 1", Point(0, 30), FONT_HERSHEY_SIMPLEX, 1, textColour);
    putText(data, "Subset 2", Point(400, 30), FONT_HERSHEY_SIMPLEX, 1, textColour);
    putText(data, "Subset 3", Point(800, 30), FONT_HERSHEY_SIMPLEX, 1, textColour);

    for (int subset_idx = 0; subset_idx < schema->local_num_subsets(); subset_idx++) {

        stringstream sx;
        stringstream sy;
        stringstream sz;
        sx << schema->local_field_value(subset_idx, MODEL_DISPLACEMENT_X_FS);
        sy << schema->local_field_value(subset_idx, MODEL_DISPLACEMENT_Y_FS);
        sz << schema->local_field_value(subset_idx, MODEL_DISPLACEMENT_Z_FS);
        putText(data, "X:", Point(subset_idx * 400, 80), FONT_HERSHEY_SIMPLEX, 1, textColour);
        putText(data, sx.str(), Point(subset_idx * 400 + 100, 80), FONT_HERSHEY_SIMPLEX, 1, textColour);
        putText(data, "Y:", Point(subset_idx * 400, 130), FONT_HERSHEY_SIMPLEX, 1, textColour);
        putText(data, sy.str(), Point(subset_idx * 400 + 100, 130), FONT_HERSHEY_SIMPLEX, 1, textColour);
        putText(data, "Z:", Point(subset_idx * 400, 180), FONT_HERSHEY_SIMPLEX, 1, textColour);
        putText(data, sz.str(), Point(subset_idx * 400 + 100, 180), FONT_HERSHEY_SIMPLEX, 1, textColour);
    }
    stringstream time_str;

    putText(data, "cross_time:", Point(0, 230), FONT_HERSHEY_SIMPLEX, 1, textColour);
    putText(data, cross_and_trian_time_str.str(), Point(300, 230), FONT_HERSHEY_SIMPLEX, 1, textColour);

    putText(data, "state_2_timer:", Point(0, 280), FONT_HERSHEY_SIMPLEX, 1, textColour);
    putText(data, state_2_total_time_str.str(), Point(300, 280), FONT_HERSHEY_SIMPLEX, 1, textColour);

    putText(data, "total_corr_loop:", Point(0, 330), FONT_HERSHEY_SIMPLEX, 1, textColour);
    time_str.str("");
    time_str << total_corr_loop.get()->totalElapsedTime();
    putText(data, time_str.str(), Point(300, 330), FONT_HERSHEY_SIMPLEX, 1, textColour);

    putText(data, "corr_time:", Point(0, 380), FONT_HERSHEY_SIMPLEX, 1, textColour);
    time_str.str("");
    time_str << corr_time.get()->totalElapsedTime();
    putText(data, time_str.str(), Point(300, 380), FONT_HERSHEY_SIMPLEX, 1, textColour);

    putText(data, "write_time:", Point(0, 430), FONT_HERSHEY_SIMPLEX, 1, textColour);
    time_str.str("");
    time_str << write_time.get()->totalElapsedTime();
    putText(data, time_str.str(), Point(300, 430), FONT_HERSHEY_SIMPLEX, 1, textColour);

    vconcat(DisplayFrame,data,DisplayFrame);

    x_button = Rect(0,960,400, 20);
    y_button = Rect(400,960,400, 20);
    z_button = Rect(800,960,400, 20);
    Scalar rectangle_color;
    rectangle_color = x_button_clicked == true?Scalar(100,127,255):Scalar(0,255,0);
    rectangle(DisplayFrame,x_button,rectangle_color,1,8,0);
    rectangle_color = y_button_clicked == true?Scalar(100,127,255):Scalar(0,255,0);
    rectangle(DisplayFrame,y_button,rectangle_color,1,8,0);
    rectangle_color = z_button_clicked == true?Scalar(100,127,255):Scalar(0,255,0);
    rectangle(DisplayFrame,z_button,rectangle_color,1,8,0);

    putText(DisplayFrame, "Show X", Point(150, 975), FONT_HERSHEY_SIMPLEX, 0.5, Scalar(0,255,0));
    putText(DisplayFrame, "Show Y", Point(550, 975), FONT_HERSHEY_SIMPLEX, 0.5, Scalar(0,255,0));
    putText(DisplayFrame, "Show Z", Point(950, 975), FONT_HERSHEY_SIMPLEX, 0.5, Scalar(0,255,0));
    imshow("Real Time DIC - Haemish Kyd",DisplayFrame);
    setMouseCallback("Real Time DIC - Haemish Kyd",dataMouseCallBack);
}

void *GetVideo(void *threadid) {
    long tid;
    Mat InputFrame;
    tid = (long)threadid;
    while (MainDataStruct.WriteThreadRunning)
    {
        if (!MainDataStruct.use_arducam)
        {
            cap2 >> frame1; // get a new frame from camera
            cap1 >> frame2;
        }
        else
        {
            if (!cap1.read(InputFrame))
            {
                std::cout << "Capture read error" << std::endl;
                break;
            }
            frame1 = InputFrame(cv::Rect(0, 0, InputFrame.cols / 2, InputFrame.rows));
            frame2 = InputFrame(cv::Rect(InputFrame.cols / 2, 0, InputFrame.cols / 2, InputFrame.rows));
        }
    }
}

void *WriteImageFiles(void *threadid) {
    long tid;
    tid = (long)threadid;
    while (MainDataStruct.WriteThreadRunning) {        
        ReadComplete.wait(tid);
        imwrite("./rdisk_images/Img_0001_0.jpeg", frame1);
        imwrite("./rdisk_images/Img_0001_1.jpeg", frame2);
        WriteComplete.notify(tid);
    }
    pthread_exit(NULL);
}

string gstreamer_pipeline (int capture_width, int capture_height, int display_width, int display_height, int framerate, int flip_method) {    
    return "nvarguscamerasrc ! video/x-raw(memory:NVMM), width=(int)" + std::to_string(capture_width) + ", height=(int)" +
           std::to_string(capture_height) + ", format=(string)NV12, framerate=(fraction)" + std::to_string(framerate) +
           "/1 ! nvvidconv flip-method=" + std::to_string(flip_method) + " ! video/x-raw, width=(int)" + std::to_string(display_width) + ", height=(int)" +
           std::to_string(display_height) + ", format=(string)BGRx ! videoconvert ! video/x-raw, format=(string)BGR ! appsink";
}

int main(int argc, char *argv[]) {
    int return_val;
    int start_calibration = false;
    float Brightness;
    int system_state = 0;
    bool only_cc = false;
    string test_file_name;
    pthread_t threads[1];
    Mat OutputFrame;
    Mat DisplayFrame;
	int capture_width = 4032 ;
    int capture_height = 3040 ;
    int display_width = 2048 ;
    int display_height = 1536 ;    
    int framerate = 13 ;
    int flip_method = 2 ;
    Mat InputFrame;
    
    MainDataStruct.use_arducam = false;
    if (argc == 3){
        if (strcmp(argv[1], "--arducam") == 0)
        {
            MainDataStruct.use_arducam = true;
        }
        if (strcmp(argv[2], "--only_cc") == 0)
        {
            only_cc = true;
            cout << "Only executing cross correlation." << endl;
        }
        else
        {
            read_script(argv[2]);
        }
    }
    else if (argc == 1){
        if (strcmp(argv[1], "--only_cc") == 0)
        {
            only_cc = true;
            cout << "Only executing cross correlation." << endl;
        }
        else
        {
            read_script(argv[1]);
        }
    }

    string port("/dev/ttyUSB0");
    SerialPort serial_port(port);
    try {
        serial_port.Open();
    }
    catch (const SerialPort::OpenFailed &) {
        cout << "Serial Port did not open correctly." << endl;
    }
    if (serial_port.IsOpen()) {
        serial_port.SetBaudRate(SerialPort::BAUD_115200);
        serial_port.SetCharSize(SerialPort::CHAR_SIZE_8);
        serial_port.SetFlowControl(SerialPort::FLOW_CONTROL_NONE);
        serial_port.SetParity(SerialPort::PARITY_NONE);
        serial_port.SetNumOfStopBits(SerialPort::STOP_BITS_1);
    }

    x_button_clicked = false;
    y_button_clicked = false;
    z_button_clicked = false;
    MainDataStruct.WriteThreadRunning = true;

    if (!MainDataStruct.use_arducam)
    {
        cap2.set(CAP_PROP_FRAME_WIDTH,1280);
        cap2.set(CAP_PROP_FRAME_HEIGHT,960);
        cap1.set(CAP_PROP_FRAME_WIDTH,1280);
        cap1.set(CAP_PROP_FRAME_HEIGHT,960);
    }

    Brightness = 0;
    MainDataStruct.FrameWidth = 640;
    MainDataStruct.FrameHeight = 480;

    cout << "====================================" << endl << endl;
    cout << "Default Brightness -------> " << Brightness << endl;
    cout << "Default Width      -------> " << MainDataStruct.FrameWidth << endl;
    cout << "Default Height     -------> " << MainDataStruct.FrameHeight << endl;
    cout << "====================================" << endl;

    if (!MainDataStruct.use_arducam)
    {
        cap1.open(1);
        cap2.open(2);
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
    }
    else{ 
        std::string pipeline = gstreamer_pipeline(capture_width,
        capture_height,
        display_width,
        display_height,
        framerate,
        flip_method);
        std::cout << "Using pipeline: \n\t" << pipeline << "\n";
        
        cap1.open(pipeline, cv::CAP_GSTREAMER);
        
        if(!cap1.isOpened()) {
            std::cout<<"Failed to open camera."<<std::endl;
            return (-1);
        } else {
            cout << "Arducam is open\n";
        }
    }

    //This is here simply to clear the file.
    std::stringstream timeFileName;
    timeFileName << "timing.txt";
    std::ofstream ofs(timeFileName.str(), std::ofstream::out);
    ofs << "Cross Correlation and Triangulations" << "," << "Entire Process" << ","
        << "Correlation Loop" << "," << "Correlation Process" << ","
        << "Write Process" << "," << "Hz" <<endl;
    ofs.close();

    ScriptRun script_obj = (ScriptRun(&MainDataStruct.script_stack, &serial_port, &schema));
    namedWindow("Real Time DIC - Haemish Kyd", WINDOW_AUTOSIZE);
    for (;;) {
        switch (system_state) {
            case 0: {
                if (!MainDataStruct.use_arducam)
                {
                    cap2 >> frame1; // get a new frame from camera
                    cap1 >> frame2;
                }
                else{
                    if (!cap1.read(InputFrame)) {
                        std::cout<<"Capture read error"<<std::endl;
                        break;
                    }
                    frame1 = InputFrame(cv::Rect(0,0,InputFrame.cols/2,InputFrame.rows));
                    frame2 = InputFrame(cv::Rect(InputFrame.cols/2,0,InputFrame.cols/2,InputFrame.rows));
                }				
                /**
                 * If we only do the cross correlation (mapping points to points) we can use
                 * this to trace the accuracy of the points.
                 */
                if (!MainDataStruct.Left_X.empty() || !MainDataStruct.Left_Y.empty()) {
                    for (int num_it = 0; num_it < MainDataStruct.Left_X.size(); num_it++) {
                        Scalar colour_choice = Scalar(255, 0 , 0);
                        circle(frame1, Point(MainDataStruct.Left_X.at(num_it), MainDataStruct.Left_Y.at(num_it)), 3,
                               colour_choice, 2);
                        circle(frame2, Point(MainDataStruct.Right_X.at(num_it), MainDataStruct.Right_Y.at(num_it)), 3,
                               colour_choice, 2);
                    }
                }
                /**
                 * This allows for one user point that allows one to check the accuracy of a
                 * point in the left mapped to a point in the right.
                 */
                if ((MainDataStruct.User_Left_X != 0) || (MainDataStruct.User_Left_Y != 0)){
                    Scalar colour_choice = Scalar(0, 255 , 0);
                    circle(frame1, Point(MainDataStruct.User_Left_X, MainDataStruct.User_Left_Y), 3,
                           colour_choice, 2);
                    circle(frame2, Point(MainDataStruct.User_Right_X, MainDataStruct.User_Right_Y), 3,
                           colour_choice, 2);
                }
                hconcat(frame1, frame2, OutputFrame);
                /**
                 * If we only do the cross correlation (mapping points to points) we can
                 * use this to trace the accuracy of the points.
                 */
                if (!MainDataStruct.Left_X.empty() || !MainDataStruct.Left_Y.empty()) {
                    for (int num_it = 0; num_it < MainDataStruct.Left_X.size(); num_it++) {
                        Scalar colour_choice = Scalar(255, 0 , 0);
                        line(OutputFrame, Point(MainDataStruct.Left_X.at(num_it), MainDataStruct.Left_Y.at(num_it)),
                             Point(MainDataStruct.Right_X.at(num_it) + 640, MainDataStruct.Right_Y.at(num_it)),
                             colour_choice);
                    }
                }
                if ((MainDataStruct.User_Right_X != 0) || (MainDataStruct.User_Right_Y != 0)){
                    Scalar colour_choice = Scalar(0, 255 , 0);
                    line(OutputFrame, Point(MainDataStruct.User_Left_X, MainDataStruct.User_Left_Y),
                         Point(MainDataStruct.User_Right_X + 640, MainDataStruct.User_Right_Y),
                         colour_choice);
                }

                resize(OutputFrame,DisplayFrame,Size(1280,480));
                imshow("Real Time DIC - Haemish Kyd", DisplayFrame);

                if (start_calibration == true){
                    run_calibration(frame1,frame2,&start_calibration);
                }
                }
                break;
            case 1:
                if (!MainDataStruct.use_arducam)
                {
                    cap2 >> frame1; // get a new frame from camera
                    cap1 >> frame2;
                }
                else
                {
                    if (!cap1.read(InputFrame)) {
                        std::cout<<"Capture read error"<<std::endl;
                        break;
                    }
                    frame1 = InputFrame(cv::Rect(0,0,InputFrame.cols/2,InputFrame.rows));
                    frame2 = InputFrame(cv::Rect(InputFrame.cols/2,0,InputFrame.cols/2,InputFrame.rows));
                }
                DICe::initialize(argc, argv);
                information_extraction();
                run_cross_correlation();

                cross_and_trian_time_str.str("");
                cross_and_trian_time_str << cross_time.get()->totalElapsedTime();

                if (only_cc == true){
                    system_state = 0;
                    setMouseCallback("Real Time DIC - Haemish Kyd",projectionMouseCallBack);
                }
                else {
                    (void) pthread_create(&threads[0], NULL, WriteImageFiles, (void *) 0);
                    (void)pthread_create(&threads[0], NULL, GetVideo, (void *)0);
                    ReadComplete.notify(0);
                    system_state = 2;
                }
                break;
            case 2: {
                Teuchos::TimeMonitor state_2_monitor(*state_2_timer);

                WriteComplete.wait(0);
                main_stereo_3d_correlation();

                outputImageInformation();
            }
                break;
                }
        state_2_total_time_str.str("");
        state_2_total_time_str << state_2_timer.get()->totalElapsedTime()<<"("<<int(1.0/state_2_timer.get()->totalElapsedTime())<<" Hz)";

        if (system_state > 0) {
            write_timing_metrics();
        }

        Teuchos::TimeMonitor::zeroOutTimers();

        char c = waitKey(5);

        /**
         * Run the calibration
         */
        if (c == 'c'){
            start_calibration = true;
        }
        /**
         * Quit the program completely
         */
        if (c == 'q') {
            MainDataStruct.WriteThreadRunning = false;
            pthread_cancel(threads[0]);
            break;
        }
        /**
         * Initiate the 3D DIC
         */
        if (c == 'i') {
            imwrite("./rdisk_images/Img_0000_0.jpeg", frame1);
            imwrite("./rdisk_images/Img_0000_1.jpeg", frame2);
            system_state = 1;
        }
        /**
         * Run the pre-loaded script
         */
        if (c == 'r') {
            MainDataStruct.myScript = &script_obj;
        }
        if (MainDataStruct.myScript != NULL){
            if (MainDataStruct.myScript->_script_loaded){
                MainDataStruct.myScript->ExecuteStep();
            }
        }
        if (serial_port.IsOpen()) {
            /**
             * Run the target forward a millimetre
             */
            if (c == '+') {
                serial_port.Write("F0.12\n");
                cout << "Forward" << endl;
            }
            /**
             * Run the target backwards a millimetre
             */
            if (c == '-') {
                serial_port.Write("B0.12\n");
                cout << "Back" << endl;
            }
            /**
             *  Run forwards at a speed
             */
            if (c == 'w')
            {
                serial_port.Write("W200\n");
                cout << "Cont. Forward" << endl;
            }
            /**
             *  Run backwards at a speed
             */
            if (c == 'x')
            {
                serial_port.Write("X200\n");
                cout << "Cont. Back" << endl;
            }
            while (serial_port.IsDataAvailable()) {
                char data_byte;
                // Specify a timeout value (in milliseconds).
                int ms_timeout = 250;

                data_byte = serial_port.ReadByte(ms_timeout);
                cout << data_byte << flush;
            }
        }
    }

    if (!MainDataStruct.use_arducam)
    {
        cap1.release();
        cap2.release();
    }
    else
    {
	    cap1.release();
    }
    cv::destroyAllWindows() ;
    return return_val;
    }
