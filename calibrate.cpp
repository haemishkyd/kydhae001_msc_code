#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>

// DICe includes for classes used below:
#include <DICe.h>
#include <DICe_Schema.h>
#include <DICe_Calibration.h>

#include "opencv2/opencv.hpp"
#include "calibrate.h"

using namespace DICe::field_enums;
using namespace DICe;
using namespace std;
using namespace cv;

#define NUM_FRAMES_FOR_CAL 16

class Rt
{
public:
    double r11;
    double r12;
    double r13;
    double r21;
    double r22;
    double r23;
    double r31;
    double r32;
    double r33;
    double t1;
    double t2;
    double t3;
    double cx_0;
    double cy_0;
    double cx_1;
    double cy_1;
    double p1_0;
    double p2_0;
    double p1_1;
    double p2_1;
};
Rt myRt;
int calibration_state = 0;
int captured_frame_count = 0;
std::chrono::time_point<std::chrono::high_resolution_clock> capture_timer;

//Mat frame1, frame2;
//
//int main(int argc, char *argv[]) {
//    int start_calibration = false;
//    int calibration_state = 0;
//    int captured_frame_count = 0;
//    int test_image_frame_count = 0;
//    float Brightness;
//    float FrameWidth;
//    float FrameHeight;
//    Mat DisplayFrame;
//
//    std::chrono::time_point<std::chrono::high_resolution_clock> main_timer;
//    std::chrono::time_point<std::chrono::high_resolution_clock> capture_timer;
//
//
//    VideoCapture cap2(0); // open the default camera
//    VideoCapture cap1(2); // open the default camera
//
//    cap2.set(CAP_PROP_FRAME_WIDTH,1280);
//    cap2.set(CAP_PROP_FRAME_HEIGHT,960);
//    cap1.set(CAP_PROP_FRAME_WIDTH,1280);
//    cap1.set(CAP_PROP_FRAME_HEIGHT,960);
//
//    Brightness = cap1.get(CV_CAP_PROP_BRIGHTNESS);
//    FrameWidth = cap1.get(CV_CAP_PROP_FRAME_WIDTH);
//    FrameHeight = cap1.get(CV_CAP_PROP_FRAME_HEIGHT);
//
//    cout << "====================================" << endl << endl;
//    cout << "Default Brightness -------> " << Brightness << endl;
//    cout << "Default Width      -------> " << FrameWidth << endl;
//    cout << "Default Height     -------> " << FrameHeight << endl;
//    cout << "====================================" << endl;
//
//    if (!cap1.isOpened()) {  // check if we succeeded
//        std::cout << "First camera cannot be found\n";
//        return -1;
//    } else {
//        cout << "Camera 1 is open\n";
//    }
//    if (!cap2.isOpened()) {  // check if we succeeded
//        std::cout << "Second camera cannot be found\n";
//        return -1;
//    } else {
//        cout << "Camera 2 is open\n";
//    }
//
//    namedWindow("Real Time DIC Camera Calibrate - Haemish Kyd", WINDOW_AUTOSIZE);
//
//    main_timer = std::chrono::high_resolution_clock::now();
//    for (;;) {
//        int ti = 0;
//        Mat OutputFrame;
//
//        cap2 >> frame1; // get a new frame from camera
//        cap1 >> frame2;
//
//        hconcat(frame1, frame2, OutputFrame);
//        resize(OutputFrame,DisplayFrame,Size(1280,480));
//        imshow("Real Time DIC Camera Calibrate - Haemish Kyd", DisplayFrame);
//
//        if (start_calibration == true) {
//            switch (calibration_state) {
//                case 0: {
//                    std::stringstream ss;
//                    ss << std::setw(4) << std::setfill('0') << captured_frame_count;
//                    std::string s = ss.str();
//
//                    string filename;
//                    std::cout << "Capturing frame: " << captured_frame_count << "\n";
//                    filename = "Cal_" + s + "_0.jpeg"; //Left suffix
//                    imwrite("cal_files/" + filename, frame1); //Frame 1 Left
//                    filename = "Cal_" + s + "_1.jpeg"; //Right suffix
//                    imwrite("cal_files/" + filename, frame2); //Frame 2 Right
//                    captured_frame_count++;
//                    if (captured_frame_count >= NUM_FRAMES_FOR_CAL) {
//                        captured_frame_count = 0;
//                        calibration_state = 2;
//                    } else {
//                        calibration_state = 1;
//                        std::cout << "Waiting for next calibration frame.....\n";
//                    }
//                }
//                    break;
//                case 1:
//                    //wait for next frame signal
//                    break;
//                case 2: {
//                    DICe::initialize(argc, argv);
//                    DICe::Calibration cal("cal_input_dots.xml");
//                    // create the intersection points and generate a calibrated camera system
//                    scalar_t rms = 0.0;
//                    Teuchos::RCP<Camera_System> cam_sys = cal.calibrate(rms);
//                    // write the calibration parameters to file which includes the intersection points
//                    cam_sys->write_camera_system_file("cal.xml");
//                    std::cout << *cam_sys.get() << std::endl;
//                    std::cout << "\nDICe_Cal complete. RMS Error:" << rms << std::endl;
//                    const Matrix<scalar_t,4> *cam_data = (cam_sys.get()->camera(1)).get()->world_cam_trans_matrix();
//                    std::cout << cam_data->cols() << "x" << cam_data->rows() << std::endl;
//                    myRt.r11 = cam_data->operator()(0,0);
//                    myRt.r12 = cam_data->operator()(0,1);
//                    myRt.r13 = cam_data->operator()(0,2);
//                    myRt.r21 = cam_data->operator()(1,0);
//                    myRt.r22 = cam_data->operator()(1,1);
//                    myRt.r23 = cam_data->operator()(1,2);
//                    myRt.r31 = cam_data->operator()(2,0);
//                    myRt.r32 = cam_data->operator()(2,1);
//                    myRt.r33 = cam_data->operator()(2,2);
//                    myRt.t1 = cam_data->operator()(0,3);
//                    myRt.t2 = cam_data->operator()(1,3);
//                    myRt.t3 = cam_data->operator()(2,3);
//                    myRt.cx_0 = (cam_sys.get()->camera(0)).get()->intrinsics()->at(0);
//                    myRt.cy_0 = (cam_sys.get()->camera(0)).get()->intrinsics()->at(1);
//                    myRt.cx_1 = (cam_sys.get()->camera(1)).get()->intrinsics()->at(0);
//                    myRt.cy_1 = (cam_sys.get()->camera(1)).get()->intrinsics()->at(1);
//                    myRt.p1_0 = (cam_sys.get()->camera(0)).get()->intrinsics()->at(11);
//                    myRt.p2_0 = (cam_sys.get()->camera(0)).get()->intrinsics()->at(12);
//                    myRt.p1_1 = (cam_sys.get()->camera(1)).get()->intrinsics()->at(11);
//                    myRt.p2_1 = (cam_sys.get()->camera(1)).get()->intrinsics()->at(12);
//                    std::cout << "R 0:" << myRt.r11 << " " << myRt.r12  << " " << myRt.r13  << " " << myRt.t1 << std::endl;
//                    std::cout << "R 1:" << myRt.r21 << " " << myRt.r22  << " " << myRt.r23  << " " << myRt.t2 << std::endl;
//                    std::cout << "R 2:" << myRt.r31 << " " << myRt.r32  << " " << myRt.r33  << " " << myRt.t3 << std::endl;
//                    std::cout << "C 0:" << myRt.cx_0 << " " << myRt.cy_0 << std::endl;
//                    std::cout << "C 1:" << myRt.cx_1 << " " << myRt.cy_1 << std::endl;
//                    std::cout << "P_Cam 0:" << myRt.p1_0 << " " << myRt.p2_0 << std::endl;
//                    std::cout << "P_Cam 1:" << myRt.p1_1 << " " << myRt.p2_1 << std::endl;
//
//                    for (int int_idx = 0; int_idx<18; int_idx++){
//                        std::cout << "Item "<<int_idx;
//                        std::cout << " Camera 1: "<<(cam_sys.get()->camera(0)).get()->intrinsics()->at(int_idx);
//                        std::cout << " Camera 2: "<<(cam_sys.get()->camera(1)).get()->intrinsics()->at(int_idx) << std::endl;
//                    }
//                    start_calibration = false;
//                    calibration_state = 0;
//                }
//                    break;
//                default: {
//
//                }
//                    break;
//            }
//            main_timer = std::chrono::high_resolution_clock::now();
//        }
//        char c = waitKey(5);
//
//        if (c == 'q') {
//            break;
//        }
//        if (c == 'c') {
//            start_calibration = true;
//            calibration_state = 0;
//            std::cout << "Calibration Started\n";
//            capture_timer = std::chrono::high_resolution_clock::now();
//        }
//        else if (c == 'n') {
//            calibration_state = 0;
//        }
//        else if (c == 'j') {
//            start_calibration = true;
//            calibration_state = 2;
//            std::cout << "Calibration Started With No Capture\n";
//        }
//        if ((start_calibration == true) &&
//            (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now()-capture_timer).count() > 7) &&
//            (calibration_state != 2)){
//            calibration_state = 0;
//            capture_timer = std::chrono::high_resolution_clock::now();
//        }
//        if (c == 's') {
//            string filename;
//            std::stringstream ss;
//            ss << std::setw(4) << std::setfill('0') << test_image_frame_count;
//            std::string s = ss.str();
//
//            filename = "Img_" + s + "_0.jpeg"; //Left suffix
//            imwrite("test_files/" + filename, frame1); //Frame 1 Left
//            filename = "Img_" + s + "_1.jpeg"; //Right suffix
//            imwrite("test_files/" + filename, frame2); //Frame 2 Right
//            test_image_frame_count++;
//        }
//    }
//    return 0;
//}

int run_calibration(Mat frame1, Mat frame2, int *step_calibration, int *start_calibration)
{
    switch (calibration_state)
    {
    case 0:
        if (*step_calibration)
        {
            *step_calibration = false;

            std::stringstream ss;
            ss << std::setw(4) << std::setfill('0') << captured_frame_count;
            std::string s = ss.str();

            string filename;
            std::cout << "Capturing frame: " << captured_frame_count << "\n";
            filename = "Cal_" + s + "_0.jpeg";        //Left suffix
            imwrite("cal_files/" + filename, frame1); //Frame 1 Left
            filename = "Cal_" + s + "_1.jpeg";        //Right suffix
            imwrite("cal_files/" + filename, frame2); //Frame 2 Right
            captured_frame_count++;
            capture_timer = std::chrono::high_resolution_clock::now();
            if (captured_frame_count >= NUM_FRAMES_FOR_CAL)
            {
                captured_frame_count = 0;
                calibration_state = 1;
            }
            else
            {
                std::cout << "Waiting for next calibration frame.....\n";
            }
        }
        break;
    case 1:
        {
            DICe::initialize(0, 0);
            //DICe::Calibration cal("cal_input_dots.xml");
            DICe::Calibration cal("cal_input_checker_board.xml");
            // create the intersection points and generate a calibrated camera system
            scalar_t rms = 0.0;
            Teuchos::RCP<Camera_System> cam_sys = cal.calibrate(rms);
            // write the calibration parameters to file which includes the intersection points
            cam_sys->write_camera_system_file("cal.xml");
            std::cout << *cam_sys.get() << std::endl;
            std::cout << "\nDICe_Cal complete. RMS Error:" << rms << std::endl;
            const Matrix<scalar_t, 4> *cam_data = (cam_sys.get()->camera(1)).get()->world_cam_trans_matrix();
            std::cout << cam_data->cols() << "x" << cam_data->rows() << std::endl;
            myRt.r11 = cam_data->operator()(0, 0);
            myRt.r12 = cam_data->operator()(0, 1);
            myRt.r13 = cam_data->operator()(0, 2);
            myRt.r21 = cam_data->operator()(1, 0);
            myRt.r22 = cam_data->operator()(1, 1);
            myRt.r23 = cam_data->operator()(1, 2);
            myRt.r31 = cam_data->operator()(2, 0);
            myRt.r32 = cam_data->operator()(2, 1);
            myRt.r33 = cam_data->operator()(2, 2);
            myRt.t1 = cam_data->operator()(0, 3);
            myRt.t2 = cam_data->operator()(1, 3);
            myRt.t3 = cam_data->operator()(2, 3);
            myRt.cx_0 = (cam_sys.get()->camera(0)).get()->intrinsics()->at(0);
            myRt.cy_0 = (cam_sys.get()->camera(0)).get()->intrinsics()->at(1);
            myRt.cx_1 = (cam_sys.get()->camera(1)).get()->intrinsics()->at(0);
            myRt.cy_1 = (cam_sys.get()->camera(1)).get()->intrinsics()->at(1);
            myRt.p1_0 = (cam_sys.get()->camera(0)).get()->intrinsics()->at(11);
            myRt.p2_0 = (cam_sys.get()->camera(0)).get()->intrinsics()->at(12);
            myRt.p1_1 = (cam_sys.get()->camera(1)).get()->intrinsics()->at(11);
            myRt.p2_1 = (cam_sys.get()->camera(1)).get()->intrinsics()->at(12);
            std::cout << "R 0:" << myRt.r11 << " " << myRt.r12 << " " << myRt.r13 << " " << myRt.t1 << std::endl;
            std::cout << "R 1:" << myRt.r21 << " " << myRt.r22 << " " << myRt.r23 << " " << myRt.t2 << std::endl;
            std::cout << "R 2:" << myRt.r31 << " " << myRt.r32 << " " << myRt.r33 << " " << myRt.t3 << std::endl;
            std::cout << "C 0:" << myRt.cx_0 << " " << myRt.cy_0 << std::endl;
            std::cout << "C 1:" << myRt.cx_1 << " " << myRt.cy_1 << std::endl;
            std::cout << "P_Cam 0:" << myRt.p1_0 << " " << myRt.p2_0 << std::endl;
            std::cout << "P_Cam 1:" << myRt.p1_1 << " " << myRt.p2_1 << std::endl;

            for (int int_idx = 0; int_idx < 18; int_idx++)
            {
                std::cout << "Item " << int_idx;
                std::cout << " Camera 1: " << (cam_sys.get()->camera(0)).get()->intrinsics()->at(int_idx);
                std::cout << " Camera 2: " << (cam_sys.get()->camera(1)).get()->intrinsics()->at(int_idx) << std::endl;
            }
            (*start_calibration) = false;
            calibration_state = 0;
        }
        break;
    default:
        break;
    }
    // if (((*start_calibration) == true) &&
    //     (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - capture_timer).count() > 30) &&
    //     (calibration_state != 2))
    // {
    //     calibration_state = 0;
    // }
}
