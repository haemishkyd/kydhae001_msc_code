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

#if DICE_MPI
#  include <mpi.h>
#endif

using namespace DICe::field_enums;
using namespace DICe;

bool read_input_data_files(Teuchos::RCP<Teuchos::ParameterList> *l_input_params,
                           Teuchos::RCP<Teuchos::ParameterList> *l_correlation_params);
void run_cross_correlation(Teuchos::RCP<DICe::Schema> l_schema, Teuchos::RCP<DICe::Schema> l_stereo_schema,
                           Teuchos::RCP<DICe::Triangulation> *l_triangulation, Teuchos::RCP<Teuchos::ParameterList> *l_input_params,
                           Teuchos::RCP<Teuchos::ParameterList> *l_correlation_params,
                           std::vector<std::string> image_files, std::vector<std::string> stereo_image_files);

Teuchos::RCP<std::ostream> outStream;
/**
 * Start a whole bunch of metrics to see how long stuff takes!!
 */
Teuchos::RCP<Teuchos::Time> total_time = Teuchos::TimeMonitor::getNewCounter("## Total Time ##");
Teuchos::RCP<Teuchos::Time> cross_time = Teuchos::TimeMonitor::getNewCounter("Cross-correlation");
Teuchos::RCP<Teuchos::Time> corr_time = Teuchos::TimeMonitor::getNewCounter("Correlation");
Teuchos::RCP<Teuchos::Time> write_time = Teuchos::TimeMonitor::getNewCounter("Write Output");


bool read_input_data_files(Teuchos::RCP<Teuchos::ParameterList> *l_input_params,
                           Teuchos::RCP<Teuchos::ParameterList> *l_correlation_params) {
    /**
     * Get all of the input parameters from the input files.
     */
    *l_input_params = Teuchos::rcp(new Teuchos::ParameterList());
    Teuchos::Ptr<Teuchos::ParameterList> inputParamsPtr(l_input_params->get());
    Teuchos::updateParametersFromXmlFile("input.xml", inputParamsPtr);
    TEUCHOS_TEST_FOR_EXCEPTION(*l_input_params == Teuchos::null, std::runtime_error, "");

    *outStream << "Input Parameters: " << std::endl;
    (*l_input_params)->print(*outStream);
    *outStream << "\n--- Input read successfully ---\n" << std::endl;

    /**
     * Get all of the correlation parameters from the input files.
     */
    bool is_error_est_run = false;
    if ((*l_input_params)->isParameter(DICe::correlation_parameters_file)) {
        const std::string paramsFileName = (*l_input_params)->get<std::string>(DICe::correlation_parameters_file);
        *l_correlation_params = DICe::read_correlation_params(paramsFileName);
        *outStream << "User specified correlation Parameters: " << std::endl;
        (*l_correlation_params)->print(*outStream);
        is_error_est_run = (*l_correlation_params)->get<bool>(DICe::estimate_resolution_error, false);
        if (is_error_est_run) {
            // force the computing of the image laplacian for the reference image:
            (*l_correlation_params)->set(DICe::compute_laplacian_image, true);
        }
        *outStream << "\n--- Correlation parameters read successfully ---\n" << std::endl;
    } else {
        *outStream << "Correlation parameters not specified by user" << std::endl;
    }

    return is_error_est_run;
}

void run_cross_correlation(Teuchos::RCP<DICe::Schema> *l_schema, Teuchos::RCP<DICe::Schema> *l_stereo_schema, Teuchos::RCP<DICe::Triangulation> *l_triangulation,Teuchos::RCP<Teuchos::ParameterList> *l_input_params,
                           Teuchos::RCP<Teuchos::ParameterList> *l_correlation_params,std::vector<std::string> image_files, std::vector<std::string> stereo_image_files) {
    /* We know this is a stereo analysis so we just assume all is correct */
    Teuchos::TimeMonitor cross_time_monitor(*cross_time);
    TEUCHOS_TEST_FOR_EXCEPTION((*l_schema)->analysis_type() == GLOBAL_DIC, std::runtime_error,
                               "Error, global stereo not enabled yet");
    *outStream << "Processing cross correlation between left and right images" << std::endl;
    (*l_schema)->initialize_cross_correlation(*l_triangulation,
                                           *l_input_params); // images don't need to be loaded by here they get loaded in this routine based on the input params
    (*l_schema)->update_extents(true);
    (*l_schema)->set_ref_image(image_files[0]);
    (*l_schema)->set_def_image(stereo_image_files[0]);
    if ((*l_schema)->use_nonlinear_projection()) {
        (*l_schema)->project_right_image_into_left_frame(*l_triangulation, false);
    }
    (*l_schema)->execute_cross_correlation();
    (*l_schema)->save_cross_correlation_fields();
    *l_stereo_schema = Teuchos::rcp(new DICe::Schema(*l_input_params, *l_correlation_params, *l_schema));
    (*l_stereo_schema)->update_extents();
    (*l_stereo_schema)->set_ref_image(stereo_image_files[0]);
    assert((*l_stereo_schema) != Teuchos::null);
    //if(stereo_schema->use_nonlinear_projection())
    //  stereo_schema->project_right_image_into_left_frame(triangulation,true);
    (*l_stereo_schema)->set_frame_range(0, 2);

    // go ahead and set up the model coordinates field
    (*l_schema)->execute_triangulation(*l_triangulation, *l_stereo_schema);
}

int main(int argc, char *argv[]) {

    try {
        DICe::initialize(argc, argv);
        Teuchos::RCP<Teuchos::ParameterList> input_params;
        Teuchos::RCP<Teuchos::ParameterList> correlation_params;
        std::string output_folder;
        bool is_error_est_run;
        int_t proc_size = 1;
        int_t proc_rank = 0;

        outStream = Teuchos::rcp(&std::cout, false);
#if DICE_MPI
        MPI_Comm_size(MPI_COMM_WORLD,&proc_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&proc_rank);
#endif
        { // scope for total time
            Teuchos::TimeMonitor total_time_monitor(*total_time);
            std::cout << "Start of process." << std::endl;

            if (proc_rank == 0) DEBUG_MSG("Parsing command line options");
            bool force_exit = false;

            /******* Get the input parameters */
            is_error_est_run = read_input_data_files(&input_params, &correlation_params);

            /******* Decipher the image file names (note: zero entry is the reference image) */
            std::vector<std::string> image_files;
            std::vector<std::string> stereo_image_files;
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
            Teuchos::RCP<DICe::Schema> schema = Teuchos::rcp(new DICe::Schema(input_params, correlation_params));
            Teuchos::RCP<DICe::Schema> stereo_schema;
            // let the schema know how many images there are in the sequence and the first frame id:
            schema->set_frame_range(first_frame_id, num_frames);

            /******* Set up the subsets */
            *outStream << "Number of global subsets: " << schema->global_num_subsets() << std::endl;
            for (int_t i = 0; i < schema->local_num_subsets(); ++i) {
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
            Teuchos::RCP<DICe::Triangulation> triangulation;
            if (input_params->isParameter(DICe::calibration_parameters_file) ||
                input_params->isParameter(DICe::camera_system_file)) {
                if (proc_rank == 0)
                    update_legacy_txt_cal_input(
                            input_params); // in case an old txt format cal input file is being used it needs to have width and height added to it
#if DICE_MPI
                MPI_Barrier(MPI_COMM_WORLD);
#endif
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


            run_cross_correlation(&schema, &stereo_schema, &triangulation, &input_params, &correlation_params, image_files, stereo_image_files);

            // iterate through the images and perform the correlation:
            bool failed_step = false;

            for (int_t image_it = 1; image_it <= num_frames; ++image_it) {
                *outStream << "Processing frame: " << image_it << " of " << num_frames << ", " << image_files[image_it]
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
            } // image loop

            schema->write_stats(output_folder, file_prefix);
            if (is_stereo)
                stereo_schema->write_stats(output_folder, stereo_file_prefix);

            if (failed_step)
                *outStream << "\n--- Failed Step Occurred ---\n" << std::endl;
            else
                *outStream << "\n--- Successful Completion ---\n" << std::endl;

            for (int subset_idx = 0; subset_idx < schema->local_num_subsets(); subset_idx++) {
                std::cout << "The DISPLACEMENT_X field value for subset " << subset_idx << " is "
                          << schema->local_field_value(subset_idx, MODEL_DISPLACEMENT_X_FS) << std::endl;
                std::cout << "The DISPLACEMENT_Y field value for subset " << subset_idx << " is "
                          << schema->local_field_value(subset_idx, MODEL_DISPLACEMENT_Y_FS) << std::endl;
                std::cout << "The DISPLACEMENT_Z field value for subset " << subset_idx << " is "
                          << schema->local_field_value(subset_idx, MODEL_DISPLACEMENT_Z_FS) << std::endl;
            }
            // output timing
        } // end scope for total time
        // print the timing data with or without verbose flag
        if (input_params->get<bool>(DICe::print_timing, false)) {
            Teuchos::TimeMonitor::summarize(*outStream, false, true, false/*zero timers*/);
        }
        //  write the time output to file:
        std::stringstream timeFileName;
        timeFileName << output_folder << "timing." << proc_size << "." << proc_rank << ".txt";
        std::ofstream ofs(timeFileName.str(), std::ofstream::out);
        Teuchos::TimeMonitor::summarize(ofs, false, true, false/*zero timers*/);
        ofs.close();
        if (proc_rank != 0) // only keep the process zero copy of the timing results
            std::remove(timeFileName.str().c_str());

        DICe::finalize();

    }
    catch (std::exception &e) {
        std::cout << e.what() << std::endl;
        return 1;
    }


    return 0;
}
