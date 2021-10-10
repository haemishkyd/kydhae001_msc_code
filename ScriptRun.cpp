//
// Created by haemish on 2020/07/04.
//
#include <string>
#include <vector>
#include <iostream>
#include <chrono>
#include <iterator>

#include <DICe_ImageIO.h>
#include <DICe_Schema.h>
#include <DICe_Triangulation.h>

#include <sciplot/sciplot.hpp>

#define throw(...)
#include "SerialPort.h"
#undef throw
#include "ScriptRun.h"

using namespace DICe::field_enums;
using namespace DICe;
using namespace std;
using namespace std::chrono;
using namespace sciplot;

Plot plot_obj;

ScriptRun::ScriptRun(vector<string> *p_stack, SerialPort *controlPort, Teuchos::RCP<DICe::Schema> *passedSchema){
    _script_stack = p_stack;
    _serial_port = controlPort;
    _schema = passedSchema;
    CurrentStackPointer = 0;
    X_Axis_Counter = 0;
    _script_loaded = true;
    DataSampleInterval = 0;
    // Change its palette
    plot_obj.palette("dark2");
    plot_obj.legend().show(false);
    cout << "Script Loaded: "<< _script_stack->size() << " lines." << endl;
    for (int p=0; p<_script_stack->size();p++){
        cout << _script_stack->at(p) << endl;
    }
}

void ScriptRun::InitStep(){
    for (int subset_idx = 0; subset_idx < NumberOfSubsets; subset_idx++){
        vector<float> single_subset;
        _y_plot.push_back(single_subset);
    }
}

void ScriptRun::ExecuteStep(){
    static int OldStackPointer=1000000;
    if (_script_loaded == true) {
        if (OldStackPointer != CurrentStackPointer){
            IteratorCount = 0;
            OldStackPointer = CurrentStackPointer;
        }
        cout << "Load Command: " << CurrentStackPointer << "/" << _script_stack->size() << endl;
        string current_command = _script_stack->at(CurrentStackPointer);
        cout << current_command << endl;
        char action = current_command.at(0);
        string data;
        switch (action)
        {
        case 'S':
        case 'F':
        case 'B':
        case 'W':        
            data = current_command.substr(2, 2);
            break;
        case 'G':
        case 'V':
        case 'T':
            data = current_command.substr(2, 4);
            break;
        default:
            data = current_command.substr(2, 2);
            break;
        }

        duration<double, std::milli> time_span = high_resolution_clock::now() - DataSampleTimer;
        if ((DataSampleInterval > 0) && (time_span.count() > DataSampleInterval)){
            _x_plot.push_back(X_Axis_Counter++);
            for (int subset_idx = 0; subset_idx < NumberOfSubsets; subset_idx++){
                switch (WhichAxisToDraw){
                    case 0:
                        _y_plot[subset_idx].push_back((*_schema)->local_field_value(subset_idx, MODEL_DISPLACEMENT_Z_FS));
                        break;
                    case 1:
                        _y_plot[subset_idx].push_back((*_schema)->local_field_value(subset_idx, MODEL_DISPLACEMENT_Y_FS));
                        break;
                    case 2:
                        _y_plot[subset_idx].push_back((*_schema)->local_field_value(subset_idx, MODEL_DISPLACEMENT_X_FS));
                        break;
                }
                
            }
            DataSampleTimer = high_resolution_clock::now();
        }

        if (action == 'S') {
            DataSampleInterval = stoi(data)*100; //stored in 100s of milliseconds in script.
            cout << "Sample Period Set" << endl;
            CurrentStackPointer++;
        }

        if (action == 'G') {
            char command_string[10];
            sprintf(command_string, "W%04d\n", stoi(data));
            _serial_port->Write(command_string);
            LastExecutionPoint = high_resolution_clock::now();
            cout << "Executing G" << endl;
            CurrentStackPointer++;
        }

        if (action == 'V')
        {
            char command_string[10];
            sprintf(command_string, "X%04d\n", stoi(data));
            _serial_port->Write(command_string);
            LastExecutionPoint = high_resolution_clock::now();
            cout << "Executing V" << endl;
            CurrentStackPointer++;
        }

        if (action == 'T')
        {
            char command_string[10];
            sprintf(command_string, "T%04d\n", stoi(data));
            _serial_port->Write(command_string);
            LastExecutionPoint = high_resolution_clock::now();
            cout << "Executing T" << endl;
            CurrentStackPointer++;
        }

        if (action == 'F')
        {
            _serial_port->Write("F0.12\n");
            IteratorCount++;
            if (IteratorCount >= stoi(data))
            {
                LastExecutionPoint = high_resolution_clock::now();
                cout << "Executing F" << endl;
                CurrentStackPointer++;
            }
        }

        if (action == 'B') {
            _serial_port->Write("B0.12\n");
            IteratorCount++;
            if (IteratorCount >= stoi(data)) {
                LastExecutionPoint = high_resolution_clock::now();
                cout << "Executing B" << endl;
                CurrentStackPointer++;
            }
        }

        if (action == 'W') {
            duration<double, std::milli> time_span = high_resolution_clock::now() - LastExecutionPoint;
            if (time_span.count() > 1000*stoi(data)){
                LastExecutionPoint = high_resolution_clock::now();
                cout << "Executing W" << endl;
                CurrentStackPointer++;
            }
        }

        if (current_command == "DRAW") {            
            for (int subset_idx = 0; subset_idx < NumberOfSubsets; subset_idx++){
                plot_obj.drawCurve(_x_plot, _y_plot[subset_idx]).label("Subset");
            }            
            plot_obj.show();
            plot_obj.save("script_output.pdf");
            cout << "Drawing Result" << endl;
            CurrentStackPointer++;
        }

        if (current_command == "END") {
            cout << "Executing END" << endl;
            _script_loaded = false;
        }
    }
}
