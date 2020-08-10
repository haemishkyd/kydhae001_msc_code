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

#include "SerialPort.h"
#include "ScriptRun.h"

using namespace DICe::field_enums;
using namespace DICe;
using namespace std;
using namespace std::chrono;
using namespace sciplot;

plot plot_obj;

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
    cout << "Script Loaded: "<< _script_stack->size() << " lines." << endl;
    for (int p=0; p<_script_stack->size();p++){
        cout << _script_stack->at(p) << endl;
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
        string data = current_command.substr(2,2);

        duration<double, std::milli> time_span = high_resolution_clock::now() - DataSampleTimer;
        if ((DataSampleInterval > 0) && (time_span.count() > DataSampleInterval)){
            _x_plot.push_back(X_Axis_Counter++);
            _y_plot.push_back((*_schema)->local_field_value(2, MODEL_DISPLACEMENT_Z_FS));
            DataSampleTimer = high_resolution_clock::now();
        }

        if (action == 'S') {
            DataSampleInterval = stoi(data)*100; //stored in 100s of milliseconds in script.
            cout << "Sample Period Set" << endl;
            CurrentStackPointer++;
        }
        if (action == 'F') {
            _serial_port->Write("F0.12\n");
            IteratorCount++;
            if (IteratorCount >= stoi(data)) {
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
            plot_obj.draw(_x_plot,_y_plot);
            plot_obj.show();
            cout << "Drawing Result" << endl;
            CurrentStackPointer++;
        }

        if (current_command == "END") {
            cout << "Executing END" << endl;
            _script_loaded = false;
        }
    }
}
