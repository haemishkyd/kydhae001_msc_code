//
// Created by haemish on 2020/07/04.
//
#include <string>
#include <vector>
#include <iostream>
#include <chrono>
#include <sstream>
#include <iterator>
#include "SerialPort.h"
#include "ScriptRun.h"

using namespace std;
using namespace std::chrono;

ScriptRun::ScriptRun(vector<string> *p_stack){
    ScriptStack = p_stack;
    CurrentStackPointer = 0;
    ScriptLoaded = true;
}

void ScriptRun::ExecuteStep(SerialPort *controlPort){
    static int OldStackPointer=1000000;
    if (ScriptLoaded == true) {
        if (OldStackPointer != CurrentStackPointer){
            IteratorCount = 0;
            OldStackPointer = CurrentStackPointer;
        }
        string current_command = ScriptStack->at(CurrentStackPointer);
        char action = current_command.at(0);
        string data = current_command.substr(2,2);
        cout << current_command << endl;

        if (action == 'F') {
            controlPort->Write("F0.17\n");
            IteratorCount++;
            if (IteratorCount >= stoi(data)) {
                LastExecutionPoint = high_resolution_clock::now();
                cout << "Executing F" << endl;
                CurrentStackPointer++;
            }
        }

        if (action == 'B') {
            controlPort->Write("B0.17\n");
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

        if (current_command == "END") {
            cout << "Executing END" << endl;
            ScriptLoaded = false;
        }
    }
}
