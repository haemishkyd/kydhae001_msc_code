//
// Created by haemish on 2020/07/04.
//

#ifndef CUSTOM_APP_SCRIPTRUN_H
#define CUSTOM_APP_SCRIPTRUN_H


using namespace std;
using namespace std::chrono;

class ScriptRun {
public:
    vector<string> *_script_stack;
    SerialPort *_serial_port;
    Teuchos::RCP<DICe::Schema> *_schema;
    std::vector<float> _x_plot;
    std::vector<vector<float>> _y_plot;
    bool _script_loaded;
    int IteratorCount;
    int X_Axis_Counter;
    high_resolution_clock::time_point LastExecutionPoint;
    high_resolution_clock::time_point DataSampleTimer;
    int DataSampleInterval;
    int CurrentStackPointer;
    int NumberOfSubsets;
    int WhichAxisToDraw;

    ScriptRun(vector<string> *p_stack, SerialPort *controlPort, Teuchos::RCP<DICe::Schema> *passedSchema);
    void ExecuteStep();
    void InitStep();

private:
};


#endif //CUSTOM_APP_SCRIPTRUN_H
