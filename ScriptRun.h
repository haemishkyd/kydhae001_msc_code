//
// Created by haemish on 2020/07/04.
//

#ifndef CUSTOM_APP_SCRIPTRUN_H
#define CUSTOM_APP_SCRIPTRUN_H

using namespace std;
using namespace std::chrono;

class ScriptRun {
public:
    vector<string> *ScriptStack;
    high_resolution_clock::time_point LastExecutionPoint;
    int CurrentStackPointer;
    bool ScriptLoaded;
    int IteratorCount;

    ScriptRun(vector<string> *p_stack);
    void ExecuteStep(SerialPort *controlPort);
};


#endif //CUSTOM_APP_SCRIPTRUN_H
