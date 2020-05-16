//
// Created by haemish on 2020/05/16.
//

#ifndef CUSTOM_APP_SEMAPHORE_H
#define CUSTOM_APP_SEMAPHORE_H
#include <mutex>
#include <condition_variable>
#include <iostream>

class Semaphore {
public:
    Semaphore (int count_ = 0)
            : count(count_)
    {
    }

    inline void notify( int tid ) {
        std::unique_lock<std::mutex> lock(mtx);
        count++;
        std::cout << "thread " << tid <<  " notify" << std::endl;
        //notify the waiting thread
        cv.notify_one();
    }
    inline void wait( int tid ) {
        std::unique_lock<std::mutex> lock(mtx);
        while(count == 0) {
            std::cout << "thread " << tid << " wait" << std::endl;
            //wait on the mutex until notify is called
            cv.wait(lock);
            std::cout << "thread " << tid << " run" << std::endl;
        }
        count--;
    }
private:
    std::mutex mtx;
    std::condition_variable cv;
    int count;
};


#endif //CUSTOM_APP_SEMAPHORE_H
