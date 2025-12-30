#include <fstream>
#include <iostream>
#include "Alignment.hpp"
#include "json.hpp"
#include "Threads.hpp"
#include "Msg.hpp"



ThreadTask::ThreadTask(void) {

}

ThreadTask::~ThreadTask(void) {

}

ReadAlignmentTask::ReadAlignmentTask(std::string fn, size_t os, 
                                     std::vector<std::vector<Alignment*>>* a, 
                                     std::vector<std::string>* an,
                                     std::atomic<int>* counter) :
    ThreadTask(), 
    fileName(fn), 
    offset(os), 
    alignments(a), 
    alignmentNames(an), 
    completedCount(counter) {

}

ReadAlignmentTask::~ReadAlignmentTask(void) {

}

void ReadAlignmentTask::readAlignment(void) {

    std::ifstream ifs(fileName);
    nlohmann::json j;
    try
        {
        j = nlohmann::json::parse(ifs);
        }
    catch (nlohmann::json::parse_error& ex)
        {
        Msg::warning("Error parsing JSON file " + fileName + " at byte " + std::to_string(ex.byte));
        (*completedCount)++;
        return;
        }

    if (!j.is_array())
        {
        Msg::error("Expected JSON array of sampled alignments in " + fileName);
        (*completedCount)++;
        return;
        }
    
    std::vector<Alignment*> fileAlignments;
    bool nameRecorded = false;
    
    for (const auto& sample : j)
        {
        std::string name = sample["Name"];
        if (nameRecorded == false)
            {
            (*alignmentNames)[offset] = name;
            nameRecorded = true;
            }
        Alignment* aln = new Alignment(sample["Data"]);
        fileAlignments.push_back(aln);
        }
    (*alignments)[offset] = fileAlignments;
    
    (*completedCount)++;
}



#if 1
// Threaded version

ThreadPool::ThreadPool(void) :
    queueHead(0),
    queueTail(0),
    queueSize(0),
    threadCount(std::thread::hardware_concurrency()),
    tasksInFlight(0),
    running(true),
    threads(new std::thread[threadCount]) {

    for (int i = 0; i < threadCount; i++)
        threads[i] = std::thread(&ThreadPool::worker, this);
}

ThreadPool::~ThreadPool(void) {

    wait();
    
    {
        std::lock_guard<std::mutex> lock(mutex);
        running = false;
    }
    
    taskAvailable.notify_all();
    
    if (threads)
        {
        for (int i = 0; i < threadCount; i++)
            threads[i].join();
        delete [] threads;
        }
}

void ThreadPool::pushTask(ThreadTask* task) {

    std::lock_guard<std::mutex> lock(mutex);
    
    if (queueSize >= queueCapacity)
        Msg::error("ThreadPool queue overflow");
    
    taskQueue[queueTail] = task;
    queueTail = (queueTail + 1) & queueMask;
    ++queueSize;
    ++tasksInFlight;
    
    taskAvailable.notify_one();
}

ThreadTask* ThreadPool::popTask(void) {

    std::unique_lock<std::mutex> lock(mutex);
    
    taskAvailable.wait(lock, [this] {
        return queueSize > 0 || !running;
    });
    
    if (!running && queueSize == 0)
        return nullptr;
    
    ThreadTask* task = taskQueue[queueHead];
    queueHead = (queueHead + 1) & queueMask;
    --queueSize;
    
    return task;
}

void ThreadPool::wait(void) {

    std::unique_lock<std::mutex> lock(mutex);
    
    allComplete.wait(lock, [this] {
        return tasksInFlight == 0;
    });
}

void ThreadPool::worker(void) {
    
    while (running)
        {
        ThreadTask* task = popTask();
        
        if (task)
            {
            task->run();
            
            if (--tasksInFlight == 0)
                allComplete.notify_all();
            }
        }
}

#else
// Serial version (for debugging or single-threaded builds)

ThreadPool::ThreadPool(void):
    queueHead(0),
    queueTail(0),
    queueSize(0),
    threadCount(1),
    tasksInFlight(0),
    running(false),
    threads(nullptr) {

}

ThreadPool::~ThreadPool(void) {

}

void ThreadPool::pushTask(ThreadTask* task) {

    task->run();
}

ThreadTask* ThreadPool::popTask(void) {

    return nullptr;
}

void ThreadPool::wait(void) {

}

void ThreadPool::worker(void) {

}

#endif
