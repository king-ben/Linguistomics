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

ReadAlignmentTask::ReadAlignmentTask(std::string fn, size_t os, std::vector<std::vector<Alignment*>>* a, std::vector<std::string>* an) :
    ThreadTask(), fileName(fn), offset(os), alignments(a), alignmentNames(an) {

}

ReadAlignmentTask::~ReadAlignmentTask(void) {
}

std::mutex ReadAlignmentTask::coutMutex;

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
        return;
        }

    if (!j.is_array())
        Msg::error("Expected JSON array of sampled alignments in " + fileName);

    {
    std::lock_guard<std::mutex> lock(coutMutex);
    std::cout << "   * Alignment " << offset+1 << ": \"" << fileName << "\"" << std::endl;
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

    // spawn worker threads
    for (int i = 0; i < threadCount; i++)
        threads[i] = std::thread(&ThreadPool::worker, this);
}

ThreadPool::~ThreadPool(void) {

    // wait for all pending work to complete
    wait();
    
    // signal shutdown
    {
        std::lock_guard<std::mutex> lock(mutex);
        running = false;
    }
    
    // wake all workers so they can exit
    taskAvailable.notify_all();
    
    // join all threads
    if (threads)
        {
        for (int i = 0; i < threadCount; i++)
            threads[i].join();
        delete [] threads;
        }
}

void ThreadPool::pushTask(ThreadTask* task) {

    std::lock_guard<std::mutex> lock(mutex);
    
    // check for queue overflow
    if (queueSize >= queueCapacity)
        Msg::error("ThreadPool queue overflow");
    
    // add to circular buffer (fast modulo via bitmask)
    taskQueue[queueTail] = task;
    queueTail = (queueTail + 1) & queueMask;
    ++queueSize;
    ++tasksInFlight;
    
    // wake one waiting worker
    taskAvailable.notify_one();
}

ThreadTask* ThreadPool::popTask(void) {

    std::unique_lock<std::mutex> lock(mutex);
    
    // wait until there's a task OR we're shutting down
    taskAvailable.wait(lock, [this] {
        return queueSize > 0 || !running;
    });
    
    // check shutdown condition
    if (!running && queueSize == 0)
        return nullptr;
    
    // pop from circular buffer
    ThreadTask* task = taskQueue[queueHead];
    queueHead = (queueHead + 1) & queueMask;
    --queueSize;
    
    return task;
}

void ThreadPool::wait(void) {

    std::unique_lock<std::mutex> lock(mutex);
    
    // block until all tasks are complete (no spinning)
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
            // execute outside any lock -- this is where parallelism happens
            task->run();
            
            // decrement in-flight count and potentially wake wait()
            if (--tasksInFlight == 0)
                allComplete.notify_all();
            }
        // no yield() here -- immediately try to get next task
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
    threads(nullptr)
{
}

ThreadPool::~ThreadPool(void) {
}

void ThreadPool::pushTask(ThreadTask* task) {

    // execute immediately in serial mode
    MathCache cache;
    task->run(cache);
}

ThreadTask* ThreadPool::popTask(void) {

    return nullptr;
}

void ThreadPool::wait(void) {

    // Nothing to wait for in serial mode - tasks execute synchronously
}

void ThreadPool::worker(void) {
}

#endif
