#pragma once

#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/interprocess/sync/interprocess_semaphore.hpp>
#include <future>



class ThreadPool {

    public:
        static ThreadPool& threadPoolInstance(void) {
        
            static ThreadPool singleThreadPool;
            return singleThreadPool;
        }

        ThreadPool(const ThreadPool& t) = delete;
        
        template <typename F>
        void post(F&& f) {
        
            queue.wait();
            boost::asio::post(*workers, [this, f = std::forward<F>(f)]
                {
                f();
                queue.post();
                });
        }

        template <typename F>
        auto submit(F&& f) -> std::future<decltype(f())> {
        
            queue.wait();
            std::promise<decltype(f())> promise;
            auto future = promise.get_future();
            boost::asio::post(*workers, [this, promise = std::move(promise), f = std::forward<F>(f)] () mutable
                {
                promise.set_value(f());
                queue.post();
                });
            return future;
        }

        void print(void) {
        
            std::cout << "Thread pool has " << threads << " threads:" << std::endl;
        }
        
        void wait(void) {
        
            workers->wait();
        }

        void reset(void) {
        
            if(workers)
                {
                wait();
                workers.reset();
                }
            workers = std::make_unique<boost::asio::thread_pool>(threads);
        }

    private:
        ThreadPool(void) : threads(std::thread::hardware_concurrency()), queue(std::thread::hardware_concurrency()) {
        
            reset();
        }

        size_t const threads;
        std::unique_ptr<boost::asio::thread_pool> workers;
        boost::interprocess::interprocess_semaphore queue;
};
