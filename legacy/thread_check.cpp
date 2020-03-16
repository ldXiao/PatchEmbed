
#include <vector>
#include <thread>
#include <future>
#include <algorithm>
#include <iostream>
#include <random>
#include <chrono>
class join_threads
{
public:
  explicit join_threads(std::vector<std::thread>& threads)
    : threads_(threads) {}

  ~join_threads()
  {
    for (size_t i = 0; i < threads_.size(); ++i)
    {
      if(threads_[i].joinable())
      {
        threads_[i].join();
      }
    }
  }

private:
  std::vector<std::thread>& threads_;
};

template<typename Iterator, typename Func>
void parallel_for_each(Iterator first, Iterator last, Func func)
{
  const auto length = std::distance(first, last);
  if (0 == length) return;

  const auto min_per_thread = 25u;
  const unsigned max_threads = (length + min_per_thread - 1) / min_per_thread;

  const auto hardware_threads = std::thread::hardware_concurrency();

  const auto num_threads= std::min(hardware_threads != 0 ?
        hardware_threads : 2u, max_threads);

  const auto block_size = length / num_threads;

  std::vector<std::future<void>> futures(num_threads - 1);
  std::vector<std::thread> threads(num_threads-1);
  join_threads joiner(threads);

  auto block_start = first;
  for (unsigned i = 0; i < num_threads - 1; ++i)
  {
    auto block_end = block_start;
    std::advance(block_end, block_size);
    std::packaged_task<void (void)> task([block_start, block_end, func]()
    {
      std::for_each(block_start, block_end, func);
    });
    futures[i] = task.get_future();
    threads[i] = std::thread(std::move(task));
    block_start = block_end;
  }

  std::for_each(block_start, last, func);

  for (size_t i = 0; i < num_threads - 1; ++i)
  {
    futures[i].get();
  }
}


using namespace std;

constexpr size_t ARRAY_SIZE = 500'000'000;
typedef std::vector<uint64_t> Array;

template <class FE, class F>
void test_for_each(const Array& a, FE fe, F f, atomic<uint64_t>& result)
{
  auto time_begin = chrono::high_resolution_clock::now();
  result = 0;
  fe(a.begin(), a.end(), f);
  auto time_end = chrono::high_resolution_clock::now();

  cout << "Result = " << result << endl;
  cout << "Time: " << chrono::duration_cast<chrono::milliseconds>(
            time_end - time_begin).count() << endl;
}

int main()
{
  random_device device;
  default_random_engine engine(device());
  uniform_int_distribution<uint8_t> distribution(0, 255);

  Array a;
  a.reserve(ARRAY_SIZE);

  cout << "Generating array ... " << endl;
  for (size_t i = 0; i < ARRAY_SIZE; ++i)
    a.push_back(distribution(engine));

  atomic<uint64_t> result;
  auto acc = [&result](uint64_t value) { result += value; };

  cout << "parallel_for_each ..." << endl;
  test_for_each(a, parallel_for_each<Array::const_iterator, decltype(acc)>, acc, result);
  cout << "for_each ..." << endl;
  test_for_each(a, for_each<Array::const_iterator, decltype(acc)>, acc, result);

  return 0;
}