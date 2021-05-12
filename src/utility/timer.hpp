#ifndef INCLUDE_TIMER
#define INCLUDE_TIMER

#include <atomic> // std::atomic_thread_fence
#include <chrono> // std::chrono
#include <sstream> // std::stringstream

namespace ppjsdm {

template<typename Clock = std::chrono::high_resolution_clock>
class Timer {
public:
  Timer():
  start_(Clock::now()) {}

  template<typename Rep = typename Clock::duration::rep>
  auto elapsed_time() const {
    std::stringstream ss;
    ss << "Elapsed time: ";
    std::atomic_thread_fence(std::memory_order_relaxed);
    auto dur(Clock::now() - start_);
    std::atomic_thread_fence(std::memory_order_relaxed);
    const auto h(std::chrono::duration_cast<std::chrono::hours>(dur));
    const auto m(std::chrono::duration_cast<std::chrono::minutes>(dur -= h));
    const auto s(std::chrono::duration_cast<std::chrono::seconds>(dur -= m));
    const auto ms(std::chrono::duration_cast<std::chrono::milliseconds>(dur -= s));
    if(h.count() > 0) {
      ss << h.count() << " hours, ";
    }
    if(m.count() > 0) {
      ss << m.count() << " minutes, ";
    }
    if(s.count() > 0) {
      ss << s.count() << " seconds, ";
    }
    ss << ms.count() << " milliseconds.\n";
    return ss.str();
  }

private:
  const typename Clock::time_point start_;
};

using PreciseTimer = Timer<>;
using SystemTimer = Timer<std::chrono::system_clock>;
using MonotonicTimer = Timer<std::chrono::steady_clock>;

} // namespace ppjsdm

#endif // INCLUDE_TIMER
