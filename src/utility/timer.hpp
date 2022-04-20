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
  start_(Clock::now()),
  current_(start_) {}

  void set_current() {
    current_ = get_current_time();
  }

  template<typename Unit = std::chrono::seconds>
  auto get_elapsed_time() const {
    return std::chrono::duration_cast<Unit>(get_current_time() - current_);
  }

  template<typename Unit = std::chrono::seconds>
  auto get_total_time() const {
    return std::chrono::duration_cast<Unit>(get_current_time() - start_);
  }

  auto print_elapsed_time() const {
    return make_printing_string(get_current_time() - current_);
  }

  auto print_total_time() const {
    return make_printing_string(get_current_time() - start_);
  }

private:
  const typename Clock::time_point start_;
  typename Clock::time_point current_;

  auto get_current_time() const {
    std::atomic_thread_fence(std::memory_order_relaxed);
    const auto now(Clock::now());
    std::atomic_thread_fence(std::memory_order_relaxed);
    return now;
  }

  auto make_printing_string(decltype(current_ - start_) dur) const {
    std::stringstream ss;
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
};

using PreciseTimer = Timer<>;
using SystemTimer = Timer<std::chrono::system_clock>;
using MonotonicTimer = Timer<std::chrono::steady_clock>;

} // namespace ppjsdm

#endif // INCLUDE_TIMER
