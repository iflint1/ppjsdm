#ifndef INCLUDE_ARRAY_SIZE
#define INCLUDE_ARRAY_SIZE

namespace ppjsdm {

template<typename T, std::size_t N>
constexpr std::size_t array_size(T (&)[N]) {
  return N;
}

} // namespace ppjsdm

#endif // INCLUDE_ARRAY_SIZE
