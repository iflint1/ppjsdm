#ifndef INCLUDE_APPROXIMATE_DRAW
#define INCLUDE_APPROXIMATE_DRAW

namespace ppjsdm {
namespace detail {

template<typename Configuration, typename Model>
struct approximate_draw_helper {
  static auto get(const Model&) {
    return Configuration{};
  }
};

} // namespace detail

template<typename Configuration, typename Model>
inline auto approximate_draw(const Model& model) {
  return detail::approximate_draw_helper<Configuration, Model>::get(model);
}

} // namespace ppjsdm

#endif // INCLUDE_APPROXIMATE_DRAW
