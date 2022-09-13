#ifndef INCLUDE_APPROXIMATE_DRAW
#define INCLUDE_APPROXIMATE_DRAW

namespace ppjsdm {
namespace detail {

template<typename Configuration, typename Generator, typename Model>
struct approximate_draw_helper {
  static auto get(Generator& generator, const Model&) {
    return Configuration{};
  }
};

} // namespace detail

template<typename Configuration, typename Generator, typename Model>
inline auto approximate_draw(Generator& generator, const Model& model) {
  return detail::approximate_draw_helper<Configuration, Generator, Model>::get(generator, model);
}

} // namespace ppjsdm

#endif // INCLUDE_APPROXIMATE_DRAW
