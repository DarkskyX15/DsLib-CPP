
#ifndef _DSL_GENERAL_HPP_
#define _DSL_GENERAL_HPP_ 1

namespace dsl {
namespace general {

namespace utils {

template <typename T>
struct is_hashable
{
    template <typename U>
    static constexpr bool __hashable(
        decltype(std::hash<U>()(std::declval<U>()))
    ) { return true; }
    
    template <typename U>
    static constexpr bool __hashable(...) { return false; }

    static constexpr bool value = __hashable<T>(
        std::numeric_limits<size_t>::max()
    );
};

template<typename T>
constexpr bool is_hashable_v = is_hashable<T>::value;

}
// namespace dsl::utils

}}
// namespace dsl::general


#endif /* _DSL_GENERAL_HPP_ */

