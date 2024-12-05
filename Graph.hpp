
#ifndef _DSL_GRAPH_HPP_
#define _DSL_GRAPH_HPP_

#if defined(__cpp_structured_bindings) && defined(__cpp_variable_templates) && defined(__cpp_if_constexpr)

#include "General.hpp"

#include <vector>
#include <cstddef>
#include <unordered_map>
#include <unordered_set>
#include <limits>
#include <queue>
#include <stdint.h>
#include <functional>

#ifdef __cpp_concepts
#include <concepts>

namespace dsl {
namespace graph {
namespace concepts {

/**
 * Concept to limit index type
 */
template<class T>
concept IndexType = std::regular<T>;

/**
 * Type can be hashed by std::hash
 */
template<class T>
concept Hashable = 
    requires(T t) {
        { std::hash<T>()(t) } -> std::same_as<size_t>;
    }
;

/**
 * Concept to limit index type for matrix storage
 * index type should satisfy IndexType and std::integral
 */
template<class T>
concept MatrixIndexType =
    IndexType<T> &&
    std::integral<T>
;

/**
 * Concept to limit index for hash list storage
 */
template<class T>
concept HashIndex =
    IndexType<T> &&
    Hashable<T>
;

/**
 * Concept for index provider class
 * Implemented by DefaultIndexProvider
 */
template<class T>
concept IndexProvider = 
    std::default_initializable<T> &&
    requires {
        typename T::key_type;
        typename T::index_type;
        typename T::value_type;
    } &&
    requires(
        T t,
        const T ct,
        T::key_type key,
        T::index_type idx,
        T::value_type val,
        size_t c,
        std::vector<typename T::index_type> contain
    ) {
        { t.remove(idx) } ;
        { t.rewind(idx) } ;
        { t.insert(val) } -> std::same_as<typename T::index_type>;
        { t.available() } -> std::same_as<typename T::index_type>;
        { t.findAll(key) } -> std::same_as<std::vector<typename T::index_type>>;
        { t.find(key) } -> std::same_as<typename T::index_type>;
        { t.size() } -> std::same_as<size_t>;
        { t.at(idx) } -> std::same_as<typename T::value_type&>;
        { ct.at(idx) } -> std::same_as<const typename T::value_type&>;
        { t.allIndexes(contain) } ;
        { t.emplace() };
    }
;

/**
 * Concept for storage provider class
 * Implemented by MatrixStorage & HashListStorage
 */
template<class T>
concept StorageProvider = 
    std::default_initializable<T> &&
    requires {
        typename T::index_type;
        typename T::weight_type;
        typename T::storage_type;
        typename T::null_weight;
        requires std::same_as<const typename T::weight_type, decltype(T::fallback)>;
    } &&
    requires(
        T t,
        T::index_type index,
        T::weight_type weight,
        size_t size,
        std::vector<std::pair<
            typename T::index_type,
            typename T::weight_type*
        >> contain
    ) {
        { T::null_weight::value() } -> std::same_as<typename T::weight_type>;
        { t.sync(size) } ;
        { t.addIndex(index) } ;
        { t.removeIndex(index) } -> std::same_as<typename T::index_type>;
        { t.addEdge(index, index, weight) } ;
        { t.removeEdge(index, index) } ;
        { t.getWeight(index, index) } -> std::same_as<const typename T::weight_type&>;
        { t.setWeight(index, index, weight) } ;
        { t.expose() } -> std::same_as<typename T::storage_type*>;
        { t.getForth(index, contain) } ;
        { t.getBack(index, contain) } ;
    }
;

}}}
// namespace dsl::graph::concepts

#define DSL_MACRO_INDEX_TYPE std::regular
#define DSL_MACRO_PREDICATE(...) std::predicate<__VA_ARGS__>
#define DSL_MACRO_INDEX_PROVIDER concepts::IndexProvider
#define DSL_MACRO_STORE_PROVIDER concepts::StorageProvider
#define DSL_MACRO_VALUE_TYPE std::regular
#define DSL_MACRO_MATRIX_INDEX concepts::MatrixIndexType
#define DSL_MACRO_WEIGHT_TYPE std::regular
#define DSL_MACRO_DEFAULT_INDEX std::integral
#define DSL_MACRO_HASH_INDEX concepts::HashIndex

#else

#define DSL_MACRO_INDEX_PROVIDER class
#define DSL_MACRO_STORE_PROVIDER class
#define DSL_MACRO_INDEX_TYPE class
#define DSL_MACRO_PREDICATE(...) class
#define DSL_MACRO_VALUE_TYPE class
#define DSL_MACRO_MATRIX_INDEX class
#define DSL_MACRO_WEIGHT_TYPE class
#define DSL_MACRO_DEFAULT_INDEX class
#define DSL_MACRO_HASH_INDEX class

#endif /* C++20 */

namespace dsl {
namespace graph {

namespace utils {

/**
  Helper struct to provide default weight value.
  Do template specialization to specify default value for a type.
  Rewrite the following function:
  ```cpp
  static constexpr T value() noexcept { ... }
  ```
 */
template<DSL_MACRO_WEIGHT_TYPE T>
struct null_weight {
    static constexpr T value() noexcept { return T(); }
};

// boolean specialization
template<>
struct null_weight<bool> {
    static constexpr bool value() noexcept { return false; }
};

/**
  Helper struct to specify a limit for specific index type.
  Use std::numeric_limits<T> by default.
  Rewrite the following functions:
  ```cpp
  static constexpr T max() noexcept { ... }
  static constexpr T min() noexcept { ... }
  ```
 */
template<DSL_MACRO_INDEX_TYPE T>
struct index_limits {
    static constexpr T max() noexcept { return std::numeric_limits<T>::max(); }
    static constexpr T min() noexcept { return std::numeric_limits<T>::min(); }
};


/**
 * @brief 自定义类型的键选择器
 * @tparam T 给定类型
 * @details 
 * 若图的值类型为自定义类型，可以特化该结构体来指定自定义类型的键以方便部分函数的使用。
 * 特化时需定义的类型如下：
 * -# key_type 键类型标识
 * 特化时需重写的函数如下：
 * -# static constexpr const key_type& key(const T&)
 * -# static constexpr key_type& key(T&)
 * -# static constexpr key_type&& arg(Args&&...)
 */
template<class T>
struct key_selector {
    typedef T key_type;
    static constexpr const key_type& key(const T& t) { return t; }
    static constexpr key_type& key(T& t) { return t; }

    static constexpr key_type&& arg(key_type&& key) {
        return static_cast<key_type&&>(key); 
    }
};

}
// namespace dsl::graph::utils

namespace defines {

// update strategy used in GraphAccessor
enum class UpdateStrategy: uint8_t { forth, both, back, none };

}
// namespace dsl::graph::defines

namespace accessors {

template<
    class _ValTp,
    class _IdxTp,
    class _WhtTp
>
struct AdjPreview {
    typedef _ValTp value_type;
    typedef _IdxTp index_type;
    typedef _WhtTp weight_type;

    index_type index;
    value_type* value_ptr;
    weight_type* weight_ptr;

    AdjPreview(
        const index_type& idx,
        value_type* vp,
        weight_type* wp
    ): index(idx), value_ptr(vp), weight_ptr(wp) { }
};

/**
 * SimpleGraph类使用的访问器（非const）
 * 提供对于访问器所指结点的值的引用访问，
 * 以及所指结点的所有邻接点的值的引用访问。
 */
template<
    class _ValTp,
    class _IdxTp,
    class _WhtTp,
    class _cvIdxProv,
    class _cvStProv
>
class GraphAccessor {
public:
    typedef _ValTp& ref_type;
    typedef std::remove_cv_t<_ValTp> value_type;
    typedef _IdxTp index_type;
    typedef _WhtTp weight_type;
    typedef AdjPreview<_ValTp, _IdxTp, _WhtTp> preview_type;
    typedef std::vector<preview_type> adjacent_list;
    typedef std::function<bool(_ValTp&)> callback_t;

private:
    typedef std::remove_cv_t<_cvIdxProv> index_prov_t;
    typedef std::remove_cv_t<_cvStProv> store_prov_t;

    typedef GraphAccessor<
        _ValTp, _IdxTp, _WhtTp, _cvIdxProv, _cvStProv
    > self;

    _cvIdxProv* index_ptr;
    _cvStProv* storage_ptr;
    index_type index;
    adjacent_list forth_list, back_list;

    void updateForth() {
        forth_list.clear();
        std::vector<std::pair<index_type, weight_type*>> following;
        storage_ptr->getForth(index, following);
        for (auto& pair: following) {
            forth_list.emplace_back(
                pair.first,
                &(index_ptr->at(pair.first)),
                pair.second
            );
        }
    }
    void updateBack() {
        back_list.clear();
        std::vector<std::pair<index_type, weight_type*>> following;
        storage_ptr->getBack(index, following);
        for (auto& pair: following) {
            back_list.emplace_back(
                pair.first,
                &(index_ptr->at(pair.first)),
                pair.second
            );
        }
    }

    void copy_from(const self& gc) {
        index_ptr = gc.index_ptr;
        storage_ptr = gc.storage_ptr;
        index = gc.index;
#ifdef DSL_GRAPH_COPY_FULL_ACCESSOR 
        forth_list = gc.forth_list;
        back_list = gc.back_list;
#endif
    }
    void move_from(self& gc) {
        index_ptr = gc.index_ptr;
        storage_ptr = gc.storage_ptr;
        index = gc.index;
        forth_list = std::move(gc.forth_list);
        back_list = std::move(gc.back_list);
    }

public:
    GraphAccessor(
        _cvIdxProv* i,
        _cvStProv* s,
        const index_type& idx
    ): index_ptr(i), storage_ptr(s), index(idx)
    { }

    GraphAccessor(const self& ga) { copy_from(ga); }
    GraphAccessor(self&& ga) { move_from(ga); }
    self& operator=(const self& ga) { copy_from(ga); return *this; }
    self& operator=(self&& ga) { move_from(ga); return *this; }
    
    inline bool invalid() const {
        return (index_ptr == nullptr || storage_ptr == nullptr);
    }

    /**
     * 返回访问器对应的结点值的引用
     */
    ref_type operator*() const { return index_ptr->at(index); }

    /**
     * 返回访问器对应的内部下标的值
     */
    index_type raw() const { return index; }

    /**
     * 将通过下标移动访问器至指定点
     * （不做连通性检查）
     */
    self& move(
        const index_type& dest,
        defines::UpdateStrategy update = defines::UpdateStrategy::none
    ) {
        index = dest;
        return updateAdjacent(update);
    }

    self next(
        const index_type& dest,
        defines::UpdateStrategy update = defines::UpdateStrategy::none
    ) const {
        return self(index_ptr, storage_ptr, dest).updateAdjacent(update);
    }

    self& updateAdjacent(
        defines::UpdateStrategy how
    ) {
        switch (how) {
        case defines::UpdateStrategy::back :
            updateBack(); break;
        case defines::UpdateStrategy::forth :
            updateForth(); break;
        case defines::UpdateStrategy::both :
            updateForth();
            updateBack();
            break;
        default:
            break;
        }
        return *this;
    }

    /**
     返回由出边连接的邻接点的下标及其值指针的对构成的列表的引用
     内含`struct AdjPreview`
     结构化绑定元素顺序：
     1. index_type
     2. value_type*
     3. weight_type*
    */
    const adjacent_list& listForth() const { return forth_list; }

    /**
     返回由入边连接的邻接点的下标及其值指针的对构成的列表的引用
     内含`struct AdjPreview`
     结构化绑定元素顺序：
     1. index_type
     2. value_type*
     3. weight_type*
    */
    const adjacent_list& listBack() const { return back_list; }
};

}
// namespace dsl::graph::accessor

template<
    DSL_MACRO_VALUE_TYPE _ValTp,
    DSL_MACRO_DEFAULT_INDEX _IdxTp,
    bool _EnableReverseHashBoost = true
>
class DefaultIndexProvider {
public:
    typedef _ValTp value_type;
    typedef _IdxTp index_type;
    typedef utils::key_selector<value_type>::key_type key_type;
private:
    typedef utils::key_selector<value_type> select_key;
    typedef utils::index_limits<_IdxTp> limit;
    typedef std::unordered_map<key_type, index_type> reverse_map_t;
    static constexpr bool enable_rhb = (
        general::utils::is_hashable_v<key_type> &&
        _EnableReverseHashBoost
    );

    index_type present_index;
    std::unordered_map<index_type, value_type> st;
    reverse_map_t* rst_ptr;
public:
#ifdef DSL_DEBUG

    void _show() const {
        for (auto& [index, val]: st) {
            std::cout << index << "=>" << &val << '\n';
        }
        std::cout << "Total: " << st.size() << '\n';
        if constexpr (enable_rhb) {
            std::cout << "Have reverse map:\n";
            for (auto& [f, v]: *rst_ptr) {
                std::cout << f << "=>" << v << '\n';
            }
        }
    }

#endif

    DefaultIndexProvider():
    present_index(limit::min()), st(), rst_ptr(nullptr) { 
        if constexpr (enable_rhb) {
            rst_ptr = new reverse_map_t();
        }
    }
    ~DefaultIndexProvider() {
        if constexpr (enable_rhb) {
            if (rst_ptr != nullptr) delete rst_ptr;
        }
    }

    void allIndexes(std::vector<index_type>& contain) const {
        for (auto& pair: st) contain.emplace_back(pair.first);
    }

    index_type insert(const value_type& node) {
        st.emplace(present_index, node);
        // create reverse map if value is hashable
        if constexpr (enable_rhb) {
            rst_ptr->emplace(
                select_key::key(node),
                present_index
            );
        }
        if (present_index < limit::max() - 1) return present_index++;
        return present_index;
    }
    template<class... Args>
    index_type emplace(Args&&... args) {
        st.emplace(
            std::piecewise_construct, 
            std::forward_as_tuple(present_index),
            std::forward_as_tuple(args...)
        );
        // create reverse map if key is hashable
        if constexpr (enable_rhb) {
            if constexpr (std::is_same_v<key_type, value_type>) {
                rst_ptr->emplace(
                    std::piecewise_construct,
                    std::forward_as_tuple(args...),
                    std::forward_as_tuple(present_index)
                );
            } else {
                rst_ptr->emplace(
                    select_key::arg(std::forward<Args>(args)...),
                    present_index
                );
            }
        }
        if (present_index < limit::max() - 1) return present_index++;
        return present_index;
    }
    void remove(index_type index) { st.erase(index); }

    size_t size() const { return st.size(); }

    std::vector<index_type> findAll(const key_type& key) const {
        std::vector<index_type> results;
        if constexpr (enable_rhb) {
            results.emplace_back(rst_ptr->at(key));
        } else {
            for (
                auto iter = st.cbegin();
                iter != st.cend();
                ++iter
            ) {
                if (select_key::key(iter->second) == key) {
                    results.emplace_back(iter->first);
                }
            }
        }
        return results;
    }

    index_type find(const key_type& key) const {
        if constexpr (enable_rhb) {
            return rst_ptr->at(key);
        } else {
            for (auto iter = st.cbegin(); iter != st.cend(); ++iter) {
                if (select_key::key(iter->second) == key) return iter->first;
            }
            return limit::max();
        }
    }

    const value_type& at(index_type index) const { return st.at(index); }
    value_type& at(index_type index) { return st.at(index); }

    index_type available() const {
        auto beg = st.cbegin();
        if (beg == st.cend()) return limit::max();
        return beg->first;
    }

    void rewind(index_type index) { present_index = index; }

};

template<
    DSL_MACRO_HASH_INDEX _IdxTp,
    DSL_MACRO_WEIGHT_TYPE _WhtTp,
    bool _Directed
>
class HashListStorage {
public:
    typedef _IdxTp index_type;
    typedef _WhtTp weight_type;
    typedef std::unordered_map<
        index_type, std::unordered_map<index_type, weight_type>*
    > storage_type;
    typedef utils::null_weight<weight_type> null_weight;

    static constexpr weight_type fallback = null_weight::value();

private:
    typedef HashListStorage<_IdxTp, _WhtTp, _Directed> self;
    typedef std::vector<std::pair<index_type, weight_type*>> contain_type;
    typedef utils::index_limits<index_type> idx_limit;

    storage_type list;
    size_t edge_count;

    void clear_up() {
        for (auto& pair: list) delete pair.second;
        edge_count = 0;
    }
    void copy_from(const self& h) {
        clear_up();
        for(auto [index, st_ptr]: h.list) {
            list.emplace(
                index,
                new std::unordered_map<index_type, weight_type>(*st_ptr)
            );
        }
        edge_count = h.edge_count;
    }
    void move_from(self& rh) {
        clear_up();
        for (auto [index, st_ptr]: rh.list) {
            list.emplace(index, st_ptr);
        }
        rh.list.clear();
        edge_count = rh.edge_count;
    }

public:
#ifdef DSL_DEBUG

    void _show() const {
        for (auto [index, st_ptr]: list) {
            std::cout << "|*" << index ;
            for (auto [idx, weight]: *st_ptr) {
                std::cout << " ->[" << idx << ": " << weight << ']';
            }
            std::cout << '\n';
        }
        std::cout << "Edge Count: " << edge_count << '\n';
    }

#endif

    HashListStorage(): list(), edge_count(0) {  }
    ~HashListStorage() { clear_up(); }
    
    HashListStorage(const self& h) { copy_from(h); }
    HashListStorage(self&& h) { move_from(h); }
    self& operator= (const self& h) { copy_from(h); return *this; }
    self& operator= (self&& h) { move_from(h); return *this; }

    /**
     * [StorageProvider.expose]
     */
    storage_type* expose() const { return &list; }

    /**
     * [StorageProvider.sync]
     */
    void sync(size_t vex_size) { return ; }

    /**
     * [StorageProvider.addIndex]
     */
    void addIndex(const index_type& idx) {
        auto iter = list.find(idx);
        if (iter != list.end()) return ;
        list.emplace(
            idx,
            new std::unordered_map<index_type, weight_type>()
        );
    }

    /**
     * [StorageProvider.removeIndex]
     */
    index_type removeIndex(const index_type& idx) {
        size_t rm_edge = 0;
        auto iter = list.find(idx);
        if (iter == list.end()) return idx_limit::max();
        rm_edge += iter->second->size();
        list.erase(iter);
        for (auto [index, st_ptr]: list) {
            if constexpr (_Directed) {
                rm_edge += st_ptr->erase(idx);
            } else {
                st_ptr->erase(idx);
            }
        }
        edge_count -= rm_edge;
        return idx_limit::max();
    }

    /**
     * Will not add edge if any of the indexes does not exist.
     * [StorageProvider.addEdge]
     */
    void addEdge(
        const index_type& from,
        const index_type& to,
        const weight_type& weight
    ) {
        auto iter_from = list.find(from);
        auto iter2_to = list.find(to);
        if (iter_from == list.end() || iter2_to == list.end()) 
            return ;
        
        auto adj_ptr = iter_from->second;
        auto iter_to = adj_ptr->find(to);
        if (iter_to == adj_ptr->end()) {
            adj_ptr->emplace(to, weight);
        } else {
            iter_to->second = weight;
        }
        ++edge_count;
        
        if constexpr (!_Directed) {
            auto adj_ptr = iter2_to->second;
            auto iter_from = adj_ptr->find(from);
            if (iter_from == adj_ptr->end()) {
                adj_ptr->emplace(from, weight);
            } else {
                iter_from->second = weight;
            }
        }
    }

    /**
     * [StorageProvider.removeEdge]
     */
    void removeEdge(
        const index_type& from,
        const index_type& to
    ) {
        auto iter_from = list.find(from);
        if (iter_from == list.end()) return ;
        iter_from->second->erase(to);
        --edge_count;
    }

    /**
     * [StorageProvider.getWeight]
     */
    const weight_type& getWeight(
        const index_type& from,
        const index_type& to
    ) const {
        auto iter = list.find(from);
        if (iter == list.cend()) {
            return fallback;
        } else {
            auto iter2 = iter->second->find(to);
            if (iter2 == iter->second->cend()) {
                return fallback;
            } else {
                return iter2->second;
            }
        }
    }

    /**
     * [StorageProvider.setWeight]
     */
    void setWeight(
        const index_type& from,
        const index_type& to,
        const weight_type& weight
    ) {
        auto iter = list.find(from);
        if (iter == list.end()) {
            return ;
        } else {
            auto iter2 = iter->second->find(to);
            if (iter2 == iter->second->end()) {
                return ;
            } else {
                iter2->second = weight;
            }
        }
        if constexpr (!_Directed) {
            auto iter3 = list.find(to);
            auto iter4 = iter3->second->find(from);
            iter4->second = weight;
        }
    }

    /**
     * [StorageProvider.getForth]
     */
    void getForth(
        const index_type& idx,
        contain_type& contain
    ) const {
        auto iter = list.find(idx);
        if (iter == list.cend()) return ;
        for (auto& pair: *(iter->second)) {
            contain.emplace_back(
                pair.first,
                &(pair.second)
            );
        }
    }

    /**
     * [StorageProvider.getBack]
     */
    void getBack(
        const index_type& idx,
        contain_type& contain
    ) const {
        for (auto [index, st_ptr]: list) {
            auto iter = st_ptr->find(idx);
            if (iter == st_ptr->cend()) continue;
            contain.emplace_back(
                index,
                &(iter->second)
            );
        }
    }
};

template<
    DSL_MACRO_MATRIX_INDEX _IdxTp,
    DSL_MACRO_WEIGHT_TYPE _WhtTp,
    bool _Directed
>
class MatrixStorage {
public:
    typedef _IdxTp index_type;
    typedef _WhtTp weight_type;
    typedef std::vector<std::vector<weight_type>*> storage_type;
    typedef utils::null_weight<weight_type> null_weight;

    static constexpr weight_type fallback = null_weight::value();

private:
    typedef MatrixStorage<_IdxTp, _WhtTp, _Directed> self;

    std::vector<std::vector<weight_type>*> matrix;
    size_t vex_size, edge_count;

    void clear_up() {
        for (auto ptr: matrix) { if (ptr != nullptr) delete ptr; }
        matrix.clear();
    }
    void copy_from(const self& mst) {
        clear_up();
        matrix.reserve(mst.matrix.size());
        for (auto ptr: mst.matrix) {
            matrix.push_back(new std::vector<weight_type>(*ptr));
        }
        vex_size = mst.vex_size;
        edge_count = mst.edge_count;
    }
    void move_from(self& mst) {
        clear_up();
        matrix.reserve(mst.matrix.size());
        for (auto ptr: mst.matrix) matrix.push_back(ptr);
        mst.matrix.clear();
        vex_size = mst.vex_size;
        edge_count = mst.edge_count;
    }

public:
/* Debug Functions */
#ifdef DSL_DEBUG

    void _show() const {
        std::cout << "Matrix:\n";
        for (auto ptr: matrix) {
            for (const auto& weight: (*ptr)) std::cout << weight << ' ';
            std::cout << '\n';
        }
        std::cout << "Vex Size: " << vex_size
            << "\tEdge Count: " << edge_count << '\n';
    }

#endif

    MatrixStorage():
        matrix(), vex_size(0), edge_count(0)
    { }
    ~MatrixStorage() { clear_up(); }

    MatrixStorage(const self& ms) { copy_from(ms); }
    MatrixStorage(self&& ms) { move_from(ms); }
    self& operator=(const self& ms) { copy_from(ms); return *this; }
    self& operator=(self&& ms) { move_from(ms); return *this; }

    /**
     * Sync storage structure with vertex count of current graph
     * [StorageProvider.sync]
     */
    void sync(size_t v_size) {
        if (vex_size <= v_size) {
            vex_size = v_size;
#ifdef DSL_DEBUG
            std::cout << "MS: sync -> update size\n";
#endif
        } else {
#ifdef DSL_DEBUG
            std::cout << "MS: sync -> squeeze size\n";
#endif  
            for (size_t i = 0; i < v_size; ++i) {
                for (size_t j = v_size; j < vex_size; ++j) {
                    matrix[i]->at(j) = null_weight::value();
                }
            }
            for (size_t i = v_size; i < vex_size; ++i) {
                for (size_t j = 0; j < vex_size; ++j) {
                    matrix[i]->at(j) = null_weight::value();
                }
            }
            vex_size = v_size;
            return ;
        }
        if (matrix.size() >= v_size) return;
#ifdef DSL_DEBUG
        std::cout << "MS: sync -> extend matrix\n";
#endif
        size_t need = v_size - matrix.size();
        for (std::vector<weight_type>* ptr: matrix) {
            for (size_t i = 0; i < need; ++i)
                ptr->push_back(null_weight::value());
        }
        for (size_t i = 0; i < need; ++i) {
            auto ptr = new std::vector<weight_type>(
                v_size, null_weight::value()
            );
            matrix.push_back(ptr);
        }
    }

    /**
     * [StorageProvider.addNode]
     */
    void addIndex(index_type idx) { }

    /**
     * Add an edge to graph
     * [StorageProvider.addEdge]
     */
    void addEdge(
        index_type from,
        index_type to,
        const weight_type& weight
    ) {
        if (from >= vex_size || to >= vex_size) return ;
        if constexpr (_Directed) {
            matrix[from]->at(to) = weight;
        } else {
            matrix[from]->at(to) = weight;
            matrix[to]->at(from) = weight;
        }
        ++edge_count;
    }

    /**
     * Remove an edge from graph
     * [StorageProvider.removeEdge]
     */
    void removeEdge(
        index_type from,
        index_type to
    ) {
        if (from >= vex_size || to >= vex_size) return ;
        if constexpr (_Directed) {
            matrix[from]->at(to) = null_weight::value();
        } else {
            matrix[from]->at(to) = null_weight::value();
            matrix[to]->at(from) = null_weight::value();
        }
        --edge_count;
    }

    /**
     * Expose storage
     * [StorageProvider.expose]
     */
    storage_type* expose() const {
        return &matrix;
    }

    /**
     * Get weight of specified edge
     * [StorageProvider.getWeight]
     */
    const weight_type& getWeight(
        index_type from,
        index_type to
    ) const {
        if (from >= vex_size || to >= vex_size) return fallback;
        return matrix[from]->at(to);
    }

    /**
     * Set weight of specified edge
     * [StorageProvider.setWeight]
     */
    void setWeight(
        index_type from,
        index_type to,
        const weight_type& weight
    ) const {
        if (from >= vex_size || to >= vex_size) return ;
        if constexpr (_Directed) {
            matrix[from]->at(to) = weight;
        } else {
            matrix[from]->at(to) = weight;
            matrix[to]->at(from) = weight;
        }
    }
    
    /**
     * [StorageProvider.getForth]
     */
    void getForth(
        index_type idx,
        std::vector<std::pair<
            index_type, weight_type*
        >>& contain
    ) const {
        if (idx >= vex_size) return ;
        for (size_t i = 0; i < vex_size; ++i) {
            if (matrix[idx]->at(i) != null_weight::value()) {
                contain.push_back(
                    std::make_pair(
                        i, &(matrix[idx]->at(i))
                    )
                );
            }
        }
    }

    /**
     * [StorageProvider.getBack]
     */
    void getBack(
        index_type idx,
        std::vector<std::pair<
            index_type, weight_type*
        >>& contain
    ) const {
        if (idx >= vex_size) return ;
        for (size_t i = 0; i < vex_size; ++i) {
            if (matrix[i]->at(idx) != null_weight::value()) {
                contain.push_back(
                    std::make_pair(
                        i, &(matrix[i]->at(idx))
                    )
                );
            }
        }
    }

    /**
     * [StorageProvider.removeIndex]
     */
    index_type removeIndex(index_type idx) {
        index_type last = static_cast<index_type>(vex_size - 1);
        size_t rm_edges = 0;
        // count rm_edge
        if constexpr (_Directed) {
            for (size_t i = 0; i < vex_size - 1; ++i) {
                if (matrix[idx]->at(i) != null_weight::value()) ++rm_edges;
                matrix[idx]->at(i) = matrix[last]->at(i);
                if (matrix[i]->at(idx) != null_weight::value()) ++rm_edges;
                matrix[i]->at(idx) = matrix[i]->at(last);
            }
            if (matrix[idx]->at(idx) != null_weight::value()) ++rm_edges;
        } else {
            for (size_t i = 0; i < vex_size - 1; ++i) {
                if (matrix[idx]->at(i) != null_weight::value()) ++rm_edges;
                matrix[idx]->at(i) = matrix[last]->at(i);
            }
            if (matrix[idx]->at(idx) != null_weight::value()) ++rm_edges;
        }
        // remove
        matrix[idx]->at(last) = null_weight::value();
        matrix[last]->at(idx) = null_weight::value();
        // dump
        for (size_t i = 0; i < vex_size - 1; ++i) {
            matrix[idx]->at(i) = matrix[last]->at(i);
            matrix[i]->at(idx) = matrix[i]->at(last);
        }
        matrix[idx]->at(idx) = matrix[last]->at(last);
        edge_count -= rm_edges;
        return last;
    }

};

/*!
 * @brief 
 * @tparam _ValTp 
 * @tparam _WhtTp 
 * @tparam _Directed 
 * @tparam _IdxTp 
 * @tparam _StProv 
 * @tparam _IdxProv 
 */
template<
    DSL_MACRO_VALUE_TYPE _ValTp,
    DSL_MACRO_WEIGHT_TYPE _WhtTp,
    bool _Directed,
    DSL_MACRO_INDEX_TYPE _IdxTp = size_t,
    DSL_MACRO_STORE_PROVIDER _StProv = HashListStorage<_IdxTp, _WhtTp, _Directed>,
    DSL_MACRO_INDEX_PROVIDER _IdxProv = DefaultIndexProvider<_ValTp, _IdxTp>
>
class SimpleGraph {
public:
    typedef _ValTp value_type;
    typedef _IdxTp index_type;
    typedef _WhtTp weight_type;
    typedef utils::null_weight<weight_type> null_weight;
    typedef utils::index_limits<index_type> idx_limit;
    typedef _IdxProv::key_type key_type;

    static constexpr index_type nindex = idx_limit::max();
    static constexpr weight_type nweight = _StProv::fallback;

private:
    typedef _IdxProv index_prov_t;
    typedef _StProv store_prov_t;
    typedef
    SimpleGraph<_ValTp, _WhtTp, _Directed, _IdxTp, _StProv, _IdxProv>
    self;

#ifdef DSL_DEBUG
public:
#else
private:
#endif
    index_prov_t index_provider;
    store_prov_t storage_provider;

private:
    void copy_from(const self& g) {
        index_provider = g.index_provider;
        storage_provider = g.storage_provider;
    }
    void move_from(self& g) {
        index_provider = std::move(g.index_provider);
        storage_provider = std::move(g.storage_provider);
    }

public:
    typedef accessors::GraphAccessor<
        value_type, index_type, weight_type,
        index_prov_t, store_prov_t
    > accessor;

    typedef accessors::GraphAccessor<
        const value_type, index_type, weight_type,
        const index_prov_t, const store_prov_t
    > const_accessor;

    SimpleGraph(){  };
    SimpleGraph(const self& g) { copy_from(g); }
    SimpleGraph(self&& rg) { move_from(rg); }
    self& operator= (const self& g) { copy_from(g); return *this; }
    self& operator= (self&& rg) { move_from(rg); return *this; }

    std::vector<index_type> findIndexes(
        const key_type& key, size_t count = 1
    ) const {
        return index_provider.findAll(key);
    }

    index_type find(
        const key_type& key
    ) const {
        return index_provider.find(key);
    }

    std::vector<index_type> allIndexes() const {
        std::vector<index_type> result;
        index_provider.allIndexes(result);
        return result;
    }

    const value_type& nodeAt(const index_type& index) const {
        return index_provider.at(index);
    }
    value_type& nodeAt(const index_type& index) {
        return index_provider.at(index);
    }

    index_type addNode(const value_type& val) {
        storage_provider.sync(index_provider.size() + 1);
        index_type idx = index_provider.insert(val);
        storage_provider.addIndex(idx);
        return idx;
    }
    template<class... Args>
    index_type emplaceNode(Args&&... args) {
        storage_provider.sync(index_provider.size() + 1);
        index_type idx = index_provider.emplace(std::forward<Args>(args)...);
        storage_provider.addIndex(idx);
        return idx;
    }
    index_type removeNode(const index_type& idx) {
        const index_type& ret = storage_provider.removeIndex(idx);
        if (ret != idx_limit::max()) {
            index_provider.at(idx) = index_provider.at(ret);
            index_provider.remove(ret);
            index_provider.rewind(ret);
        } else {
            index_provider.remove(idx);
        }
        storage_provider.sync(index_provider.size());
        return ret;
    }

    size_t countVertex() const { return index_provider.size(); }
    size_t countEdge() const { return storage_provider.size(); }

    self& addEdge(
        const index_type& from,
        const index_type& to,
        const weight_type& weight = null_weight::value()
    ) {
        if constexpr (std::is_same_v<weight_type, bool>) {
            storage_provider.addEdge(from, to, true);
        } else {
            storage_provider.addEdge(from, to, weight);
        }
        return *this;
    }
    self& addEdgeByKey(
        const key_type& key_from,
        const key_type& key_to,
        const weight_type& weight = null_weight::value()
    ) {
        return addEdge(
            index_provider.find(key_from),
            index_provider.find(key_to),
            weight
        );
    }

    self& removeEdge(
        const index_type& from,
        const index_type& to
    ) {
        storage_provider.removeEdge(from, to);
        return *this;
    }
    self& removeEdgeByKey(
        const key_type& key_from,
        const key_type& key_to
    ) {
        storage_provider.removeEdge(
            index_provider.find(key_from),
            index_provider.find(key_to)
        );
        return *this;
    }

    const weight_type& getWeight(
        const index_type& from,
        const index_type& to
    ) const {
        return storage_provider.getWeight(from, to);
    }
    const weight_type& getWeightByKey(
        const key_type& key_from,
        const key_type& key_to
    ) const {
        return storage_provider.getWeight(
            index_provider.find(key_from),
            index_provider.find(key_to)
        );
    }

    void setWeight(
        const index_type& from,
        const index_type& to,
        const weight_type& weight
    ) {
        storage_provider.setWeight(from, to, weight);
    }
    void setWeightByKey(
        const key_type& key_from,
        const key_type& key_to,
        const weight_type& weight
    ) {
        storage_provider.setWeight(
            index_provider.find(key_from),
            index_provider.find(key_to),
            weight
        );
    }

    value_type& operator[] (const index_type& index) {
        return index_provider.at(index);
    }
    const value_type& operator[] (const index_type& index) const {
        return index_provider.at(index);
    }

    accessor access() {
        index_type idx = index_provider.available();
        if (idx != idx_limit::max()) {
            return accessor(
                &index_provider,
                &storage_provider,
                idx
            );
        } else {
            return accessor(
                nullptr,
                nullptr,
                idx
            );
        }
    }
    accessor access(const index_type& init_idx) {
        if (init_idx != idx_limit::max()) {
            return accessor(
                &index_provider,
                &storage_provider,
                init_idx
            );
        } else {
            return accessor(
                nullptr,
                nullptr,
                init_idx
            );
        }
    }

    const_accessor const_access() const {
        index_type idx = index_provider.available();
        if (idx != idx_limit::max()) {
            return const_accessor(
                &index_provider,
                &storage_provider,
                idx
            );
        } else {
            return const_accessor(
                nullptr,
                nullptr,
                idx
            );
        }
    }
    const_accessor const_access(const index_type& init_idx) const {
        if (init_idx != idx_limit::max()) {
            return const_accessor(
                &index_provider,
                &storage_provider,
                init_idx
            );
        } else {
            return const_accessor(
                nullptr,
                nullptr,
                init_idx
            );
        }
    }

};

namespace algorithms {

template<
    class _ValTp, class _WhtTp, class _IdxTp,
    class _StProv, class _IdxProv
>
void BFS(
    const accessors::GraphAccessor<
        _ValTp, _IdxTp, _WhtTp, _IdxProv, _StProv 
    >& accessor,
    const typename accessors::GraphAccessor<
        _ValTp, _IdxTp, _WhtTp, _IdxProv, _StProv 
    >::callback_t& on_node
) {
    typedef accessors::GraphAccessor<
        _ValTp, _IdxTp, _WhtTp, _IdxProv, _StProv 
    > __accessor;
    std::unordered_set<_IdxTp> visited;
    std::queue<__accessor> q;
    if (accessor.invalid()) return ;
    q.push(accessor);
    visited.insert(accessor.raw());
    while (!q.empty()) {
        __accessor acc = q.front();
        q.pop();
        if (!on_node(*acc)) return ;
        acc.updateAdjacent(defines::UpdateStrategy::forth);
        for (auto [idx, vp, wp]: acc.listForth()) {
            auto iter = visited.find(idx);
            if (iter == visited.end()) {
                q.push(acc.next(idx));
                visited.insert(idx);
            }
        }
    }
}

template<
    class _ValTp, class _WhtTp, class _IdxTp,
    class _StProv, class _IdxProv
>
void DFS(
    const accessors::GraphAccessor<
        _ValTp, _IdxTp, _WhtTp, _IdxProv, _StProv 
    >& accessor,
    const typename accessors::GraphAccessor<
        _ValTp, _IdxTp, _WhtTp, _IdxProv, _StProv 
    >::callback_t& on_node
) {
    typedef accessors::GraphAccessor<
        _ValTp, _IdxTp, _WhtTp, _IdxProv, _StProv 
    > __accessor;
    std::unordered_set<_IdxTp> visited;
    if (accessor.invalid()) return ;
    visited.insert(accessor.raw());

    std::function<void(__accessor)> _rec = 
    [&visited, &on_node, &_rec](__accessor acc) -> void {
        if (!on_node(*acc)) return ;
        acc.updateAdjacent(defines::UpdateStrategy::forth);
        for (auto [idx, vp, wp]: acc.listForth()) {
            auto iter = visited.find(idx);
            if (iter == visited.end()) {
                visited.insert(idx);
                _rec(acc.next(idx));
            }
        }
    };
    _rec(accessor);
}

// TODO Single source shortest path
void Dijkstra() {

}

// TODO Multiple source shortest path
void Floyd() {

}

}
// namespace dsl::graph::algorithms

}}
// namespace dsl::graph

// operators for accessor

// equality compare
template<
    class _ValTpL, class _ValTpR,
    class _IdxProvL, class _IdxProvR,
    class _StProvL, class _StProvR,
    class _IdxTp, class _WhtTp
>
constexpr bool operator== (
    const dsl::graph::accessors::GraphAccessor<
        _ValTpL, _IdxTp, _WhtTp, _IdxProvL, _StProvL
    >& lga,
    const dsl::graph::accessors::GraphAccessor<
        _ValTpR, _IdxTp, _WhtTp, _IdxProvR, _StProvR
    >& rga
) { return lga.raw() == rga.raw(); }
// not equal
template<
    class _ValTpL, class _ValTpR,
    class _IdxProvL, class _IdxProvR,
    class _StProvL, class _StProvR,
    class _IdxTp, class _WhtTp
>
constexpr bool operator!= (
    const dsl::graph::accessors::GraphAccessor<
        _ValTpL, _IdxTp, _WhtTp, _IdxProvL, _StProvL
    >& lga,
    const dsl::graph::accessors::GraphAccessor<
        _ValTpR, _IdxTp, _WhtTp, _IdxProvR, _StProvR
    >& rga
) { return lga.raw() != rga.raw(); }

// cancel macros

#undef DSL_MACRO_INDEX_TYPE
#undef DSL_MACRO_INDEX_PROVIDER
#undef DSL_MACRO_STORE_PROVIDER
#undef DSL_MACRO_PREDICATE
#undef DSL_MACRO_VALUE_TYPE
#undef DSL_MACRO_WEIGHT_TYPE
#undef DSL_MACRO_DEFAULT_INDEX
#undef DSL_MACRO_HASH_INDEX

#endif
#endif /* _DSL_GRAPH_HPP_ */
