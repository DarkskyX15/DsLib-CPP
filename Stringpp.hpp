
#ifndef _DSL_STRINGPP_HPP_
#define _DSL_STRINGPP_HPP_

#if __cplusplus >= 202002L

#include <cstddef>
#include <stdint.h>
#include <iterator>
#include <vector>
#include <limits>
#include <iostream>
#include <memory>
#include <functional>
#include <unordered_map>

namespace dsl {

namespace concepts {

template<class T>
concept StringLike = 
    std::regular<T> &&
    requires (T str) {
        // typedef requirements
        typename T::iterator;
        typename T::const_iterator;
        typename T::value_type;
        // hash specialization
        typename std::hash<typename T::value_type>;
        // random access iterator
        requires std::random_access_iterator<typename T::iterator>;
        // char like value_type
        requires std::equality_comparable<typename T::value_type>;
        requires std::totally_ordered<typename T::value_type>;
        // member limitations
        { str.begin() } -> std::same_as<typename T::iterator> ;
        { str.end() } -> std::same_as<typename T::iterator> ;
        { str.cbegin() } -> std::same_as<typename T::const_iterator>;
        { str.cend() } -> std::same_as<typename T::const_iterator>;
        { T::npos - 1 } -> std::same_as<size_t>;
    } && requires (const typename T::value_type chr) {
        // value_type hash compatible
        { std::hash<typename T::value_type>()(chr) } -> std::same_as<size_t>;
    }
;

template<class T>
concept DStringCompatibleChar = 
    std::copyable<T> && std::equality_comparable<T> &&
    std::totally_ordered<T> &&
    requires (const char* &p, const T chr, std::string& s) {
        { T::eat(p) } -> std::same_as<T> ;
        { chr.byteLength() } -> std::same_as<size_t>;
        { chr.isSpace() } -> std::same_as<bool>;
        { chr.put(s) } ;
        { chr.show() } ;
    }
;

} // namespace concepts

enum class MatchAlgo: uint8_t {
    BruteForce, KR,
    KMP, BMHBNFS_T, Sunday,
    BMHBNFS_S
};

template<concepts::StringLike S>
class MatchPattern {
public:
    virtual size_t match(const S& str, size_t offset = 0) const = 0;
    virtual std::vector<size_t> matchAll(const S& str, bool digest) const = 0;
    virtual ~MatchPattern() { }
};

namespace inner {

template<concepts::StringLike S>
class BruteForceMatch: public MatchPattern<S> {
private:
    S pattern;
public:
    BruteForceMatch(const S& pat): pattern(pat) { }
    size_t match(const S& str, size_t offset = 0) const override {
        typename S::const_iterator 
            front = str.cbegin(),
            beg1 = str.cbegin(), end1 = str.cend(),
            beg2 = pattern.cbegin(), end2 = pattern.cend();
        size_t pat_size = end2 - beg2;

        beg1 += offset;
        if (beg1 >= end1) return S::npos;
        for (auto iter = beg2; beg1 != end1; ++beg1, iter = beg2) {
            if ( (end1 - beg1) < pat_size ) break;
            bool found = true;
            for (auto copy = beg1; copy != end1, iter != end2; ++copy, ++iter) {
                if (*copy != *iter) { found = false; break; }
            }
            if (found) return beg1 - front;
        }
        return S::npos;
    }
    std::vector<size_t> matchAll(const S& str, bool digest) const override {
        std::vector<size_t> results;
        size_t pos = 0, pat_size = pattern.cend() - pattern.cbegin();
        if (!digest) pat_size = 1;
        pos = match(str);
        while (pos != S::npos) {
            results.push_back(pos);
            pos = match(str, pos + pat_size);
        }
        return results;
    }
};

template<concepts::StringLike S>
class KMPMatch: public MatchPattern<S> {
private:
    S pattern;
    std::vector<size_t> next;
public:
    KMPMatch(const S& pat): pattern(pat), next() {
        typename S::const_iterator beg1 = pat.cbegin();
        size_t m = pat.cend() - pat.cbegin();
        next.resize(m);
        next[0] = S::npos;
        for (size_t i = 1; i < m - 1; ++i) {
            auto j = next[i - 1];
            while (j != S::npos && (*(beg1 + i) != *(beg1 + j + 1))) {
                j = next[j];
            }
            if (*(beg1 + i) == *(beg1 + j + 1)) {
                next[i] = j + 1;
            } else {
                next[i] = S::npos;
            }
        }
    }
    size_t match(const S& str, size_t offset = 0) const override {
        size_t size_pat = pattern.cend() - pattern.cbegin();
        size_t size_str = str.cend() - str.cbegin();
        if (size_pat + offset > size_str) return S::npos;

        size_t pos = S::npos;
        typename S::const_iterator
            front = pattern.cbegin(),
            beg1 = pattern.cbegin(), end1 = pattern.cend(),
            beg2 = str.cbegin(), end2 = str.cend();
        beg2 += offset;
        
        size_t i = 0, j = offset;
        while (beg1 != end1 && beg2 != end2) {
            if (*beg1 == *beg2) {
                ++i; ++j; ++beg1; ++beg2;
            } else if (i > 0) {
                i = next[i - 1] + 1;
                beg1 = front + i;
            } else {
                ++beg2; ++j;
            }
        }
        if (i == size_pat) pos = j - size_pat;
        return pos;
    }
    std::vector<size_t> matchAll(const S& str, bool digest) const override {
        std::vector<size_t> results;
        size_t pos = 0, pat_size = pattern.cend() - pattern.cbegin();
        if (!digest) pat_size = 1;
        pos = match(str);
        while (pos != S::npos) {
            results.push_back(pos);
            pos = match(str, pos + pat_size);
        }
        return results;
    }
};

template<concepts::StringLike S>
class SundayMatch: public MatchPattern<S> {
private:
    S pattern;
    std::unordered_map<typename S::value_type, size_t> offset_table;
public:
    SundayMatch(const S& pat): pattern(pat), offset_table() {
        typename S::const_iterator
            beg = pat.cbegin(), end = pat.cend();
        size_t pat_size = end - beg;
        size_t index = 0;
        for (; beg != end; ++beg, ++index) {
            offset_table[*beg] = pat_size - index;
        }
    }
    size_t match(const S& str, size_t offset = 0) const override {
        typename S::const_iterator
            beg1 = pattern.cbegin(), end1 = pattern.cend(), iter = beg1,
            beg2 = str.cbegin(), end2 = str.cend(), p = beg2;
        auto not_found = offset_table.cend();
        size_t pat_size = end1 - beg1, str_size = end2 - beg2;
        if (str_size < pat_size + offset) return S::npos;
        beg2 += offset; end2 -= pat_size - 1;

        for (; beg2 < end2; ) {
            iter = beg1; p = beg2;
            for (; (iter != end1) && (*iter == *p); ++iter, ++p) ;
            if (iter == end1) {
                return beg2 - str.cbegin();
            } else {
                auto offset = offset_table.find(*(beg2 + pat_size));
                if (offset != not_found) beg2 += offset->second;
                else beg2 += pat_size + 1;
            }
        }
        return S::npos;
    }
    std::vector<size_t> matchAll(const S& str, bool digest) const override {
        std::vector<size_t> results;
        size_t pos = 0, pat_size = pattern.cend() - pattern.cbegin();
        if (!digest) pat_size = 1;
        pos = match(str);
        while (pos != S::npos) {
            results.push_back(pos);
            pos = match(str, pos + pat_size);
        }
        return results;
    }
};

template<concepts::StringLike S>
class BMMatch: public MatchPattern<S> {
private:
    S pattern;
    std::unordered_map<typename S::value_type, size_t> bm_bc;
    size_t k;

    void build() {
        typename S::const_iterator
            front = pattern.cbegin(),
            beg = pattern.cbegin(), end = pattern.cend();
        size_t pat_size = end - beg;

        // build k
        std::vector<size_t> pi(pat_size, 0);
        size_t index = 1; ++beg;
        for (; beg != end; ++beg, ++index) {
            size_t j = pi[index - 1];
            while (j > 0 && *beg != *(front + j)) j = pi[j - 1];
            if (*beg == *(front + j)) ++j;
            pi[index] = j;
        }
        k = pat_size - pi.back();

        // build bm_bc
        beg = pattern.cbegin();
        size_t loop_size = pat_size - 1;
        auto last = *(end - 1);
        bool flag = false;
        for (; loop_size > 0; ++beg, --loop_size) {
            bm_bc[*beg] = loop_size;
            if (!flag && last == *beg) flag = true;
        }
        if (!flag) bm_bc[last] = pat_size;
    }

public:
    BMMatch(const S& pat):
    pattern(pat), bm_bc(), k(0) { build(); }

    size_t match(const S& str, size_t offset = 0) const override {
        typename S::const_iterator
            pat_begin = pattern.cbegin(), pat_end = pattern.cend(),
            str_begin = str.cbegin(), str_end = str.cend();
        auto bc_end = bm_bc.cend();
        size_t str_size = str_end - str_begin;
        size_t pat_size = pat_end - pat_begin;
        size_t pat_last_pos = pat_size - 1;
        // Invalid Offset
        if (pat_size + offset > str_size) return S::npos;
        str_begin += offset;

        auto str_iter = str_begin + pat_last_pos;
        for (; str_iter < str_end; ) {
            // Check last char
            auto pat_iter = pat_begin + pat_last_pos;
            auto str_last_pos = str_iter;
            if (*pat_iter == *str_iter) {
                size_t loop_size = pat_size - 1;
                bool found = true;
                --pat_iter; --str_iter;
                for (; loop_size > 0; --loop_size, --pat_iter, --str_iter) {
                    if (*pat_iter != *str_iter) {
                        found = false; break;
                    }
                }
                if (found) return (++str_iter) - str.cbegin();
            }

            if (str_end - str_last_pos <= 1) return S::npos;
            
            auto bc_iter = bm_bc.find(*(str_last_pos + 1));
            str_iter = str_last_pos;
            if (bc_iter == bc_end) {
                str_iter += pat_size + 1; // Sunday
            } else {
                // Horspool
                auto m_iter = bm_bc.find(*str_last_pos);
                if (m_iter == bc_end) str_iter += pat_size;
                else str_iter += m_iter->second;
            }
        }
        return S::npos;
    }

    std::vector<size_t> matchAll(const S& str, bool digest) const override {
        std::vector<size_t> results;
        typename S::const_iterator
            pat_begin = pattern.cbegin(), pat_end = pattern.cend(),
            str_begin = str.cbegin(), str_end = str.cend();
        auto bc_end = bm_bc.cend();
        size_t str_size = str_end - str_begin;
        size_t pat_size = pat_end - pat_begin;
        size_t pat_last_pos = pat_size - 1;
        size_t offset_0 = pat_last_pos;
        size_t offset = offset_0;

        auto str_iter = str_begin + pat_last_pos;
        for (; str_iter < str_end; ) {
            // Check last char
            auto pat_iter = pat_begin + pat_last_pos;
            auto str_last_pos = str_iter;
            if (*pat_iter == *str_iter) {
                size_t loop_size = offset;
                bool found = true;
                --pat_iter; --str_iter;
                for (; loop_size > 0; --loop_size, --pat_iter, --str_iter) {
                    if (*pat_iter != *str_iter) {
                        found = false; break;
                    }
                }
                if (found) {
                    results.push_back(
                        (str_last_pos - str_begin) - pat_last_pos
                    );
                    offset = k - 1;
                    str_iter = str_last_pos;
                    str_iter += (digest ? pat_size : k);
                    continue;
                }
            }

            if (str_end - str_last_pos <= 1) return results;
            offset = pat_last_pos;
            
            auto bc_iter = bm_bc.find(*(str_last_pos + 1));
            str_iter = str_last_pos;
            if (bc_iter == bc_end) {
                str_iter += pat_size + 1; // Sunday
            } else {
                // Horspool
                auto m_iter = bm_bc.find(*str_last_pos);
                if (m_iter == bc_end) str_iter += pat_size;
                else str_iter += m_iter->second;
            }
        }
        return results;
    }
};

template<concepts::StringLike S>
class KRMatch: public MatchPattern<S> {
private:
    S pattern;
    size_t pattern_hash;
    size_t code_size, mod;
    size_t front_scaler;

    inline static size_t p_mod(int64_t a, int64_t b) {
        return a >= 0 ? a % b : ((a * -1 / b) + 1) * b + a;
    }
public:
    KRMatch(
        const S& pat,
        size_t __code_size = 0x10FFFF,
        size_t __mod = 10000019
    ): pattern(pat), code_size(__code_size), mod(__mod) {
        typename S::const_iterator
            beg = pat.cbegin(), end = pat.cend();
        size_t pat_size = (end - beg) - 1, base = code_size % mod;
        front_scaler = 1;
        while (pat_size > 0) {
            if (pat_size & 1)
                front_scaler = (front_scaler * base) % mod;
            pat_size >>= 1;
            base = (base * base) % mod;
        }
        pattern_hash = 0;
        auto chr_hasher = std::hash<typename S::value_type>();
        for (; beg != end; ++beg) {
            pattern_hash *= code_size;
            pattern_hash %= mod;
            pattern_hash += chr_hasher(*beg);
            pattern_hash %= mod;
        }
    }
    size_t match(const S& str, size_t offset = 0) const override {
        typename S::const_iterator
            beg1 = pattern.cbegin(), end1 = pattern.cend(),
            beg2 = str.cbegin(), end2 = str.cend();
        size_t pat_size = end1 - beg1, str_size = end2 - beg2;
        if (str_size < pat_size + offset) return S::npos;
        beg2 += offset; end2 -= pat_size - 1;
        auto chr_hasher = std::hash<typename S::value_type>();
        
        auto iter = beg2;
        size_t str_hash = 0;
        int64_t hash_temp = 0;
        for (size_t i = pat_size; i > 0; --i, ++iter) {
            str_hash *= code_size;
            str_hash %= mod;
            str_hash += chr_hasher(*iter);
            str_hash %= mod;
        }

        for (; beg2 != end2; ++beg2) {
            if (str_hash == pattern_hash) {
                iter = beg2;
                for (
                    auto pat_iter = beg1;
                    pat_iter != end1;
                    ++pat_iter, ++iter
                ) {
                    if (*pat_iter != *iter) break;
                    else return beg2 - str.cbegin();
                }
            } else {
                hash_temp = str_hash;
                hash_temp -= (chr_hasher(*beg2) * front_scaler) % mod;
                hash_temp *= code_size;
                hash_temp += chr_hasher(*(beg2 + pat_size));
                str_hash = p_mod(hash_temp, mod);
            }
        }
        return S::npos;
    }
    std::vector<size_t> matchAll(const S& str, bool digest) const override {
        std::vector<size_t> results;
        size_t pos = 0, pat_size = pattern.cend() - pattern.cbegin();
        if (!digest) pat_size = 1;
        pos = match(str);
        while (pos != S::npos) {
            results.push_back(pos);
            pos = match(str, pos + pat_size);
        }
        return results;
    }
};

template<concepts::StringLike S>
class InvalidMatch: public MatchPattern<S> {
public:
    InvalidMatch() = default;
    size_t match(const S& str, size_t offset = 0) const override {
        return S::npos;
    }
    std::vector<size_t> matchAll(const S& str, bool digest) const override {
        return std::vector<size_t>();
    }
};

template<concepts::StringLike S>
MatchPattern<S>* rawPatternFactory(
    const S& pattern,
    MatchAlgo algo = MatchAlgo::BruteForce
) {
    if (pattern.cend() - pattern.cbegin() <= 0)
        return new InvalidMatch<S>();
    switch (algo) {
    case MatchAlgo::BruteForce :
        return new BruteForceMatch<S>(pattern);
    case MatchAlgo::KMP :
        return new KMPMatch<S>(pattern);
    case MatchAlgo::Sunday :
        return new SundayMatch<S>(pattern);
    case MatchAlgo::BMHBNFS_T :
        return new BMMatch<S>(pattern);
    case MatchAlgo::KR :
        return new KRMatch<S>(pattern);
    default:
        return new BruteForceMatch<S>(pattern);
    }
}

} // namespace dsl::inner

template<concepts::StringLike S>
std::shared_ptr<MatchPattern<S>> patternFactory(
    const S& pattern,
    MatchAlgo algo = MatchAlgo::BruteForce
) {
    if (pattern.cend() - pattern.cbegin() <= 0)
        return std::make_shared<inner::InvalidMatch<S>>();
    switch (algo) {
    case MatchAlgo::BruteForce :
        return std::make_shared<inner::BruteForceMatch<S>>(pattern);
    case MatchAlgo::KMP :
        return std::make_shared<inner::KMPMatch<S>>(pattern);
    case MatchAlgo::Sunday :
        return std::make_shared<inner::SundayMatch<S>>(pattern);
    case MatchAlgo::BMHBNFS_T :
        return std::make_shared<inner::BMMatch<S>>(pattern);
    case MatchAlgo::KR :
        return std::make_shared<inner::KRMatch<S>>(pattern);
    default:
        return std::make_shared<inner::BruteForceMatch<S>>(pattern);
    }
}

template<concepts::StringLike S>
size_t stringMatch(
    const S& str,
    const S& pat,
    size_t offset = 0,
    MatchAlgo algo = MatchAlgo::BruteForce
) {
    MatchPattern<S>* mat = inner::rawPatternFactory(pat, algo);
    auto res = mat->match(str, offset);
    delete mat; return res;
}

template<concepts::StringLike S>
std::vector<size_t> stringMatchAll(
    const S& str,
    const S& pat,
    bool digest = true,
    MatchAlgo algo = MatchAlgo::BruteForce
) {
    MatchPattern<S>* mat = inner::rawPatternFactory(pat, algo);
    auto res = mat->matchAll(str, digest);
    delete mat; return res;
}

class UTF8Char{
    
    friend std::ostream& operator<<(
        std::ostream& out,
        const UTF8Char& chr
    ) {
        chr.show(); return out;
    }

private:
    uint32_t data;
    inline static bool onBit(char n, char pos) { return n & (1 << pos); }

    /// @brief private default construction
    UTF8Char() = default;
public:
    ~UTF8Char() = default;

    static UTF8Char eat(const char* &ptr) {
        auto chr = UTF8Char();
        if (*ptr == 0) return chr;
        if (*ptr < 0) {
            auto t = *ptr;
            char size = 6;
            while (onBit(t, size)) --size;
            size = 8 - size;
            uint32_t uc = (static_cast<uint8_t>(t << size)) >> size; ++ptr;
            for (char i = 1; i < size - 1; ++i, ++ptr) {
                uc <<= 6;
                uc |= (*ptr) & 0b111111;
            }
            uint32_t ps = size - 2;
            ps <<= 30; uc |= ps;
            chr.data = uc;
        } else {
            chr.data |= *ptr; 
            ++ptr;
        }
        return chr;
    }

    UTF8Char(const UTF8Char& chr) = default;
    UTF8Char& operator= (const UTF8Char& chr) = default;

    UTF8Char(UTF8Char&& r_chr) = default;
    UTF8Char& operator= (UTF8Char&& r_chr) = default;

    size_t byteLength() const { return ((data >> 30) & 0b11) + 1; }
    int32_t codePoint() const { return data & 0x3FFFFFFF; }
    void put(std::string& s) const {
        uint32_t uc = codePoint();
        char size = byteLength();
        if (size == 1) {
            s.push_back(uc);
            return ;
        }
        uint32_t cc = 0;
        switch (size) {
        case 2:
            cc = (0b110 << 5);
            cc |= (0x7C0 & uc) >> 6;
            break;
        case 3:
            cc = (0b1110 << 4);
            cc |= (0xF000 & uc) >> 12;
            break;
        case 4:
            cc = (0b11110 << 3);
            cc |= (0x1C0000 & uc) >> 18;
            break;
        }
        s.push_back(cc);
        cc = (0b10 << 6);
        if (size >= 4) s.push_back(cc | ((0x3F000 & uc) >> 12));
        if (size >= 3) s.push_back(cc | ((0xFC0 & uc) >> 6));
        if (size >= 2) s.push_back(cc | (0x3F & uc));
    }
    void show() const {
        uint32_t uc = codePoint();
        char size = byteLength();
        if (size == 1) {
            std::cout.put(uc);
            return ;
        }
        uint32_t cc = 0;
        switch (size) {
        case 2:
            cc = (0b110 << 5);
            cc |= (0x7C0 & uc) >> 6;
            break;
        case 3:
            cc = (0b1110 << 4);
            cc |= (0xF000 & uc) >> 12;
            break;
        case 4:
            cc = (0b11110 << 3);
            cc |= (0x1C0000 & uc) >> 18;
            break;
        }
        std::cout.put(cc);
        cc = (0b10 << 6);
        if (size >= 4) std::cout.put(cc | ((0x3F000 & uc) >> 12));
        if (size >= 3) std::cout.put(cc | ((0xFC0 & uc) >> 6));
        if (size >= 2) std::cout.put(cc | (0x3F & uc));
    }

    bool operator==(const UTF8Char& chr) const { return data == chr.data; }
    bool operator!=(const UTF8Char& chr) const { return !operator==(chr); }
    bool operator<(const UTF8Char& chr) const { return codePoint() < chr.codePoint(); }
    bool operator>=(const UTF8Char& chr) const { return !operator<(chr); }
    bool operator>(const UTF8Char& chr) const { return codePoint() > chr.codePoint(); }
    bool operator<=(const UTF8Char& chr) const { return !operator>(chr); }

    bool isSpace() const {
        return 
            (data == '\t')  ||
            (data == ' ')   ||
            (data == '\n')  ||
            (data == '\r')  ;
    }
};

template<concepts::DStringCompatibleChar CharType>
class DynamicString {

    friend std::ostream& operator<<(
        std::ostream& out,
        const DynamicString<CharType>& str
    ) {
        for (auto chr : str.storage) chr.show();
        return out;
    }

private:
    typedef DynamicString<CharType> self;

    std::vector<CharType> storage;
    size_t byte_size;

    constexpr static auto v_is_space = 
    [](const CharType& chr) -> bool { return chr.isSpace(); };

    // assignments

    void move_from(self&& r_str) {
        storage = std::move(r_str.storage);
        byte_size = r_str.byte_size;
    }
    void copy_from(const self& str) {
        storage = str.storage;
        byte_size = str.byte_size;
    }

    // comparision

    enum class CompareResult: uint8_t { more, equal, less };
    CompareResult compare(const self& str) const {
        auto beg1 = storage.begin(), end1 = storage.end();
        auto beg2 = str.storage.begin(), end2 = str.storage.end();
        while (beg1 != end1 && beg2 != end2) {
            if (*beg1 < *beg2) return CompareResult::less;
            if (*beg1 > *beg2) return CompareResult::more;
            ++beg1; ++beg2;
        }
        if (beg1 == end1 && beg2 == end2) return CompareResult::equal;
        if (beg1 == end1) return CompareResult::less;
        return CompareResult::more;
    }

public:
    // types

    typedef CharType value_type;
    typedef std::vector<CharType>::iterator iterator;
    typedef std::vector<CharType>::const_iterator const_iterator;

    constexpr static size_t npos = std::numeric_limits<size_t>::max();

    // default construction
    DynamicString(): storage(), byte_size(0) {}

    // construct from c string
    explicit
    DynamicString(const char* c_str):
    storage(), byte_size(0) {
        auto ptr = c_str;
        while (*ptr != '\000') {
            storage.push_back(CharType::eat(ptr));
            byte_size += storage.back().byteLength();
        }
    }

    // construct from c++ string
    explicit
    DynamicString(const std::string& str): DynamicString(str.c_str()) { }

    // copy assignment
    DynamicString(const self& str):
    storage(), byte_size(0) { copy_from(str); }
    self& operator=(const self& str) {
        copy_from(str); return *this;
    }
    // move assignment
    DynamicString(self&& r_str):
    storage(), byte_size(0) { move_from(std::move(r_str)); }
    self& operator=(self&& r_str) {
        move_from(std::move(r_str));
        return *this;
    }

    // states

    size_t byteSize() const { return byte_size; }
    size_t size() const { return storage.size(); }

    // match

    size_t find(
        const self& pat,
        size_t offset = 0, 
        MatchAlgo algo = MatchAlgo::BruteForce
    ) const {
        return stringMatch<self>(*this, pat, offset, algo);
    }

    std::vector<size_t> findAll(
        const self& pat,
        bool digest = true,
        MatchAlgo algo = MatchAlgo::BruteForce
    ) const {
        return stringMatchAll<self>(*this, pat, digest, algo);
    }

    // functions

    const CharType& at(size_t index) const { return storage.at(index); }
    const CharType& operator[] (size_t index) const { return storage.at(index); }

    self substr(size_t start, size_t end = npos, size_t step = 1) const {
        end = (end == npos ? storage.size() : end);
        self copy;
        if (end <= start) return copy;
        for (size_t index = start; index < end; index += step) {
            copy.storage.push_back(storage[index]);
            copy.byte_size += storage[index].byteLength();
        }
        return copy;
    }

    std::vector<self> split(
        std::function<bool(const CharType& c)> pred = v_is_space
    ) const {
        std::vector<self> splitted;
        size_t index = 0, front = 0;
        for (; index < storage.size();) {
            if (pred(storage[index])) {
                if (index > front) {
                    splitted.push_back(substr(front, index));
                }
                front = index + 1;
                index = front;
            } else ++index;
        }
        if (index > front)
            splitted.push_back(substr(front, index));
        return splitted;
    }

    std::vector<self> split(
        const self& sep,
        bool digest = true,
        MatchAlgo method = MatchAlgo::BruteForce
    ) const {
        std::vector<self> splitted;
        auto pos_list = findAll(sep, digest, method);
        size_t sep_size = sep.size();
        pos_list.push_back(npos);

        if (!pos_list.size()) {
            splitted.push_back(*this);
            return splitted;
        }
        if (pos_list[0] > 0)
            splitted.push_back(substr(0, pos_list[0]));
        size_t front = 0, back = 0;
        for (size_t i = 0; i < pos_list.size() - 1; ++i) {
            front = pos_list[i] + sep_size;
            back = pos_list[i + 1];
            if (front >= back) continue;
            splitted.push_back(substr(front, back));
        }
        return splitted;
    }

    // compatibility

    std::string toString() const {
        std::string res;
        for (auto chr : storage) chr.put(res);
        return res;
    }

    // operators

    bool operator==(const self& str) const {
        return compare(str) == CompareResult::equal;
    }
    bool operator!=(const self& str) const { return !operator==(str); }
    bool operator>(const self& str) const {
        return compare(str) == CompareResult::more;
    }
    bool operator<=(const self& str) const { return !operator>(str); }
    bool operator<(const self& str) const {
        return compare(str) == CompareResult::less;
    }
    bool operator>=(const self& str) const { return !operator<(str); }

    // iterators

    iterator begin() { return storage.begin(); }
    iterator end() { return storage.end(); }
    const_iterator cbegin() const { return storage.cbegin(); }
    const_iterator cend() const { return storage.cend(); }
};

typedef DynamicString<UTF8Char> U8String;

namespace literals {

    U8String operator ""_u8(const char* str, std::size_t) {
        return U8String(str);
    }

} // namespace dsl::literals

} // namespace dsl;

// hash specialization

namespace std {

    template<>
    struct hash<dsl::UTF8Char> {
        inline size_t operator() (dsl::UTF8Char chr) const noexcept {
            return hash<uint32_t>()(chr.codePoint());
        }
    };

}

#endif // cplusplus 20

#endif /* _DSL_STRINGPP_HPP_ */
