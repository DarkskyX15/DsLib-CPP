/*
 * @Date: 2024-09-21 13:30:58
 * @Author: DarkskyX15
 * @LastEditTime: 2024-10-28 16:39:12
 */

#ifndef _NUMPP_HPP_
#define _NUMPP_HPP_

#include <iostream>
#include <stdint.h>
#include <string>
#include <vector>

namespace dsl {

    struct DigitUnit {
        uint8_t data;
        DigitUnit() { data = 0; }
        DigitUnit(uint8_t h, uint8_t l) { hset(h); lset(l); }
        uint8_t hget() const { return (data >> 4); }
        uint8_t lget() const { return (data & 0b1111); }
        void hset(uint8_t val) {
            data &= 0b1111;
            data |= (val << 4);
        }
        void lset(uint8_t val) {
            data &= 0b11110000;
            data |= val;
        }
    };

    class BigNum {

    friend std::istream& operator>> (std::istream& in, BigNum& bn) {
        std::string raw;
        in >> raw;
        bn.operator=(BigNum(raw));
        return in;
    }
    friend std::ostream& operator<< (std::ostream& out, const BigNum& bn) {
        if (!bn.sign) out << '-';
        auto iter = bn.storage.rbegin();
        auto end = bn.storage.rend();
        if (bn.high) { out << static_cast<uint16_t>(iter->hget()); }
        out << static_cast<uint16_t>(iter->lget());
        ++iter;
        while (iter != end) {
            out << static_cast<uint16_t>(iter->hget())
                << static_cast<uint16_t>(iter->lget());
            ++iter;
        }
        return out;
    }

    private:
        std::vector<DigitUnit> storage;
        bool sign, is_zero, high;

        void set_zero() {
            storage.clear();
            sign = true;
            is_zero = true;
            high = false;
            storage.push_back(DigitUnit(0, 0));
        }

        enum class CompareResult: uint8_t { greater, equal, lesser };
        CompareResult signed_compare(const BigNum& num) const {
            // sign judge
            if (sign && !num.sign) return CompareResult::greater;
            if (!sign && num.sign) return CompareResult::lesser;
            bool positive = sign;
            CompareResult abs_result = abs_compare(num);
            if (abs_result != CompareResult::equal) {
                if (positive) return abs_result;
                else {
                    if (abs_result == CompareResult::greater)
                        return CompareResult::lesser;
                    if (abs_result == CompareResult::lesser)
                        return CompareResult::greater;
                }
            }
            return abs_result;
        }
        CompareResult abs_compare(const BigNum& num) const {
            CompareResult abs_result = CompareResult::equal;
            // size judge
            if (storage.size() < num.storage.size())
                abs_result = CompareResult::lesser;
            else if (storage.size() > num.storage.size())
                abs_result = CompareResult::greater;
            else {
                // same length
                if (!high && num.high) abs_result = CompareResult::lesser;
                else if (high && !num.high) abs_result = CompareResult::greater;
            }
            if (abs_result != CompareResult::equal) return abs_result;
            // same size
            auto iter1 = storage.rbegin(), iter2 = num.storage.rbegin();
            auto end1 = storage.rend(), end2 = num.storage.rend();
            uint8_t digit1, digit2;
            while (iter1 != end1) {
                // high
                digit1 = iter1->hget();
                digit2 = iter2->hget();
                if (digit1 < digit2) {
                    abs_result = CompareResult::lesser;
                    break;
                } else if (digit1 > digit2) {
                    abs_result = CompareResult::greater;
                    break;
                }
                // low
                digit1 = iter1->lget();
                digit2 = iter2->lget();
                if (digit1 < digit2) {
                    abs_result = CompareResult::lesser;
                    break;
                } else if (digit1 > digit2) {
                    abs_result = CompareResult::greater;
                    break;
                }
                ++iter1; ++iter2;
            }
            return abs_result;
        }

        void unsigned_minus(const BigNum& num, CompareResult cmp) {
            if (cmp == CompareResult::greater) {
                // bigger
                auto iter1 = storage.begin(), end1 = storage.end();
                auto iter2 = num.storage.begin(), end2 = num.storage.end();
                uint8_t digit1, digit2;
                bool borrow = false;
                // both forward
                while (iter2 != end2) {
                    // low
                    digit1 = iter1->lget();
                    digit2 = iter2->lget();
                    digit2 += borrow;
                    if (digit1 < digit2) {
                        digit1 += 10 - digit2;
                        borrow = true;
                    } else {
                        digit1 -= digit2;
                        borrow = false;
                    }
                    iter1->lset(digit1);
                    // high
                    digit1 = iter1->hget();
                    digit2 = iter2->hget();
                    digit2 += borrow;
                    if (digit1 < digit2) {
                        digit1 += 10 - digit2;
                        borrow = true;
                    } else {
                        digit1 -= digit2;
                        borrow = false;
                    }
                    iter1->hset(digit1);
                    ++iter1; ++iter2;
                }
                // top forward
                bool half_front;
                while (iter1 != end1) {
                    // half front
                    half_front = false;
                    if (iter1->data <= 15 && iter1->data != 0) half_front = true;
                    // low
                    digit1 = iter1->lget();
                    if (digit1 < borrow) {
                        digit1 += 10 - borrow;
                        borrow = true;
                    } else {
                        digit1 -= borrow;
                        borrow = false;
                    };
                    iter1->lset(digit1);
                    if (half_front) break;
                    // high
                    digit1 = iter1->hget();
                    if (digit1 < borrow) {
                        digit1 += 10 - borrow;
                        borrow = true;
                    } else {
                        digit1 -= borrow;
                        borrow = false;
                    };
                    iter1->hset(digit1);
                    ++iter1;
                }
                // backward
                while (true) {
                    digit1 = storage.back().data;
                    if (digit1 == 0) {
                        storage.pop_back();
                        continue;
                    } else if (digit1 <= 15) {
                        high = false;
                        break;
                    } else {
                        high = true;
                        break;
                    }
                }
            } else {
                // smaller
                sign = false;
                // padding
                size_t offset = num.storage.size() - storage.size();
                while(offset--) storage.push_back(DigitUnit(0, 0));
                // both forward
                auto iter1 = storage.begin(), end1 = storage.end();
                auto iter2 = num.storage.begin(), end2 = num.storage.end();
                uint8_t digit1, digit2;
                bool borrow = false;
                // both forward
                while (iter2 != end2) {
                    // low
                    digit1 = iter1->lget();
                    digit2 = iter2->lget();
                    digit1 += borrow;
                    if (digit2 < digit1) {
                        digit2 += 10 - digit1;
                        borrow = true;
                    } else {
                        digit2 -= digit1;
                        borrow = false;
                    }
                    iter1->lset(digit2);
                    // high
                    digit1 = iter1->hget();
                    digit2 = iter2->hget();
                    digit1 += borrow;
                    if (digit2 < digit1) {
                        digit2 += 10 - digit1;
                        borrow = true;
                    } else {
                        digit2 -= digit1;
                        borrow = false;
                    }
                    iter1->hset(digit2);
                    ++iter1; ++iter2;
                }
                // backward
                while (true) {
                    digit1 = storage.back().data;
                    if (digit1 == 0) {
                        storage.pop_back();
                        continue;
                    } else if (digit1 <= 15) {
                        high = false;
                        break;
                    } else {
                        high = true;
                        break;
                    }
                }
            }
        }

        void unsigned_add(const BigNum& num) {
            // padding
            if (storage.size() < num.storage.size()) {
                size_t offset = num.storage.size() - storage.size();
                while (offset--) storage.push_back(DigitUnit(0, 0));
            }
            // both forward
            auto iter1 = storage.begin(), end1 = storage.end();
            auto iter2 = num.storage.begin(), end2 = num.storage.end();
            uint8_t carry = 0;
            while (iter2 != end2) {
                // low
                carry += iter1->lget() + iter2->lget();
                iter1->lset(carry % 10);
                carry /= 10;
                // high
                carry += iter1->hget() + iter2->hget();
                iter1->hset(carry % 10);
                carry /= 10;
                ++iter1; ++iter2;
            }
            // top forward
            while (iter1 != end1) {
                // low
                carry += iter1->lget();
                iter1->lset(carry % 10);
                carry /= 10;
                // high
                carry += iter1->hget();
                iter1->hset(carry % 10);
                carry /= 10;
                ++iter1;
            }
            // last forward
            if (carry > 0) storage.push_back(DigitUnit(0, carry));
            // backward
            auto riter = storage.rbegin();
            if (riter->hget()) high = true;
            else high = false;
        }

        void unsigned_multiply(const BigNum& num) {
            // padding
            BigNum result;
            result.storage.resize(storage.size() + num.storage.size());
            // iteration
            auto iter1 = storage.begin(), end1 = storage.end();
            auto iter2 = num.storage.begin(), end2 = num.storage.end();
            auto writer = result.storage.begin(), wt = writer;
            uint8_t digit1, digit2, carry = 0;
            for (;iter1 != end1; ++iter1, ++writer) {
                // low
                digit1 = iter1->lget();
                wt = writer;
                carry = 0;
                for (; iter2 != end2; ++iter2, ++wt) {
                    // low
                    digit2 = iter2->lget();
                    carry += digit2 * digit1;
                    carry += wt->lget();
                    wt->lset(carry % 10);
                    carry /= 10;
                    // high
                    digit2 = iter2->hget();
                    carry += digit2 * digit1;
                    carry += wt->hget();
                    wt->hset(carry % 10);
                    carry /= 10;
                }
                wt->lset(carry);
                // reset iterator
                iter2 = num.storage.begin();
                // high
                digit1 = iter1->hget();
                wt = writer;
                carry = 0;
                for (; iter2 != end2; ++iter2) {
                    // low
                    digit2 = iter2->lget();
                    carry += digit2 * digit1;
                    carry += wt->hget();
                    wt->hset(carry % 10);
                    carry /= 10;
                    // high
                    ++wt;
                    digit2 = iter2->hget();
                    carry += digit2 * digit1;
                    carry += wt->lget();
                    wt->lset(carry % 10);
                    carry /= 10;
                }
                wt->hset(carry);
                // reset iterator
                iter2 = num.storage.begin();
            }
            // backward
            while (true) {
                digit1 = result.storage.back().data;
                if (digit1 == 0) {
                    result.storage.pop_back();
                    continue;
                } else if (digit1 <= 15) {
                    result.high = false;
                    break;
                } else {
                    result.high = true;
                    break;
                }
            }
            result.is_zero = false;
            operator=(std::move(result));
        }

    public:
        /// @brief 空构造
        BigNum(): storage(1), high(false), sign(true), is_zero(true) {}
        ~BigNum() = default;
        BigNum(const BigNum& num) {
            this->storage = num.storage;
            this->high = num.high;
            this->is_zero = num.is_zero;
            this->sign = num.sign;
        };
        BigNum& operator= (const BigNum& num) {
            this->storage = num.storage;
            this->high = num.high;
            this->is_zero = num.is_zero;
            this->sign = num.sign;
            return *this;
        };
        BigNum(BigNum&& num) {
            this->storage = std::move(num.storage);
            this->high = num.high;
            this->is_zero = num.is_zero;
            this->sign = num.sign;
        }
        BigNum& operator= (BigNum&& num) {
            this->storage = std::move(num.storage);
            this->high = num.high;
            this->is_zero = num.is_zero;
            this->sign = num.sign;
            return *this;
        }

        /// @brief 由数字构造
        /// @param num 数字
        explicit BigNum(int64_t num) {
            // set sign
            if (num < 0) {
                sign = false;
                num *= -1;
            } else sign = true;

            // if zero
            if (num == 0) {
                is_zero = true;
                high = false;
                storage.push_back(DigitUnit(0, 0));
                return ;
            }

            is_zero = false;
            high = true;
            while (num) {
                uint8_t digit1 = num % 10;
                num /= 10;
                if (num) {
                    uint8_t digit2 = num % 10;
                    num /= 10;
                    storage.push_back(DigitUnit(digit2, digit1));
                } else {
                    high = false;
                    storage.push_back(DigitUnit(0, digit1));
                }
            }
        }

        /// @brief 由字符串构造
        /// @param str 字符串
        explicit BigNum(const std::string& str) {
            // empty str
            if (!str.size()) {
                storage.push_back(DigitUnit(0, 0));
                high = false;
                is_zero = true;
                sign = true;
                return;
            }
            // sign
            size_t offset = 0;
            if (str[0] == '-') {
                offset = 1;
                sign = false;
                // empty digit
                if (str.size() == 1) {
                    storage.push_back(DigitUnit(0, 0));
                    sign = true;
                    is_zero = true;
                    high = false;
                    return ;
                }
            } else sign = true;

            // erase front zero
            while (offset < str.size() && str[offset] == '0') ++offset;

            // all zero
            if (offset >= str.size()) {
                storage.push_back(DigitUnit(0, 0));
                sign = true;
                is_zero = true;
                high = false;
                return ;
            }

            // get digits
            size_t index = str.size() - 1;
            is_zero = false;
            high = true;
            while(index >= offset && index < str.size()) {
                uint8_t digit1 = str[index] - '0';
                --index;

                if (index < offset || index >= str.size()) {
                    high = false;
                    storage.push_back(DigitUnit(0, digit1));
                    break;
                } else {
                    uint8_t digit2 = str[index] - '0';
                    --index;
                    storage.push_back(DigitUnit(digit2, digit1));
                }
            }
        }

        /// @brief 返回字符串表示
        /// @return 字符串
        std::string str() const {
            std::string value;
            if (!sign) value.push_back('-');
            if (is_zero) {
                value.push_back('0');
                return value;
            }
            auto iter = storage.rbegin(), end = storage.rend();
            if (high) {
                value.push_back(iter->hget() + '0');
                value.push_back(iter->lget() + '0');
            } else {
                value.push_back(iter->lget() + '0');
            }
            ++iter;
            while (iter != end) {
                value.push_back(iter->hget() + '0');
                value.push_back(iter->lget() + '0');
                ++iter;
            }
            return value;
        }
        /// @brief 反转数字符号
        /// @return 自身引用
        BigNum& flip() {
            sign = sign ? false : true;
            return *this;
        }

        operator bool() const { return !is_zero; }

        bool operator== (const BigNum& num) const { return !operator!=(num); }
        bool operator!= (const BigNum& num) const {
            if (num.high != this->high) return true;
            if (num.is_zero != this->is_zero) return true;
            if (num.sign != this->sign) return true;
            if (num.storage.size() != this->storage.size()) return true;

            auto iter1 = storage.rbegin(), iter2 = num.storage.rbegin();
            auto end1 = storage.rend(), end2 = num.storage.rend();
            while (iter1 != end1 && iter2 != end2) {
                if (iter1->data != iter2->data) return true;
                ++iter1;
                ++iter2;
            }
            return false;
        }

        bool operator< (const BigNum& num) const {
            auto result = signed_compare(num);
            if (result == CompareResult::lesser) return true;
            return false;
        }
        bool operator> (const BigNum& num) const {
            auto result = signed_compare(num);
            if (result == CompareResult::greater) return true;
            return false;
        }
        bool operator<= (const BigNum& num) const { return !operator>(num); }
        bool operator>= (const BigNum& num) const { return !operator<(num); }

        BigNum operator- () const {
            BigNum value(*this);
            value.sign = value.sign ? false : true;
            return value;
        }
        
        BigNum& operator+= (const BigNum& num) {
            // zero check
            if (num.is_zero) return *this;
            if (this->is_zero) return operator=(num);
            // sign check
            // a + b
            if (sign && num.sign) unsigned_add(num);
            // -a + -b
            else if (!sign && !num.sign) unsigned_add(num);
            else {
                // |a| == |b|
                CompareResult cr = abs_compare(num);
                if (cr == CompareResult::equal) { set_zero(); return *this; }
                // a + -b
                if (sign && !num.sign) unsigned_minus(num, cr);
                // -a + b = -(a - b)
                else {
                    unsigned_minus(num, cr);
                    sign = sign ? false : true;
                }
            }
            return *this;
        }
        BigNum operator+ (const BigNum& num) const {
            BigNum value(*this);
            return value.operator+=(num);
        }

        BigNum& operator-= (const BigNum& num) {
            // zero check
            if (num.is_zero) return *this;
            if (is_zero) { operator=(num); return flip(); }
            // sign check
            // a - -b = a + b
            if (sign && !num.sign) unsigned_add(num);
            // -a - b = -(a + b)
            else if (!sign && num.sign) unsigned_add(num);
            else {
                CompareResult cr = abs_compare(num);
                if (cr == CompareResult::equal) { set_zero(); return *this; }
                // a - b
                if (sign && num.sign) unsigned_minus(num, cr);
                // -a - -b = - (a - b)
                else {
                    unsigned_minus(num, cr);
                    sign = sign ? false : true;
                }
            }
            return *this;
        }
        BigNum operator- (const BigNum& num) const {
            BigNum value(*this);
            return value.operator-=(num);
        }

        BigNum& operator*= (const BigNum& num) {
            // zero check
            if (is_zero) return *this;
            if (num.is_zero) { set_zero(); return *this; }
            // sign check
            if ((sign && num.sign) || (!sign && !num.sign)) {
                unsigned_multiply(num);
            } else {
                unsigned_multiply(num);
                sign = false;
            }
            return *this;
        }
        BigNum operator* (const BigNum& num) const {
            BigNum value(*this);
            return value.operator*=(num);
        }

    };

namespace literals {
    BigNum operator ""_bn(const char* str, std::size_t) {
        return BigNum(std::string(str));
    }
} // namespace dsl::literals

} // namespace dsl

#endif /* _NUMPP_HPP_ */
