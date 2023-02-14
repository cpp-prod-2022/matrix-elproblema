#include <cassert>
#include <complex>
#include <iostream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <algorithm>
#include <math.h>
#include <compare>
#include <array>

const int double_accuracy = 52;

class BigInteger {

    std::vector<int> digits;
        
public:
    enum Sign: int{
        Neg = -1,
        Zero,
        Pos
    };

    static constexpr int BASE = 10;

private:
    static constexpr double PI = M_PI;


    static constexpr double eps = 1e-5;

    Sign sign = Sign::Zero;

    void clear();

    void subtract(const BigInteger& n);
    //function calculating abs(this) - abs(n)
    void add(const BigInteger& n);
    //function calculating abs(this) + abs(n)
    void erase_lead_zeros();

    void divide(BigInteger n, bool div);
    //function divides abs(this) by abs(n)

    int reverse_number(int index, int bit_sz) {
        int res = 0;
        for (int i = 0; i < bit_sz; ++i) {
            if (1 << i & index) res += 1 << (bit_sz - 1 - i);
        }
        return res;
    }
    
    void fft(std::vector<std::complex<double>>&, bool);

    void multiply(const BigInteger & n);
    
    void increase();

    void decrease();

public:
    Sign get_sign() const {
        return sign;
    }

    bool abs_less(const BigInteger & n) const;
    //function returns abs(*this) < abs(n)

    BigInteger() = default;

    BigInteger(long long n);

    BigInteger(std::string s);

    BigInteger operator-();

    explicit operator int(); 

    explicit operator bool();

    BigInteger& operator += (BigInteger);
    
    BigInteger& operator -= (BigInteger);

    BigInteger& operator /= (BigInteger);

    BigInteger& operator *= (const BigInteger&);
    
    BigInteger& operator ++ ();
    
    BigInteger operator ++ (int);

    BigInteger operator -- (int);

    BigInteger& operator -- ();

    BigInteger& operator %= (BigInteger);
    
    std::string toString() const;

    void abs();

    void add_zeros(std::size_t);

// friend functions
};


bool operator <(const BigInteger& a, const BigInteger & b) {
    if (int(a.get_sign()) < int(b.get_sign())) return true;
    if (int(a.get_sign()) > int(b.get_sign())) return false;
    if (a.get_sign() == BigInteger::Sign::Zero) return false;
    if (a.get_sign() == BigInteger::Sign::Neg) return b.abs_less(a);
    return a.abs_less(b);
}

bool operator >(const BigInteger& a, const BigInteger& b) {
    return b < a;
}

bool operator ==(const BigInteger& a, const BigInteger& b) {
    return !(a < b) && !(b < a);
}

bool operator !=(const BigInteger& a, const BigInteger& b) {
    return !(a == b);
}

bool operator <=(const BigInteger& a, const BigInteger& b) {
    return !(a > b);
}

bool operator >=(const BigInteger& a, const BigInteger& b) {
    return !(a < b);
}

void BigInteger::add_zeros(std::size_t n) {
    for (size_t i = 0; i < n; ++i) digits.push_back(0);
    erase_lead_zeros();
}

BigInteger operator / (BigInteger, const BigInteger&);

BigInteger operator * (BigInteger, const BigInteger&);

BigInteger operator % (BigInteger, const BigInteger&);

void BigInteger::abs() {
    sign = Sign(::abs(int(sign)));
}

std::ostream& operator<<(std::ostream& out, const BigInteger& n) {
    out << n.toString();
    return out;
}

std::string BigInteger::toString() const {
    if (sign == Sign::Zero) return "0";
    std::string res;
    if (sign == Sign::Neg) res.push_back('-');
    for (size_t i = 0; i < digits.size(); ++i) res.push_back(char(digits[i] + int('0')));
    return res;
}

BigInteger::BigInteger(std::string s) {
    if (s[0] == '-') sign = Sign::Neg, s.erase(s.begin());
    else sign = Sign::Pos;
    for (size_t i = 0; i < s.size(); ++i) {
        digits.push_back(s[i] - '0');        
    }
    erase_lead_zeros();
}    

BigInteger operator "" _bi(unsigned long long n) {
    return BigInteger(static_cast<long long>(n));    
}

BigInteger& BigInteger::operator %= (BigInteger n) {
    divide(n, false); 
    return *this;
}

BigInteger operator% (BigInteger a, const BigInteger &b) {
    a %= b;
    return a;
}

BigInteger BigInteger::operator -- (int) {
    BigInteger tmp = *this;
    --*this;
    return tmp;
}

BigInteger& BigInteger::operator -- () {
    if (sign == Sign::Pos) decrease();
    else if (sign == Sign::Zero) increase(), sign = Sign::Neg;
    else increase();
    return *this;    
}    

BigInteger BigInteger::operator ++ (int) {
    BigInteger tmp = *this;
    ++*this;
    return tmp;    
}

void BigInteger::increase() {
    for (int i = int(digits.size() - 1); i > -1; --i) {
        if (digits[size_t(i)] == BASE - 1) {
            digits[size_t(i)] = 0;
        } else {
            ++digits[size_t(i)];
            return;
        }       
    }
    digits.insert(digits.begin(), 1);
}

void BigInteger::decrease() {
    if (sign == Sign::Zero) {
        digits.clear();
        digits = {1};
        sign = Sign::Neg;
        return;
    }
    for (int i = int(digits.size() - 1); i > -1; --i) {
        if (digits[size_t(i)] == 0) {
            digits[size_t(i)] = BASE - 1;
        } else {
            --digits[size_t(i)];
            break;
        }
    }
    erase_lead_zeros();
    if (digits.empty()) sign = Sign::Zero;
}

BigInteger& BigInteger::operator ++ () {
    if (sign == Sign::Neg) decrease();
    else if (sign == Sign::Zero) increase(), sign = Sign::Pos; 
    else increase();
    return *this;
}

BigInteger operator* (BigInteger a, const BigInteger &b) {
    a *= b;
    return a;
}

BigInteger& BigInteger::operator *= (const BigInteger &n) {
    if (n.sign == Sign::Zero) {
        clear();
        return *this;
    }
    sign = Sign(int(sign) * int(n.sign));
    multiply(n);
    return *this;
}

void BigInteger::fft(std::vector<std::complex<double>> &polynom, bool invert) {
    int sz = int(polynom.size());
    int lg = 0;
    while ((1 << lg) < sz) ++lg;
    for (int i = 0; i < sz; ++i) {
        if (i < reverse_number(i, lg)) std::swap(polynom[size_t(i)], polynom[size_t(reverse_number(i, lg))]);
    }
    for (size_t len = 2; len <= polynom.size(); len <<= 1) {
        double ang = 2 * PI / double(len);
        if (invert) ang = -ang;
        std::complex<double> w(cos(ang), sin(ang));
        for (size_t i = 0; i < polynom.size(); i += len) {
            std::complex<double> base(1);
            for (size_t j = 0; j < len / 2; ++j) {
                auto f = polynom[i + j], s = polynom[i + len / 2 + j] * base;
                polynom[i + j] = f + s;
                polynom[i + len / 2 + j] = f - s;
                base *= w;
            }
        }
    }
    if (invert) {
        for (int i = 0; i < sz; ++i) {
            polynom[size_t(i)] /= sz;
        }
    }
}

std::vector<std::complex<double>> fit_in_base(int cnt, int size, const std::vector<int> &a) {
    std::vector<std::complex<double>> res(size_t((size + cnt - 1) / cnt));
    for (size_t i = 0; i < a.size(); ++i) {
        res[(a.size() - 1U - i) / size_t(cnt)] *= BigInteger::BASE;
        res[(a.size() - 1U - i) / size_t(cnt)] += a[i];
    }
    return res;
}

void BigInteger::multiply(const BigInteger &n) {
    using complex = std::complex<double>;
    int sz_me = int(digits.size()), sz_n = int(n.digits.size());
    int log_res = 0;
    while ((1 << log_res) < sz_me + sz_n) log_res++;
    std::vector<complex> my_digits = fit_in_base(2, 1 << log_res, digits);
    std::vector<complex> n_digits = fit_in_base(2, 1 << log_res, n.digits);
    fft(my_digits, false);    
    fft(n_digits, false);
    for (size_t i = 0; i < my_digits.size(); ++i) my_digits[i] *= n_digits[i];
    fft(my_digits, true);
    digits.clear();
    int i = 0, carry = 0;
    while (i < int(2 * my_digits.size()) || carry) {
        if (int(digits.size()) <= i) digits.push_back(0);
        if (i % 2 == 0) carry += int(round(my_digits[size_t(i) / 2].real()));
        digits[size_t(i)] += carry;
        carry = digits[size_t(i)] / BASE;
        digits[size_t(i)] %= BASE;
        ++i;
    }
    std::reverse(digits.begin(), digits.end());
    erase_lead_zeros();
}

BigInteger& BigInteger::operator /= (BigInteger n) {
    sign = Sign(int(sign) * int(n.sign));    
    divide(n, true);
    return *this;
}

BigInteger operator / (BigInteger a, const BigInteger &b) {
    a /= b;
    return a;
}

void BigInteger::divide(BigInteger n, bool div) {
    int sz_n = int(n.digits.size());
    int s = sign;
    while (n.digits.size() < digits.size()) n.digits.push_back(0);
    std::vector<int> new_digits;
    while (int(n.digits.size()) >= sz_n) {
        new_digits.push_back(0);
        while (!abs_less(n)) subtract(n), ++new_digits.back();
        n.digits.pop_back();
    }
    if (div) digits = new_digits;
    sign = Sign(s);
    erase_lead_zeros();
}

void BigInteger::clear() {
    digits.clear();
    sign = Sign::Zero;
}

void BigInteger::erase_lead_zeros() {
    std::reverse(digits.begin(), digits.end());
    while (!digits.empty() && digits.back() == 0) digits.pop_back();
    if (digits.empty()) sign = Sign::Zero;
    std::reverse(digits.begin(), digits.end());
}

BigInteger operator +(BigInteger a, const BigInteger &b) {
    a += b;
    return a;
}

BigInteger operator -(BigInteger a, const BigInteger &b) {
    a -= b;
    return a;
}

BigInteger& BigInteger::operator-=(BigInteger n) {
    *this += -n;
    return *this;
}

bool BigInteger::abs_less(const BigInteger &n) const {
    if (sign == Sign::Zero) {
        if (n.sign == Sign::Zero) return false;
        return true;
    }
    if (digits.size() > n.digits.size()) return false;
       if (digits.size() < n.digits.size()) return true;
    for (size_t i = 0; i < digits.size(); ++i) {
        if (digits[i] < n.digits[i]) return true;
        if (digits[i] > n.digits[i]) return false;
    }    
    return false;
}

void BigInteger::subtract(const BigInteger &n) {
    if (n.sign == Sign::Zero) return;
    int i = int(digits.size()) - 1, j = int(n.digits.size()) - 1;
    while (i > -1) {
        if (j > -1) digits[size_t(i)] -= n.digits[size_t(j)];
        if (digits[size_t(i)] < 0) {
            digits[size_t(i)] += BASE;
            --digits[size_t(i - 1)];
        }
        --i, --j;
    }    
    erase_lead_zeros();
    if (digits.empty()) {sign = Sign::Zero;}
}

void BigInteger::add(const BigInteger &n) {
    if (n.sign == Sign::Zero) return;
    std::reverse(digits.begin(), digits.end());
    int i = 0, j = int(n.digits.size() - 1);
    int carry = 0;
    while (j > -1 || carry) {
        while (i >= int(digits.size())) digits.push_back(0);
        digits[size_t(i)] += carry;
        if (j > -1) digits[size_t(i)] += n.digits[size_t(j)];
        carry = digits[size_t(i)] / BASE;
        digits[size_t(i)] %= BASE;
        ++i, --j;
    }    
    std::reverse(digits.begin(), digits.end());
}

BigInteger& BigInteger::operator+=(BigInteger n) {
    if (sign == Sign::Zero) {
        *this = n;
        return *this;
    }
    if (n.sign == Sign::Zero) return *this;
    if (sign == n.sign) add(n);
    else {
        if (!abs_less(n)) subtract(n);
        else {
            n.subtract(*this);
            *this = n;
        }
    }
    return *this;
}


std::istream& operator >>(std::istream& in, BigInteger& n) {
    std::string s;
    in >> s;
    n = BigInteger(s);
    return in;
}

BigInteger::BigInteger(long long n) {
    if (n < 0) sign = Sign::Neg;
    else if (n > 0) sign = Sign::Pos;
    else sign = Sign::Zero;
    n = ::abs(n);
    while (n) {
        digits.push_back(int(n % BASE));
        n /= BASE;
    }
    std::reverse(digits.begin(), digits.end());
}

BigInteger BigInteger::operator-() {
    BigInteger n(*this);
    n.sign = Sign(int(n.sign) * -1);
    return n;
}

BigInteger::operator int() {
    if (sign == Sign::Zero) return 0;
    int res = 0;
    for (size_t i = 0; i < digits.size(); ++i) res <<= 1, res += digits[i];
    return res;
}

BigInteger::operator bool() {
    return sign != Sign::Zero;
}

////////////////////////
////////Rational////////
////////////////////////

class Rational {
    using Sign = BigInteger::Sign;
    
    Sign sign = Sign::Zero;

    BigInteger num;

    BigInteger denom;

    void simplify();

    void reverse();

    BigInteger gcd(BigInteger x, BigInteger y);

    public:

    Rational(): num(0), denom(1) {}

    Rational(const BigInteger&);

    Rational(long long);

    Rational& operator += (const Rational&);

    Rational& operator -= (Rational);
    
    Rational& operator *= (const Rational&);

    Rational& operator /= (Rational);

    Rational operator - ();
    
    std::string toString() const;

    std::string asDecimal(std::size_t);

    explicit operator double ();

    //friend functions
    
    friend bool operator < (const Rational &a, const Rational &b);
};

bool operator < (const Rational &a, const Rational &b) {
    if (int(a.sign) < int(b.sign)) return true;
    if (int(a.sign) > int(b.sign)) return false;
    if (a.sign == BigInteger::Sign::Zero) return false;
    if (a.sign == BigInteger::Sign::Neg) {
        return a.num * b.denom > b.num * a.denom;
    }
    return a.num * b.denom < b.num * a.denom;
}

bool operator > (const Rational &a, const Rational &b) {
    return b < a;
}

bool operator == (const Rational &a, const Rational &b) {
    return !(a < b) && !(b < a);
}

bool operator != (const Rational &a, const Rational &b) {
    return !(a == b);
}

bool operator <= (const Rational &a, const Rational &b) {
    return !(a > b);
}

bool operator >= (const Rational &a, const Rational &b) {
    return !(a < b);
}

Rational operator + (Rational, const Rational&);

Rational operator - (Rational, const Rational&);

Rational operator * (Rational, const Rational&);

Rational operator / (Rational, const Rational&);

std::ostream& operator <<(std::ostream& out, const Rational &n) {
    out << n.toString();
    return out;
}

void Rational::reverse() {
    std::swap(num, denom);
}

Rational Rational::operator - () {
    Rational n(*this);
    n.sign = Sign(int(n.sign) * -1);
    return n;
}

Rational& Rational::operator += (const Rational &n) {
    num = int(sign) * num * n.denom + int(n.sign) * n.num * denom;
    denom *= n.denom;            
    simplify();
    return *this;
}

Rational& Rational::operator -= (Rational n) {
    *this += (-n);
    return *this;
}

Rational& Rational::operator *= (const Rational &n) {
    num = int(sign) * int(n.sign) * num * n.num;
    denom *= n.denom;
    simplify();
    return *this;
}

Rational& Rational::operator /= (Rational n) {
    n.reverse();
    *this *= n;
    return *this;
}

Rational operator + (Rational a, const Rational &b) {
    a += b;
    return a;
}

Rational operator - (Rational a, const Rational &b) {
    a -= b;
    return a;
}

Rational operator / (Rational a, const Rational &b) {
    a /= b;
    return a;
}

Rational operator * (Rational a, const Rational &b) {
    a *= b;
    return a;
}

void Rational::simplify() {
    if (num.get_sign() == BigInteger::Sign::Zero) sign = Sign::Zero;
    else {
        int s = 1;
        s *= int(num.get_sign());
        s *= int(denom.get_sign());
        sign = Sign(s);
    }
    num.abs(), denom.abs();
    BigInteger res = gcd(num, denom);
    num /= res, denom /= res;        
}

BigInteger Rational::gcd(BigInteger x, BigInteger y) {
    if (y < x) std::swap(x, y);
    while (x) {
        y %= x;
        std::swap(x, y);
    }
    return y;    
}

Rational::Rational(const BigInteger &n): num(n), denom(1) {
    simplify();
}

Rational::Rational(long long n): num(n), denom(1) {
    simplify();
}

std::string Rational::toString() const {
    std::string res;
    if (sign == Sign::Neg) res = "-";
    res += num.toString();
    if (denom != BigInteger(1)) res += "/", res += denom.toString();
    return res;
}

std::string Rational::asDecimal(std::size_t size) {
    BigInteger cpy = num;
    num.add_zeros(size);
    std::string s = (num / denom).toString();
    if (s.size() > size) {
        s = s.substr(0, s.size() - size) + "." + s.substr(s.size() - size, size);
    } else {
        s = "0." + std::string(size - s.size(), '0') + s;    
    }
    num = cpy;
    if (sign == Sign::Neg) s = "-" + s;
    return s;    
}

Rational::operator double() {
    std::string s = asDecimal(double_accuracy);
    return std::stod(s);    
}

std::istream& operator >> (std::istream& in, Rational &r) {
    std::string s; 
    std::cout << s;
    in >> s;
    r = Rational(BigInteger(s));
    return in;
}

template<size_t D, size_t N>
struct Iterator {
    static const bool value = (D * D > N? true: Iterator<D + 1, N>::value && N % D != 0);
};

template<size_t N>
struct Iterator<N, N> {
    static const bool value = true;
};

template<size_t D>
struct Iterator<D, 1> {
    static const bool value = false;
};

template<size_t N>
struct is_prime {
    static const bool value = Iterator<2, N>::value;
};

template <size_t N>
class Residue {
    size_t x = 0;
    
    static int fast_pow(size_t x, size_t y) {
        unsigned long long cnt = static_cast<unsigned long long>(x), ans = 1;
        for (int i = 0; (1 << i) <= y; ++i) {
            if ((1 << i) & y) ans *= cnt, ans %= N;
            cnt *= cnt; cnt %= N;
        }
        return ans % N;
    }

public:

    explicit operator int() {
        return x;
    }
    
    Residue& operator /= (const Residue<N>& r) {
        static_assert(is_prime<N>::value);
        unsigned long long t = fast_pow(r.x, N - 2);
        x = (static_cast<unsigned long long>(x) * t) % N;
        return *this;
    }

    Residue& operator += (const Residue<N>& r) {
        x += r.x;
        x %= N;
        return *this;
    }    
    
    Residue& operator -= (const Residue<N>& r) {
        x += N - r.x;
        x %= N;
        return *this;
    }

    Residue& operator *= (const Residue<N>& r) {
        x = (static_cast<unsigned long long>(r.x) * static_cast<unsigned long long>(x)) % N;
        return *this;
    }

    bool operator == (const Residue<N>&) const = default;
    bool operator != (const Residue<N>&) const = default;

    explicit Residue(long long n) {
        long long NLL = static_cast<long long>(N);
        x = int((n % NLL + NLL) % NLL);
    }
    
    Residue() = default;
};

template<size_t N>
Residue<N> operator + (Residue<N> x, const Residue<N>& y) {
    x += y;
    return x;
}

template<size_t N>
Residue<N> operator - (Residue<N> x, const Residue<N>& y) {
    x -= y;
    return x;
}

template<size_t N>
Residue<N> operator * (Residue<N> x, const Residue<N>& y) {
    x *= y;
    return x;
}

template<size_t N>
Residue<N> operator / (Residue<N> x, const Residue<N>& y) {
    x /= y;
    return x;
}

template<size_t N>
std::istream& operator >> (std::istream& in, Residue<N>& x) {
    long long t = 0; 
    in >> t;
    x = Residue<N>(t);
    return in;
}

template<size_t N>
std::ostream& operator << (std::ostream& out, Residue<N> x) {
    out << int(x);
    return out;
}



template<size_t M, size_t N, typename Field=Rational>
class Matrix {
public:
    
    size_t find_leader(const std::array<Field, N>& x) const {
        for (size_t j = 0; j < x.size(); ++j) {
            if (x[j] != 0) return j;
        }
        return x.size();
    }

private:

    std::array<std::array<Field, N>, M> m_data;
    using vec = std::vector<std::vector<Field>>;    
    using arr = std::array<std::array<Field, N>, M>;
    
    template<size_t K>
    std::array<std::array<Field, K>, M> trivial_multiplication(const Matrix<N, K, Field>& to_mul) const {
        std::array<std::array<Field, K>, M> res;
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                for (size_t k = 0; k < K; ++k) {
                    res[i][k] += m_data[i][j] * to_mul[j][k];
                }
            }
        }
        return res;
    }
    
    Field gauss() {
        Field res(1);
        bool sign = false;
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = M - 1; j > i; --j) {
                if (find_leader(m_data[j - 1]) > find_leader(m_data[j])) {
                    swap(m_data[j - 1], m_data[j]);
                    sign ^= 1;
                }
            }
            size_t j = find_leader(m_data[i]);
            if (j == N) break;
            for (size_t k = i + 1; k < M; ++k) {
                Field c = m_data[k][j] / m_data[i][j];
                for (size_t l = 0; l < N; ++l) {
                    m_data[k][l] -= m_data[i][l] * c;
                }
            }
        }
        if (sign) res *= Field(-1);
        for (size_t i = 0; i < M; ++i) {
            size_t j = find_leader(m_data[i]);
            if (j == N) {
                res = 0;
                break;
            }
            Field x = m_data[i][j];
            res *= x;
            for (size_t k = 0; k < N; ++k) m_data[i][k] /= x;
        }
        for (size_t i = M - 1; i < M; --i) {
            size_t j = find_leader(m_data[i]);
            if (j == N) continue;
            for (size_t k = 0; k < i; ++k) {
                Field c = m_data[k][j] / m_data[i][j];    
                for (size_t l = 0; l < N; ++l) {
                    m_data[k][l] -= m_data[i][l] * c;
                }
            }
        }
        return res;
    }

public:

    Matrix() = default;

    Matrix(const vec& mat)  {
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                m_data[i][j] = mat[i][j];
            }
        }
    }

    Matrix(const arr& mat): m_data(mat) {}

    Matrix(std::initializer_list<std::vector<Field>> mat): 
            Matrix(vec(mat)) {}

    Field trace() const {
        static_assert(M == N);
        Field res = 0;
        for (size_t i = 0; i < M; ++i) {
            res += m_data[i][i];
        }
        return res;
    }

    std::array<Field, M> getColumn(size_t j) const {
        std::array<Field, M> res;
        for (size_t i = 0; i < M; ++i) res[i] = m_data[i][j];
        return res;
    }    

    std::array<Field, N> getRow(size_t i) const {
        std::array<Field, N> res;
        for (size_t j = 0; j < N; ++j) res[j] = m_data[i][j];
        return res;
    }
    
    std::array<Field, N>& operator[](size_t i) {
        return m_data[i];
    }

    const std::array<Field, N>& operator[](size_t i) const {
        return m_data[i];
    }

    Matrix<N, M, Field> transposed() const {
        Matrix<N, M, Field> trans;
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                trans[j][i] = m_data[i][j];
            }
        }
        return Matrix<N, M, Field>(trans);
    }

    Matrix<M, N, Field>& operator *= (const Field& lambda) {
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                m_data[i][j] *= lambda;
            }
        }
        return *this;
    }

    Matrix<M, N, Field>& operator *= (const Matrix<N, N, Field>& x) {
        trivial_multiplication(x);
        return *this;
    }
    
    Matrix<M, N, Field>& operator += (const Matrix<M, N, Field>& b) {
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                m_data[i][j] += b[i][j];
            }
        }
        return *this;
    }    


    Matrix<M, N, Field>& operator -= (const Matrix<M, N, Field>& b) {
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                m_data[i][j] -= b[i][j];
            }
        }
        return *this;
    }

    size_t rank() const {
        Matrix<M, N, Field> cpy = *this;
        cpy.gauss();
        for (size_t i = 0; i < M; ++i) {
            if (find_leader(cpy[i]) == N) return i;
        }
        return M;
    }

    Field det() const {
        static_assert(N == M);
        Matrix<M, N, Field> cpy = *this;
        Field res = cpy.gauss();
        return res;
    }

    void invert() {
        static_assert(M == N);
        Matrix<M, 2 * N, Field> cpy;
        for (size_t i = 0; i < M; ++i) {
            m_data[i][N + i] = 1;
            for (size_t j = 0; j < N; ++j) {
                cpy[i][j] = m_data[i][j];
            }
        }
        cpy.gauss();
        assert(find_leader(cpy[M - 1]) == M - 1);
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                m_data[i][j] = cpy[i][N + j];
            }
        }
    }

    auto inverted() {
        Matrix<M, N, Field> cpy = *this;
        cpy.invert();
        return cpy;
    }
};

template<size_t M, size_t N, typename Field>
bool operator == (const Matrix<M, N, Field>& a, const Matrix<M, N, Field>& b) {
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            if (a[i][j] != b[i][j]) return false;
        }    
    }
    return true;
}

template<size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator - (const Matrix<M, N, Field>& a, const Matrix<M, N, Field>& b) {
    Matrix<M, N, Field> copy = a;
    copy -= b;
    return copy;
}

template<size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator + (const Matrix<M, N, Field>& a, const Matrix<M, N, Field>& b) {
    Matrix<M, N, Field> copy = a;
    copy += b;
    return copy;
}

template<size_t M_, size_t N_, size_t K_, typename F>
Matrix<M_, K_, F> operator * (const Matrix<M_, N_, F>& A, const Matrix<N_, K_, F>& B) {
    return Matrix<M_, K_, F>(A.trivial_multiplication(B));
}

template<size_t M, size_t N, typename F>
Matrix<M, N, F> operator * (F x, Matrix<M, N, F> m) {
    m *= x;
    return m;
}

template<size_t N, typename Field = Rational>
using SquareMatrix = Matrix<N, N, Field>;

