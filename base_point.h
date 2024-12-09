#pragma once

#include <iostream>
#include <concepts>
#include <numeric>
#include <algorithm>
#include <math.h>

namespace gmtr {

template<std::floating_point T>
//requires std::is_floating_point<T>::value
struct Point3D {
private:
    T v[3];
public:
    typedef T value_type;
    constexpr explicit Point3D( T x = T(0), T y = T(0), T z = T(0) )
        : v{x, y, z} {}

    constexpr T x() const noexcept {return v[0];}
    constexpr T y() const noexcept {return v[1];}
    constexpr T z() const noexcept {return v[2];}

    bool isNull() const noexcept {
        return std::all_of(std::begin(v), std::end(v), [](const auto& el){ return el == T(0); });
    }
    bool isNotValid() const noexcept {
        return std::any_of(std::begin(v), std::end(v), [](const auto& el){ return std::isnan(el); });
    }
    constexpr T normQuad() const noexcept {
        return std::accumulate(std::begin(v), std::end(v), T(0), [](T val, T x){ return val+=x*x;});
    }
    constexpr T norm() const noexcept {
        return sqrt(normQuad());
    }
    bool normalize() {
        auto nrm = normQuad();
        if(nrm > 0) {
            nrm = sqrt(nrm);
            *this /= nrm;
            return true;
        }
        return false;
    }
    constexpr Point3D normalized() const{
        auto nrm = normQuad();
        if(nrm > 0)
            return *this/norm();
        return *this;
    }

    constexpr auto operator+(const Point3D & point) const noexcept{
        return Point3D( x() + point.x(), y() + point.y(), z() + point.z() );
    }
    constexpr auto operator-(const Point3D & point) const noexcept{
        return Point3D( x() - point.x(), y() - point.y(), z() - point.z() );
    }
    constexpr auto operator*(const T & mult) const noexcept{
        return Point3D( x()*mult, y()*mult, z()*mult );
    }
    constexpr auto operator/(const T & mult) const {
        return Point3D( x()/mult, y()/mult, z()/mult );
    }
    constexpr auto operator*(const Point3D & point) const noexcept{
        return Point3D( x()*point.x(), y()*point.y(), z()*point.z() );
    }
    constexpr auto& operator+=(const Point3D & point) noexcept{
        v[0] += point.x();
        v[1] += point.y();
        v[2] += point.z();
        return *this;
    }
    constexpr auto& operator-=(const Point3D & point) noexcept{
        v[0] -= point.x();
        v[1] -= point.y();
        v[2] -= point.z();
        return *this;
    }
    constexpr auto& operator*=(const T & mult) noexcept{
        for(auto &i : v)
            i *= mult;
        return *this;
    }
    constexpr auto& operator/=(const T & div) {
        for(auto &i : v)
            i /= div;
        return *this;
    }
    constexpr auto operator-() const noexcept{
        return Point3D( -v[0], -v[1], -v[2] );
    }

    static constexpr T distanceQuad(const Point3D &p1, const Point3D &p2) noexcept{
        return (p2 - p1).normQuad();
    }

    inline static constexpr T distance(const Point3D &p1, const Point3D &p2) noexcept{
        return sqrt( distanceQuad(p1, p2) );
    }
    ///< выбрать удобную функцию dot
    static constexpr T dot(const Point3D &p1, const Point3D &p2) noexcept{
        return p1.x()*p2.x() + p1.y()*p2.y() + p1.z()*p2.z();
    }
    constexpr T dot(const Point3D &p2) const noexcept{
        return x()*p2.x() + y()*p2.y() + z()*p2.z();
    }

    static constexpr Point3D cross(const Point3D &p1, const Point3D &p2) noexcept{
        return Point3D( p1.y()*p2.z() - p1.z()*p2.y(),
                        p1.z()*p2.x() - p1.x()*p2.z(),
                        p1.x()*p2.y() - p1.y()*p2.x() );
    }
    ///< смешанное произведение векторов
    static constexpr T mix(const Point3D &p1, const Point3D &p2, const Point3D &p3) noexcept{
        return dot(cross(p1, p2), p3);
    }

    constexpr bool operator==(const Point3D &) const noexcept = default;
    constexpr bool operator!=(const Point3D &) const noexcept = default;

    template<typename U>
    constexpr Point3D<U> cast() const {return Point3D<U>(v[0], v[1], v[2]);}

    constexpr auto& shiftLeft() {
        std::rotate( std::begin(v),std::begin(v)+1, std::end(v) );
        return *this;
    }
    constexpr auto& shiftRight() {
        std::rotate( std::begin(v),std::begin(v)+2, std::end(v) );
        return *this;
    }
};

template<typename U, typename T>
constexpr auto operator*(const U& mult, const Point3D<T> & point) noexcept{
    return point*mult;
}

template<typename T>
std::ostream &operator<<(std::ostream &out, const Point3D<T> &point) {
    out << "(" << point.x() << ", " << point.y() << ", " << point.z() << ")";
    return out;
}

} // namespace gmtr
