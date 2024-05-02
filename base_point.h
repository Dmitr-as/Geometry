#ifndef BASE_POINT_H
#define BASE_POINT_H

#include <iostream>
#include <concepts>
#include <numeric>
#include <math.h>

template<typename T>
requires std::is_floating_point<T>::value
struct Point3D {
private:
    T v[3];
public:
    constexpr explicit Point3D( T x = T(0), T y = T(0), T z = T(0) )
        : v{x, y, z} {}

    const T &x() const {return v[0];}
    const T &y() const {return v[1];}
    const T &z() const {return v[2];}

    bool isNull() const {
        return std::all_of(std::begin(v), std::end(v), [](const auto& el){ return el == 0; });
        //return v[0] == 0. && v[1] == 0. && v[2] == 0.;
    }
    constexpr T normQuad() const {
        return std::accumulate(std::begin(v), std::end(v), 0, [](T val, T x){ return val+=x*x;});
        //return v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
    }
    constexpr T norm() const {
        return sqrt(normQuad());
    }
    bool normalaze() {
        auto nrm = normQuad();
        if(nrm > 0) {
            nrm = sqrt(nrm);
            v[0] /= nrm;
            v[1] /= nrm;
            v[2] /= nrm;
            return true;
        }
        return false;
    }
    constexpr Point3D normalazed() const {
        return Point3D(*this)/norm();
    }

    auto operator+(const Point3D & point) const {
        return Point3D( x() + point.x(), y() + point.y(), z() + point.z() );
    }
    auto operator-(const Point3D & point) const {
        return Point3D( x() - point.x(), y() - point.y(), z() - point.z() );
    }
    auto operator*(const T & mult) const {
        return Point3D( x()*mult, y()*mult, z()*mult );
    }
    auto operator/(const T & mult) const {
        return Point3D( x()/mult, y()/mult, z()/mult );
    }
    auto operator*(const Point3D & point) const {
        return Point3D( x()*point.x(), y()*point.y(), z()*point.z() );
    }
    auto& operator+=(const Point3D & point) {
        v[0] += point.x();
        v[1] += point.y();
        v[2] += point.z();
        return *this;
    }
    auto& operator-=(const Point3D & point) {
        v[0] -= point.x();
        v[1] -= point.y();
        v[2] -= point.z();
        return *this;
    }
    auto& operator*=(const T & mult) {
        v[0] *= mult;
        v[1] *= mult;
        v[2] *= mult;
        return *this;
    }
    auto& operator/=(const T & div) {
        v[0] /= div;
        v[1] /= div;
        v[2] /= div;
        return *this;
    }
    auto operator-() const {
        return Point3D( -v[0], -v[1], -v[2] );
    }

    static T distanceQuad(const Point3D &p1, const Point3D &p2) {
        return (p2 - p1).normQuad();
    }

    static T distance(const Point3D &p1, const Point3D &p2) {
        return sqrt( distanceQuad(p1, p2) );
    }

    static T dot(const Point3D &p1, const Point3D &p2) {
        return p1.x()*p2.x() + p1.y()*p2.y() + p1.z()*p2.z();
    }
    static Point3D cross(const Point3D &p1, const Point3D &p2) {
        return Point3D( p1.y()*p2.z() - p1.z()*p2.y(),
                        p1.z()*p2.x() - p1.x()*p2.z(),
                        p1.x()*p2.y() - p1.y()*p2.x() );
    }
    static T mix(const Point3D &p1, const Point3D &p2, const Point3D &p3) { // смешанное произведение векторов
        return dot(cross(p1, p2), p3);
    }

    bool operator==(const Point3D &) const = default;
    bool operator!=(const Point3D &) const = default;

    template<typename U>
    Point3D<U> cast() const {return Point3D<U>(v[0], v[1], v[2]);}

    friend std::ostream &operator<<(std::ostream &out, const Point3D &point) {
        out << "(" << point.x() << ", " << point.y() << ", " << point.z() << ")";
        return out;
    }

    void shiftLeft() {
        std::rotate( std::begin(v),std::begin(v)+1, std::end(v) );
    }
    void shiftRight() {
        std::rotate( std::begin(v),std::begin(v)+2, std::end(v) );
    }
};

template<typename T>
auto operator*(const T& mult, const Point3D<T> & point) {
    return point*mult;
}

#endif // BASE_POINT_H
