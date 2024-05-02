#ifndef BASE_GEOMETRY_H
#define BASE_GEOMETRY_H

//#include <math.h>
#include "base_point.h"

template<typename T>
struct Plane {
    using P3D = Point3D<T>;
private:
    P3D normal_;
    T distance_;

public:
    explicit Plane(const P3D& norm = P3D(), T dist = T(0) )
        : normal_(norm), distance_(dist) {}

    auto normal() const {return normal_;}
    auto distance() const {return distance_;}

    bool isValid() const {
        return !normal_.isNull();
    }

    T vectorDistance(const P3D& point) const { // векторизованное расстояние до точки, normal должен быть нормированным
        return P3D::dot(point, normal_) + distance_;
    }
    P3D projection(const P3D &point) const {
        return point - normal_ * vectorDistance(point);
    }
    static Plane createPlane(const P3D& pt1, const P3D& pt2, const P3D& pt3) {
        auto v1 = pt1 - pt2;
        auto v2 = pt1 - pt3;
        Plane plane(P3D::cross(v1, v2));
        if( plane.normal_.normalaze() ) {
            plane.distance_ = - P3D::dot(plane.normal_, pt1);
        }
        return plane;
    }
    static Plane createPlane( const P3D& pt1, const P3D& pt2, const P3D& pt3, const P3D& ptPositive ) {
        auto v1 = pt1 - pt2;
        auto v2 = pt1 - pt3;
        Plane plane(P3D::cross(v1, v2));
        if( plane.normal_.normalaze() ) {
            plane.distance_ = - P3D::dot(plane.normal_, pt1);
            if( !plane.isPositive( ptPositive ) ) {
                plane.normal_ *= -1.;
                plane.distance_ *= -1.;
            }
        }
        return plane;
    }

    bool isPositive(const P3D& point) const { // точка находится с положительной стороны от плоскости
        return vectorDistance(point) > 0;
    }
    bool isNegative(const P3D& point) const { // точка находится с положительной стороны от плоскости
        return vectorDistance(point) < 0;
    }

    friend std::ostream &operator<<(std::ostream &out, const Plane &plane) {
        out << "nrm = " << plane.normal_ << ", d = " << plane.distance_;
        return out;
    }
};

template<typename T>
struct Edge {
    using P3D = Point3D<T>;
private:
    P3D pt1_, pt2_, dir_;
public:
    explicit Edge(const P3D &point1, const P3D &point2)
        : pt1_(point1), pt2_(point2), dir_((pt2_ - pt1_).normalazed()) {}

    P3D point1() const {return pt1_;}
    P3D point2() const {return pt2_;}
    P3D direction() const {return dir_;}

    bool isValid() const {
        return pt1_ != pt2_;
    }

    P3D projection(const P3D &point) const {
        return pt1_ - dir_ * P3D::dot(dir_, pt1_ - point);
    }
    // проекция точки в системе координат ребра
    T alpha(const P3D &point) const { // prPt = pt1*alpha + pt2*(1-alpha)
        return alphaPrivate( projection(point) );
    }
    T alphaHalf(const P3D &point) const { // prPt = pt1*(alpha+0.5) + pt2*(0.5-alpha)
        return alphaPrivate( projection(point) ) - 0.5;
    }

    T distanceQuad(const P3D &point) const {
        auto prPt = projection(point);
        auto a = alpha(prPt);
        if( a >= 1 )
            return (pt1_ - point).normQuad();
        if( a <= 0 )
            return (pt2_ - point).normQuad();
        return (prPt - point).normQuad();
    }
    T distance(const P3D &point) const {
        return sqrt( distanceQuad(point) );
    }
    // дописать!!!
    T distance(const Edge &edge) {
        return T(0);
    }

private:
    T alphaPrivate(const P3D &projPoint) const {
        auto dirProject = pt2_ - projPoint;
        return P3D::dot(dirProject, dir_);
    }

public:
    friend std::ostream &operator<<(std::ostream &out, const Edge &edge) {
        out << "pt1 = " << edge.pt1_ << ", pt2 = " << edge.pt1_;
        return out;
    }
};

template<typename T>
struct Triangle {
    using P3D = Point3D<T>;
private:
    P3D pt1_, dir1_, dir2_, normal_;
public:
    explicit Triangle(const P3D &point1, const P3D &point2, const P3D &point3)
        : pt1_(point1), dir1_(point2 - point1), dir2_(point3 - point1) {
        normal_ = P3D::cross(dir1_, dir2_);
        normal_.normalaze();
    }

    P3D point1() const {return pt1_;}
    P3D point2() const {return dir1_ + pt1_;}
    P3D point3() const {return dir2_ + pt1_;}
    P3D normal() const {return normal_;}

    bool isValid() const {
        return !normal_.isNull();
    }

    P3D projection(const P3D &point) const {
        return point - normal_ * P3D::dot(normal_, point - pt1_);
    }

    std::pair<T,T> alpha(const P3D &point) const {
        return alphaPrivate( projection(point) );
    }
private:
    std::pair<T,T> alphaPrivate(const P3D &projPoint) const {
        auto delPt = projPoint - pt1_;
        //auto v2 = pt2 - pt1;
        //auto v3 = pt3 - pt1;
        auto ax2 = P3D::cross(normal_, dir2_);
        auto ax3 = P3D::cross(normal_, dir1_);
        auto alpha1 = P3D::dot(ax2, delPt) / P3D::dot(ax2, dir1_);
        auto alpha2 = P3D::dot(ax3, delPt) / P3D::dot(ax3, dir2_);
        return {alpha1, alpha2};
    }
};

template<typename T>
struct Sphere {
    using P3D = Point3D<T>;
private:
    P3D center_;
    T radius_;
public:
    explicit Sphere(const P3D& cntr, const T& rds)
        : center_(cntr), radius_(rds) {}

    P3D center() const {return center_;}
    T radius() const {return radius_;}

    bool isValid() const {
        return (center_ == center_) && radius_ > 0;
    }

    T vectorDistance(const P3D& point) const { // векторизованное расстояние до сферы
        return P3D::distance(center_, point) - radius_;
    }
    // сфера, описанная около тетраэдра
    static Sphere createSphere( const P3D& pt1, const P3D& pt2, const P3D& pt3, const P3D& pt4 ) {
        auto v2 = pt2 - pt1;
        auto v3 = pt3 - pt1;
        auto v4 = pt4 - pt1;
        auto pt1lengthQuad = pt1.normQuad();
        auto delta = T(0.5)/P3D::mix(v2,v3,v4);
        P3D rQuad( pt2.normQuad() - pt1lengthQuad, pt3.normQuad() - pt1lengthQuad, pt4.normQuad() - pt1lengthQuad );
        P3D ptdx( v2.x(), v3.x(), v4.x() );
        P3D ptdy( v2.y(), v3.y(), v4.y() );
        P3D ptdz( v2.z(), v3.z(), v4.z() );
        auto delta1 = P3D::mix(rQuad,ptdy,ptdz);
        auto delta2 = P3D::mix(ptdx,rQuad,ptdz);
        auto delta3 = P3D::mix(ptdx,ptdy,rQuad);
        auto pt0 = P3D(delta1, delta2, delta3)*delta;
        return Sphere( pt0, (pt0 - pt1).norm() );
    }
    // окружность, описанная около треугольника
    static Sphere createSphere( const P3D& pt1, const P3D& pt2, const P3D& pt3 ) {
        auto v2 = pt2 - pt1;
        auto v3 = pt3 - pt1;
        auto v2lengthQuad = v2.normQuad();
        auto v3lengthQuad = v3.normQuad();
        auto v2v3 = P3D::dot(v2, v3);
        auto delta = T(2)*(v2lengthQuad*v3lengthQuad - v2v3*v2v3);
        auto deltaA = v3lengthQuad*( v2lengthQuad - v2v3 );
        auto deltaB = v2lengthQuad*( v3lengthQuad - v2v3 );
        auto a = deltaA / delta;
        auto b = deltaB / delta;
        auto ptd = v3*b + v2*a;
        return Sphere( pt1 + ptd, ptd.norm() );
    }

    friend std::ostream &operator<<(std::ostream &out, const Sphere &sphere) {
        out << "cntr = " << sphere.center_ << ", r = " << sphere.radius_;
        return out;
    }
};

template<typename T>
bool goodTriangleNormals( const Point3D<T> &point1, const Point3D<T> &point2, const Point3D<T> &edgePoint1, const Point3D<T> &edgePoint2 ) {
    using P3D = Point3D<T>;
    auto edgeDir = edgePoint2 - edgePoint1;
    return P3D::dot( P3D::cross(edgeDir, edgePoint1 - point1), P3D::cross(edgePoint1 - point2, edgeDir) ) > 0;
}

#define Point3DReal Point3D<double>
#define PlaneReal Plane<double>
#define EdgeReal Edge<double>
#define TriangleReal Triangle<double>
#define SphereReal Sphere<double>


#endif // BASE_GEOMETRY_H
