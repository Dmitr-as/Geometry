#pragma once

#include "base_point.h"
#include "base_geometry.h"
#include <tuple>
#include <execution>

namespace gmtr {

template<typename T>
class CorrelationMatrix {
    using P3D = Point3D<T>;
    P3D m_e;
    P3D m_c;
public:
    CorrelationMatrix(const P3D& energy, const P3D& covariance)
        : m_e(energy)
        , m_c(covariance) // (Rxy, Ryz, Rxz)
    {}

    T determinant() const noexcept {
        return m_e.x()*m_e.y()*m_e.z() + 2*m_c.x()*m_c.y()*m_c.z()
               - m_e.x()*m_c.y()*m_c.z() - m_e.y()*m_c.x()*m_c.z() - m_e.z()*m_c.x()*m_c.y();
    }

    T minor(int i) const noexcept {
        switch (i) {
        case 0:
            return m_e.y()*m_e.z() - m_c.y()*m_c.y();
        case 1:
            return m_e.x()*m_e.z() - m_c.z()*m_c.z();
        case 2:
            return m_e.x()*m_e.y() - m_c.x()*m_c.x();
        }
        return 0;
    }

    static T determinant(const P3D& energy, const P3D& covariance) noexcept {
        return energy.x()*energy.y()*energy.z() + 2*covariance.x()*covariance.y()*covariance.z()
               - energy.x()*covariance.y()*covariance.z()
               - energy.y()*covariance.x()*covariance.z()
               - energy.z()*covariance.x()*covariance.y();
    }

    static T minor1(const P3D& energy, const P3D& covariance) noexcept {
        return energy.y()*energy.z() - covariance.y()*covariance.y();
    }

    static T minor2(const P3D& energy, const P3D& covariance) noexcept {
        return energy.x()*energy.z() - covariance.z()*covariance.z();
    }

    static T minor3(const P3D& energy, const P3D& covariance) noexcept {
        return energy.x()*energy.y() - covariance.x()*covariance.x();
    }

};

template<typename Iterator>
auto statistic(Iterator begin, Iterator end)
{
    using P3D = std::iterator_traits<Iterator>::value_type;
    auto pointMiddle = std::accumulate( /*std::execution::par_unseq,*/ begin, end, P3D() ) / std::distance(begin, end);
    auto pointEnergy = std::accumulate( /*std::execution::par_unseq,*/ begin, end, P3D(), [pointMiddle](const P3D &s, const P3D &v) {
        auto vm = v - pointMiddle;
        return s + vm*vm;
    } );
    auto covariance = std::accumulate( /*std::execution::par_unseq,*/ begin, end, P3D(), [pointMiddle](const P3D &s, const P3D &v) {
        auto vm = v - pointMiddle;
        return s + P3D(vm.x()*vm.y(), vm.y()*vm.z(), vm.z()*vm.x());
    } );
    return std::tuple{pointMiddle, pointEnergy, covariance};
}

template<typename Iterator>
auto line_regression(const Iterator begin, const Iterator end)
{
    using P3D = std::iterator_traits<Iterator>::value_type;
    using T = P3D::value_type;
    auto [pointM, pointEnergy, covariance] = statistic(begin, end);
    auto Rxy = covariance.x(), Ryz = covariance.y(), Rxz = covariance.z();
    auto det_xy = pointEnergy.x() * pointEnergy.y() - Rxy*Rxy;
    auto det_yz = pointEnergy.z() * pointEnergy.y() - Ryz*Ryz;
    auto det_xz = pointEnergy.x() * pointEnergy.z() - Rxz*Rxz;
    /*
    auto det_xy_norm = det_xy / pointEnergy.x() / pointEnergy.y(),
        det_yz_norm = det_yz / pointEnergy.z() / pointEnergy.y(),
        det_xz_norm = det_xz / pointEnergy.x() / pointEnergy.z();
    auto del_x_rxy = pointEnergy.x() * det_xy_norm * (pointEnergy.y() + pointEnergy.z()),
        del_y_rxy = pointEnergy.y() * det_xy_norm * (pointEnergy.x() + pointEnergy.z()),

        del_y_ryz = pointEnergy.y() * det_yz_norm * (pointEnergy.x() + pointEnergy.z()),
        del_z_ryz = pointEnergy.z() * det_yz_norm * (pointEnergy.x() + pointEnergy.y()),

        del_x_rxz = pointEnergy.x() * det_xz_norm * (pointEnergy.y() + pointEnergy.z()),
        del_z_rxz = pointEnergy.z() * det_xz_norm * (pointEnergy.x() + pointEnergy.y());
    */
    CorrelationMatrix<T> matrix(pointEnergy, covariance);
    auto det = matrix.determinant();
    std::vector errv{det/matrix.minor(0), det/matrix.minor(1), det/matrix.minor(2)};
    auto min = std::distance(errv.begin(), std::min_element(errv.begin(), errv.end()));
    P3D normal;
    switch (min) {
    case 0: // del_x_rxy >= del_z_rxz && del_y_rxy >= del_z_ryz
        normal = P3D((Ryz * Rxy - Rxz * pointEnergy.y()) / det_xy, (Rxz * Rxy - Ryz * pointEnergy.x()) / det_xy, 1);
        //normal.setY();
        //normal.setZ(1); // z == 1
        break;
    case 1: // del_x_rxz >= del_y_ryz && del_z_rxz >= del_y_rxy
        normal = P3D((Ryz * Rxz - Rxy * pointEnergy.z()) / det_xz, 1, (Rxy * Rxz - Ryz * pointEnergy.x()) / det_xz);
        //normal.setY(1); // y == 1
        //normal.setZ();
        break;
    case 2: // if( del_y_ryz > del_x_rxz && del_z_ryz > del_x_rxy )
        normal = P3D(1, (Rxz * Ryz - Rxy * pointEnergy.z()) / det_yz, (Rxy * Ryz - Rxz * pointEnergy.y()) / det_yz); // x == 1
        //normal.setY();
        //normal.setZ();
        break;
    }
    normal.normalize();
    return Plane<T>( normal, -P3D::dot(normal, pointM) );
}

template<typename Iterator>
auto vector_regression(Iterator begin, Iterator end)
{
    using P3D = std::iterator_traits<Iterator>::value_type;
    using T = P3D::value_type;

    auto [pointM, pointEnergy, covariance] = statistic(begin, end);
    if(pointEnergy.isNull()) {
        //qDebug() << "vectorRegression" << "energy is null";
        return Line<T>();
    }
    auto Rxy = covariance.x();
    auto Ryz = covariance.y();
    auto Rxz = covariance.z();
    auto covarianceQuad = covariance*covariance;
    auto RxyQuad = covarianceQuad.x();
    auto RyzQuad = covarianceQuad.y();
    auto RxzQuad = covarianceQuad.z();
    auto pointEnergyQuad = pointEnergy*pointEnergy;

    auto d_sum = pointEnergy.x() + pointEnergy.y() + pointEnergy.z();
    auto del_x_y = pointEnergyQuad.x()*(pointEnergy.x()+pointEnergy.y()) + (pointEnergy.x() + d_sum)*RxzQuad;
    auto del_x_z = pointEnergyQuad.x()*(pointEnergy.x()+pointEnergy.z()) + (pointEnergy.x() + d_sum)*RxyQuad;

    auto del_y_x = pointEnergyQuad.y()*(pointEnergy.y()+pointEnergy.x()) + (pointEnergy.y() + d_sum)*RyzQuad;
    auto del_y_z = pointEnergyQuad.y()*(pointEnergy.y()+pointEnergy.z()) + (pointEnergy.y() + d_sum)*RxyQuad;

    auto del_z_x = pointEnergyQuad.z()*(pointEnergy.z()+pointEnergy.x()) + (pointEnergy.z() + d_sum)*RyzQuad;
    auto del_z_y = pointEnergyQuad.z()*(pointEnergy.z()+pointEnergy.y()) + (pointEnergy.z() + d_sum)*RxzQuad;
    P3D dir;
    if(del_z_x >= del_x_z && del_z_y >= del_y_z) { // z = 1
        auto del = (pointEnergy.x() + pointEnergy.z())*(pointEnergy.y() + pointEnergy.z()) - Rxy*Rxy;
        auto del1 = Rxz*(pointEnergy.x() + pointEnergy.z()) + Rxy*Ryz;
        auto del2 = Ryz*(pointEnergy.y() + pointEnergy.z()) + Rxy*Rxz;
        dir = P3D( del1/del, del2/del, 1 );
    }
    else if(del_y_x >= del_x_y && del_y_z >= del_z_y) { // y = 1
        auto del = (pointEnergy.x() + pointEnergy.y())*(pointEnergy.z() + pointEnergy.y()) - Rxz*Rxz;
        auto del1 = Rxy*(pointEnergy.x() + pointEnergy.y()) + Rxz*Ryz;
        auto del2 = Ryz*(pointEnergy.z() + pointEnergy.y()) + Rxy*Rxz;
        dir = P3D( del1/del, 1, del2/del );
    }
    else { // if(del_x_y >= del_y_x && del_x_z >= del_z_x) { // x = 1
        auto del = (pointEnergy.y() + pointEnergy.x())*(pointEnergy.z() + pointEnergy.x()) - Ryz*Ryz;
        auto del1 = Rxy*(pointEnergy.y() + pointEnergy.x()) + Ryz*Rxz;
        auto del2 = Rxz*(pointEnergy.z() + pointEnergy.x()) + Ryz*Rxy;
        dir = P3D( 1, del1/del, del2/del );
    }
    //else {
        //qDebug() << "vectorRegression" << "other variant";
        //return Line<T>();
    //}

    dir.normalize();
    return Line<T>(pointM, dir);
}

// Sphere
template<typename Iterator>
auto sphere_regression(Iterator begin, Iterator end)
{
    using P3D = std::iterator_traits<Iterator>::value_type;
    using T = P3D::value_type;

    struct RES3 {
        P3D energy;
        P3D m3;
        P3D r;
    };

    auto l = std::distance(begin, end);
    auto pointMiddle = std::accumulate( begin, end, P3D() ) / l;
    auto params = std::accumulate( begin, end, RES3{P3D(), P3D(), P3D()}, [&pointMiddle](const RES3 &s, const P3D &v) {
        auto vm = v - pointMiddle;
        auto e = vm*vm;
        return RES3{
            s.energy + e,
            s.m3 + vm*e + vm*Point3f(e.y(), e.z(), e.x()) + vm*Point3f(e.z(), e.x(), e.y()),
            s.r + vm*Point3f(vm.y(), vm.z(), vm.x())
        };
    } );
    auto Rxy = params.r.x();
    auto Ryz = params.r.y();
    auto Rxz = params.r.z();
    P3D p1(params.energy.x(), Rxy, Rxz), p2(Rxy, params.energy.y(), Ryz), p3(Rxz, Ryz, params.energy.x());
    auto del = 2 * P3D::mix(p1, p2, p3);

    P3D r1(params.energy.y()*params.energy.z() - Ryz*Ryz, Rxz*Ryz - params.energy.z()*Rxy, Rxy*Ryz - params.energy.y()*Rxz),
        r2(Rxz*Ryz - params.energy.z()*Rxy, params.energy.x()*params.energy.z() - Rxz*Rxz, Rxy*Rxz - params.energy.x()*Ryz),
        r3(Rxy*Ryz - params.energy.y()*Rxz, Rxy*Rxz - params.energy.x()*Ryz, params.energy.x()*params.energy.y() - Rxy*Rxy);
    auto P0 = P3D(dot(params.m3, r1), dot(params.m3, r2), dot(params.m3, r3)) / del;

    return Sphere<T>( pointMiddle + P0, sqrtf( ( params.energy.x() + params.energy.y() + params.energy.z() ) / l + P0.lengthSquared() ) );
}

}
