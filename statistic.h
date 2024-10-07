#pragma once

#include "base_point.h"
#include <tuple>

namespace gmtr {

template<typename T, typename Iterator>
std::tuple< Point3D<T>, Point3D<T>, Point3D<T> > statistic(const Iterator begin, const Iterator end) {
    using P3D = Point3D<T>;
    auto pointMiddle = std::accumulate( begin, end, P3D() ) / std::distance(begin, end);
    auto pointEnergy = std::accumulate( begin, end, P3D(), [pointMiddle](const P3D &s, const P3D &v) {
        auto vm = v - pointMiddle;
        return s + vm*vm;
    } );
    auto covariance = std::accumulate( begin, end, P3D(), [pointMiddle](const P3D &s, const P3D &v) {
        auto vm = v - pointMiddle;
        return s + Point3f(vm.x()*vm.y(), vm.y()*vm.z(), vm.z()*vm.x());
    } );
    return {pointMiddle, pointEnergy, covariance};
}
}
