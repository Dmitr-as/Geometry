#ifndef ROTATION_H
#define ROTATION_H

#include "base_point.h"

namespace gmtr {

template <std::floating_point T>
class Rotation3D
{
    using P3D = Point3D<T>;
    P3D strv[3];
public:
    Rotation3D(T alpha, T beta, T gamma)
    {

    }
};

}
#endif // ROTATION_H
