#ifndef ROTATION_H
#define ROTATION_H

#include "base_point.h"

namespace gmtr {
/**
 * @brief Последовательный поворот вокруг 3 осей Z, X, Z''
 */
template <std::floating_point T>
class Rotation3D
{
    using P3D = Point3D<T>;
    P3D strMatrix[3]; // строки матрицы поворота
public:
    constexpr Rotation3D(T alpha, T beta, T gamma) {
        auto sinAlpha = sin(alpha);
        auto cosAlpha = cos(alpha);
        auto sinBeta = sin(beta);
        auto cosBeta = cos(beta);
        auto sinGamma = sin(gamma);
        auto cosGamma = cos(gamma);

        strMatrix[0] = {cosAlpha*cosGamma - sinAlpha*cosBeta*sinGamma, -cosAlpha*sinGamma - sinAlpha*cosBeta*cosGamma, sinAlpha*sinBeta};
        strMatrix[1] = {sinAlpha*cosGamma + cosAlpha*cosBeta*sinGamma, -sinAlpha*sinGamma - cosAlpha*cosBeta*cosGamma, -cosAlpha*sinBeta};
        strMatrix[2] = {sinBeta*sinGamma, sinBeta*cosGamma, cosBeta};
    }

    P3D rotate(const P3D& point) const {
        return P3D{P3D::dot(point, strMatrix[0]), P3D::dot(point, strMatrix[1]), P3D::dot(point, strMatrix[2])};
    }
};

using Rotation3DReal = Rotation3D<double>;

}
#endif // ROTATION_H
