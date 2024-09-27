#include <cassert>
#include "../base_geometry.h"

using namespace std;
using namespace gmtr;

int main()
{
    {
        Point3DReal pt;
        assert(pt.isNull());
        Point3DReal pt1(1, 1, 1), pt2(1, 1, 1);
        assert(!pt1.isNull());
        assert(!pt2.isNull());
        assert(pt != pt1);
        assert(pt != pt2);
        assert(pt1 == pt2);
        assert(pt1 + pt2 == pt2 + pt1);
        {
            pt = pt1 * 0;
            assert(pt.isNull());
            pt = 0. * pt1;
            assert(pt.isNull());
        }
        {
            assert(pt1 == pt1 * 1);
            assert(pt1 == 1 * pt1);
            assert((pt1 * 0).isNull());
            assert((0 * pt1).isNull());
        }
        assert(Point3DReal::dot(pt1, pt2) == 3.);
        assert(Point3DReal::cross(pt1, pt2).isNull());
    }
    {
        Point3DReal pt1(2., 0., 0.), pt2(0., 3., 0.);
        auto pt3 = Point3DReal::cross(pt1, pt2);
        assert(Point3DReal(0., 0., 6.) == pt3);
    }
    {
        Point3DReal pt1(1, 0, 0), pt2(0, 1, 0), pt3(0, 0, 1);
        assert(1 == Point3DReal::mix(pt1, pt2, pt3));
    }
    {
        Point3DReal pt1(1, 2, 3);
        auto pt2 = pt1;
        assert(pt2 == pt1.shiftLeft().shiftLeft().shiftLeft());
        assert(pt2 == pt1.shiftRight().shiftRight().shiftRight());
        assert(pt2 != pt1.shiftLeft());
        assert(pt2 != pt1.shiftLeft());
        assert(pt2 != pt1.shiftRight());
    }

    return 0;
}
