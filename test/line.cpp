#include <cassert>
#include <concepts>
#include "../base_geometry.h"

using namespace std;
using namespace gmtr;

template <std::floating_point T>
static inline bool fuzzyCompare(T p1, T p2)
{
    return (std::fabs(p1 - p2) * 1000000000000. <= std::min(std::fabs(p1), std::fabs(p2)));
}

int main()
{
    Point3DReal pt1(1,1,5), pt2(1,1,5), pt3(1,5,0);
    LineReal line1(pt1, pt2), line2(pt2, pt3);
    auto lineTrue = LineReal::createLine(pt1, pt3);
    auto lineFalse = LineReal::createLine(pt1, pt2);
    {
        assert(line1.isValid());
        assert(line2.isValid());

        assert(lineTrue.isValid());
        assert(!lineFalse.isValid());
    }
    {

    }
    {
        Point3DReal p1(10,10,10);
        auto prPt = line1.projection(p1);
        auto dist = Point3DReal::distance(prPt, p1);
        assert( fuzzyCompare(line1.distance(p1), dist) );
    }
    {
        double z1 = -4, z2 = 6;
        Point3DReal p1(1,1,z1), p2(1,10,z1), p3(5,-1,z2), p4(1,10,z2), p5(1,1,z2), p6(1,10,z2);
        auto line1 = LineReal::createLine(p1, p2), line2 = LineReal::createLine(p3, p4), line3 = LineReal::createLine(p6, p5);
        //LineReal line1(p1, Point3DReal{2,1,0}), line2(p3, Point3DReal{1,5,0});
        //cout << line1.distance(line2) << endl;
        assert( line1.distance(line2) == std::fabs(z2-z1) );
        assert( line1.distance(line3) == std::fabs(z2-z1) );
    }
    return 0;
}
