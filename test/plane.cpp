#include <cassert>
#include "../base_geometry.h"

using namespace std;
using namespace gmtr;

int main()
{
    PlaneReal pln, pln1(Point3DReal(0,0,1)), pln2(Point3DReal(0,1,0)), pln3(Point3DReal(1,0,0));
    {
        assert(!pln.isValid());
        assert(pln1.isValid() && pln2.isValid() && pln3.isValid());
    }
    {
        Point3DReal p1(10,10,10);
        assert( pln1.vectorDistance(p1) == p1.z() );

        assert( pln2.projection(p1) == Point3DReal(10, 0, 10) );
        assert( pln3.projection(p1) == Point3DReal(0, 10, 10) );
    }
    {
        Point3DReal p1(10,10,10), p2(10,10,10), p3(10,10,10), p4(10,-10,10),p5(10,10,-10);
        assert(!PlaneReal::createPlane(p1,p2,p3).isValid());
        assert(!PlaneReal::createPlane(p1,-p2,-p3).isValid());
        assert(PlaneReal::createPlane(p1,p4,p5).isValid());
    }
    {
        Point3DReal p1(0,0,121), p2(100,0,0), p3(0,100,0), p4(0,0,0);
        auto pln1 = PlaneReal::createPlane(p1,p2,p3);
        auto pln2 = PlaneReal::createPlane(p1,p2,p4);
        auto pln3 = PlaneReal::createPlane(p1,p4,p3);
        auto ptCross = PlaneReal::cross(pln1, pln2, pln3);
        assert( Point3DReal::distance(ptCross, p1) == 0 );
    }
    {
        assert(pln1.isPositive(Point3DReal{-4,4,10}));
        assert(pln1.isNegative(Point3DReal{2,-2,-10}));
    }

    return 0;
}
