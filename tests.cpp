#include <cassert>
#include "base_geometry.h"

using namespace std;
using namespace gmtr;

int main(int argc, char *argv[])
{
    return 0; // TODO
    {
        cout << "tests Point3D" << endl;
        Point3DReal pt1(1,1,0), pt2(0,-1,1), pt3;

        cout << "normalized " << pt1.normalized() << endl;
        cout << "normalized " << pt3.normalized() << endl;
        assert( !pt1.isNull() );
        assert( !pt2.isNull() );
        assert( pt3.isNull() );
        pt3 = pt1 + pt2;

        cout << "point1 = " << pt1 << " point2 = " << pt2 << " point3 = " << pt3 << endl;
        assert(pt1 + pt2 == Point3DReal(1,0,1));
        assert(Point3DReal::dot(pt1, pt2) == -1);
        cout << Point3DReal::distance(pt1, pt2) << endl;

        auto pt11 = pt1;
        pt11 += pt2;
        assert(pt11 == pt1 + pt2);
        pt11 -= pt2;

        pt11 = 5. * pt11;
        pt11 *= 5.;
        cout << pt11 << endl;
        pt11 /= 5.;
        cout << pt11 << endl;

        auto pt11m = -pt11;
        assert(pt11m == -pt11);

        auto ptFloat = pt11m.cast<float>();
        cout << "point float " << ptFloat << endl;

        auto ptLDouble = pt11.cast<long double>();
        cout << "point long double " << ptLDouble << endl;

        cout << sizeof(ptFloat) << " " << sizeof(pt11) << " " << sizeof(ptLDouble) << endl;

        Point3DReal ptForShift(1,2,4);
        ptForShift.shiftLeft();
        cout << " shift left " << ptForShift << endl;
        ptForShift.shiftRight();
        cout << " shift right " << ptForShift << endl;
    }

    {
        cout << "----- ----- ----- " << endl;
        Point3DReal pt1(1.5,0.3,0), pt2(0,1.5,0.3), pt3(0.3,0,-1.5);
        PlaneReal plane( PlaneReal::createPlane(pt1, pt2, pt3) );
        cout << "plane: "<< plane << endl;
        cout << "plane dist = "<< plane.vectorDistance( pt1 + (pt1 - pt2)*2.5 + (pt1 - pt3)*2.5 ) << endl;
        cout << "plane dist = "<< plane.vectorDistance( pt1 + (pt1 - pt2)*5.5 - (pt1 - pt3)*2.5 ) << endl;
        cout << "plane dist = "<< plane.vectorDistance( pt1 - (pt1 - pt2)*7.5 + (pt1 - pt3)*5.5 ) << endl;
        //static_assert( plane.vectorDistance(pt1) == 0, "v" );
        cout << "eps " << numeric_limits<double>::epsilon() << endl;
        //assert(plane.vectorDistance(pt1) == 0);
        //assert(plane.vectorDistance(pt2) == 0);
        //assert(plane.vectorDistance(pt3) == 0);
    }

    {
        cout << "----- ----- ----- " << endl;
        Point3DReal pt1(1,0,0), pt2(0,1,0), pt3(0,0,1);
        PlaneReal plane = PlaneReal::createPlane(pt1, pt2, pt3);
        cout << "plane "<< plane << endl;
        if(plane.isPositive( pt2 ))
            cout << "positive" << endl;
        else
            cout << "negative" << endl;
        assert(plane.isValid());
    }

    {
        Point3DReal pt1(0,0,0), pt2(10,0,0);
        EdgeReal edge1(pt1, pt2);
        Point3DReal pt3(pt1);
    }

    {
        Point3DReal pt1(5,0,0), pt2(5,3,0), ptAngle90Mines(5,3,9.99), ptAngle90(5,3,10), ptAngle90Plus(5,3,10.01), ptNext(5,3,15);
        Point3DReal ept1(0,0,10), ept2(10,0,10);
        assert( !goodTriangleNormals(pt1, pt1, ept1, ept2) );
        assert( !goodTriangleNormals(pt1, pt2, ept1, ept2) );
        assert( !goodTriangleNormals(pt1, ptAngle90, ept1, ept2) );
        assert( !goodTriangleNormals(pt1, ptAngle90Mines, ept1, ept2) );
        assert( goodTriangleNormals(pt1, ptAngle90Plus, ept1, ept2) );
        assert( goodTriangleNormals(pt1, ptNext, ept1, ept2) );
    }

    {
        cout << "----- ----- ----- " << endl;
        Point3DReal pt1(5.,0,0), pt2(0,5.,0), pt3(0,0,50.), pt4(5.,15.,10.);
        auto plane = PlaneReal::createPlane(pt1, pt2, pt3);
        auto arc = SphereReal::createSphere(pt1, pt2, pt3);
        cout << "arc: " << arc << endl;
        auto err = plane.vectorDistance(arc.center());
        cout << "plane = " << plane << endl;
        cout << "plane value = " << err << endl;
        //assert( err == 0 );
        auto r1q = (pt1 - arc.center()).norm();
        auto r2q = (pt2 - arc.center()).norm();
        auto r3q = (pt3 - arc.center()).norm();
        //assert( arc.radius() == r1q );
        //assert( arc.radius() == r2q );
        //assert( arc.radius() == r3q );
        //Point3DReal pt3_1(0,3,7);
        auto arc2 = SphereReal::createSphere(pt1, pt2, pt1*0.5 + pt2*0.5 + Point3DReal(0.001,0.001,0.001)*0);
        cout << "arc2 " << arc2.radius() << " " << arc2.isValid() << endl;
        auto sphere = SphereReal::createSphere(pt1, pt2, pt3, pt4);
        cout << "sphere center " << sphere.center() << ", sphere rs " << sphere.radius()
             << " " << (pt4 - sphere.center()).norm()
             << " " << (pt3 - sphere.center()).norm()
             << " " << (pt2 - sphere.center()).norm()
             << " " << (pt1 - sphere.center()).norm() << endl;
    }

    PlaneReal pln, pln1(Point3DReal(0,0,1));
    Point3DReal p1(10,10,10);
    //assert( pln1.projection(p1) == Point3DReal(10, 0, 10) );
    cout << pln1.projection(p1) << endl;

    return EXIT_SUCCESS;
}
