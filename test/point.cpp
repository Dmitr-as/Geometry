#include <cassert>
#include "../base_geometry.h"

using namespace std;
using namespace gmtr;

int main()
{
    Point3DReal pt;
    assert( pt.isNull() );
    Point3DReal pt1(1,1,1), pt2(1,1,1);
    assert( !pt1.isNull() );
    assert( !pt2.isNull() );
    assert(pt1 == pt2);
    assert(pt1 + pt2 == pt2 + pt1);

    pt = pt1*0;
    assert( pt.isNull() );
    pt = 0.*pt1;
    assert( pt.isNull() );

    assert( Point3DReal::dot(pt1, pt2) == 3. );
    assert( Point3DReal::cross(pt1, pt2).isNull() );
    {
        Point3DReal pt1(1.6, 0.2, 0.5), pt2(0.3, 1.1, 0.9);
        auto pt3 = Point3DReal::cross(pt1, pt2);
        //auto pt1cross = Point3DReal::cross(pt2, pt3);
        //cout << pt1cross << endl;
        //assert( pt1cross == pt1 );
    }
    return 0;
}
