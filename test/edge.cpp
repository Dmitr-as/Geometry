#include <cassert>
#include "../base_geometry.h"

using namespace std;
using namespace gmtr;

int main()
{
    Point3DReal pt1(1,1,5), pt2(1,1,5), pt3(1,5,0);
    EdgeReal edge1(pt1, pt2), edge2(pt2, pt3);
    {
        assert(!edge1.isValid());
        assert(edge2.isValid());
    }
    {
        //Point3DReal p1(10,10,10);
        //assert( pln1.vectorDistance(p1) == p1.z() );
    }
    {

    }
    {

    }

    return 0;
}
