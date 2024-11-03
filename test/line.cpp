#include <cassert>
#include "../base_geometry.h"

using namespace std;
using namespace gmtr;

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
        //Point3DReal p1(10,10,10);
        //assert( pln1.vectorDistance(p1) == p1.z() );
    }
    {

    }
    {

    }

    return 0;
}
