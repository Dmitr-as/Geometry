#include <cassert>
#include "../statistic.h"
#include <math.h>

using namespace std;
using namespace gmtr;

int main()
{
    Point3DReal dir1(5, 1, 50), dir2(1, 5, 1), start(1,1,1);
    dir1.normalize();
    dir2.normalize();
    auto normal = Point3DReal::cross(dir1, dir2);
    cout << "normal0 " << normal.normalized() << endl;
    cout << "dir0 " << dir1.normalized() << endl;

    vector<Point3DReal> points;
    points.reserve(1000);
    for(int i = 0; i < 100; ++i) {
        points.push_back( start + i * dir1 );
    }
    for(int i = 0; i < 100; ++i) {
        points.push_back( start + dir2 + i * dir1 );
    }
    auto lineNorm = line_regression(points.begin(), points.end());
    auto lineLongitudinal = vector_regression(points.begin(), points.end());
    cout << "normal " << lineNorm.normal() << endl;
    cout << "vector " << lineLongitudinal.direction() << endl;
    cout << "cross " << acos(Point3DReal::dot( lineLongitudinal.direction(), lineNorm.normal() )) * 180 / M_PI - 90 << endl;
    return 0;
}
