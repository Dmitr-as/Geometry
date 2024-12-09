#include <cassert>
#include "../statistic.h"

using namespace std;
using namespace gmtr;




int main()
{
    Point3DReal dir1(5, 1, 50), dir2(1, 5, 1), start(1,1,1);
    auto normal = Point3DReal::cross(dir1, dir2);
    cout << "normal0 " << normal.normalized() << endl;

    vector<Point3DReal> points;
    points.reserve(100);
    for(int i = 0; i < 10; ++i) {
        for(int j = 0; j < 10; ++j) {
            points.push_back( start + i * dir1 + j *dir2 );
        }
    }
    auto line = line_regression(points.begin(), points.end());
    cout << "normal " << line.normal() << endl;
    return 0;
}
