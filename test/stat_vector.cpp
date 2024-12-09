#include <cassert>
#include "../statistic.h"

using namespace std;
using namespace gmtr;

int main()
{
    Point3DReal dir(5, 1, 5), start(1,1,1);
    cout << "dir0 " << dir.normalized() << endl;

    vector<Point3DReal> points;
    points.reserve(25);
    for(int i = 0; i < 25; ++i) {
        points.push_back( start + i * dir );
    }
    auto line = vector_regression(points.begin(), points.end());
    cout << "dir " << line.direction() << endl;
    return 0;
}
