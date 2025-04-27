#include <cassert>
#include "../base_geometry.h"
#include "../rotation.h"

using namespace std;
using namespace gmtr;
using namespace numbers;

int main()
{
    Point3DReal simpleVector{1,0,0};
    ;
    cout << Rotation3DReal(0,0,0).rotate(simpleVector);
    cout << Rotation3DReal(pi/2,0,0).rotate(simpleVector);
    cout << Rotation3DReal(0,pi/2,0).rotate(simpleVector);
    cout << Rotation3DReal(0,0,pi/2).rotate(simpleVector) << endl;
    int n = 10;
    for(int ix = -n; ix < n; ++ix) {
        auto alp = numbers::pi * ix / n;
        for(int iy = -n; iy < n; ++iy) {
            auto bet = numbers::pi * iy / n;
            for(int iz = -n; iz < n; ++iz) {
                auto gam = numbers::pi * iz / n;
                Rotation3DReal rt(alp, bet, gam);
                auto resVector = rt.rotate(simpleVector);

            }
        }
    }


    return 0;
}
