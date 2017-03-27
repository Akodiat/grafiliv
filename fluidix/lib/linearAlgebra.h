#ifndef LINALG_H
#define LINALG_H

using namespace std;

class Matrix3
{
public:
    Matrix3(xyz a, xyz b, xyz c) : a(a), b(b), c(c)
    {}

    xyz dot(xyz v) {
        return make_xyz(
            a.x*v.x + a.y*v.y + a.z*v.z,
            b.x*v.x + b.y*v.y + b.z*v.z,
            c.x*v.x + c.y*v.y + c.z*v.z
        );
    }
private:
    xyz a, b, c;
};

#endif
