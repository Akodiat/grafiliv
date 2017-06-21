#ifndef LINALG_H
#define LINALG_H

using namespace std;

//Three-dimensional vector cross product
xyz cross(xyz a, xyz b){
    return make_xyz(
        -a.z*b.y + a.y*b.z,
        a.z*b.x - a.x*b.z,
        -a.y*b.x + a.x*b.y
    );
}

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
