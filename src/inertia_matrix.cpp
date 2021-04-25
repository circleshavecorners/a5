#include <inertia_matrix.h>
#include <cassert>

//compute inertia matrix and volume by integrating on surfaces
void inertia_matrix(Eigen::Matrix3d &I, Eigen::Vector3d &center, double &mass, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, double density)
{
    mass = 0;
    center.setZero();
    double inv6 = (1. / 6.);
    for (int f = 0; f < F.rows(); ++f)
    {
        Eigen::Vector3d v0 = V.row(F(f, 0));
        Eigen::Vector3d v1 = V.row(F(f, 1));
        Eigen::Vector3d v2 = V.row(F(f, 2));
        Eigen::Vector3d n = (v1 - v0).cross(v2 - v0);
        mass += inv6 * density * n(0) * (v0(0) + v1(0), v2(0));
    }
    double invmass = (1. / mass);
    double inv12 = 1. / 12.;
    double inv20 = 1. / 20.;
    double inv60 = 1. / 60.;
    double Ixx = 0;
    double Iyy = 0;
    double Izz = 0;
    double Ixy = 0;
    double Iyz = 0;
    double Izx = 0;
    for (int f = 0; f < F.rows(); ++f)
    {
        Eigen::Vector3d v0 = V.row(F(f, 0));
        Eigen::Vector3d v1 = V.row(F(f, 1));
        Eigen::Vector3d v2 = V.row(F(f, 2));
        //double x0 = v0(0), y0 = v0(1), z0 = v0(2);
        Eigen::Vector3d n = (v1 - v0).cross(v2 - v0);
        center(0) += 0.5 * invmass * density * n(0) * (inv12 * (v0(0) * v0(0) + v1(0) * v1(0) + v2(0) * v2(0)) + inv12 * (v0(0) * v1(0) + v0(0) * v2(0) + v1(0) * v2(0)));
        center(1) += 0.5 * invmass * density * n(1) * (inv12 * (v0(1) * v0(1) + v1(1) * v1(1) + v2(1) * v2(1)) + inv12 * (v0(1) * v1(1) + v0(1) * v2(1) + v1(1) * v2(1)));
        center(2) += 0.5 * invmass * density * n(2) * (inv12 * (v0(2) * v0(2) + v1(2) * v1(2) + v2(2) * v2(2)) + inv12 * (v0(2) * v1(2) + v0(2) * v2(2) + v1(2) * v2(2)));
        double tmpx2 = density * n(0) * inv60 * (std::pow(v0(0), 3) + std::pow(v1(0), 3) + std::pow(v2(0), 3) + std::pow(v0(0), 2) * v1(0) + std::pow(v0(0), 2) * v2(0) + std::pow(v1(0), 2) * v0(0) + std::pow(v1(0), 2) * v2(0) + std::pow(v2(0), 2) * v0(0) + std::pow(v2(0), 2) * v1(0) + v0(0) * v1(0) * v2(0));
        double tmpy2 = density * n(1) * inv60 * (std::pow(v0(1), 3) + std::pow(v1(1), 3) + std::pow(v2(1), 3) + std::pow(v0(1), 2) * v1(1) + std::pow(v0(1), 2) * v2(1) + std::pow(v1(1), 2) * v0(1) + std::pow(v1(1), 2) * v2(1) + std::pow(v2(1), 2) * v0(1) + std::pow(v2(1), 2) * v1(1) + v0(1) * v1(1) * v2(1));
        double tmpz2 = density * n(2) * inv60 * (std::pow(v0(2), 3) + std::pow(v1(2), 3) + std::pow(v2(2), 3) + std::pow(v0(2), 2) * v1(2) + std::pow(v0(2), 2) * v2(2) + std::pow(v1(2), 2) * v0(2) + std::pow(v1(2), 2) * v2(2) + std::pow(v2(2), 2) * v0(2) + std::pow(v2(2), 2) * v1(2) + v0(2) * v1(2) * v2(2));
        Ixx += tmpy2 + tmpz2;
        Iyy += tmpx2 + tmpz2;
        Izz += tmpx2 + tmpy2;
        Ixy += density * n(0) * 0.5 * (inv20 * (std::pow(v0(0), 2) * v0(1) + std::pow(v1(0), 2) * v1(1) + std::pow(v2(0), 2) * v2(1)) + inv60 * (std::pow(v0(0), 2) * v1(1) + 2 * v0(0) * v1(0) * v0(1) + std::pow(v0(0), 2) * v2(1) + 2 * v0(0) * v2(0) * v0(1) + std::pow(v1(0), 2) * v0(1) + 2 * v0(0) * v1(0) * v1(1) + std::pow(v1(0), 2) * v2(1) + 2 * v1(0) * v2(0) * v1(1) + std::pow(v2(0), 2) * v0(1) + 2 * v0(0) * v2(0) * v2(1) + std::pow(v2(0), 2) * v1(1) + 2 * v1(0) * v2(0) * v2(1) + v0(0) * v1(0) * v2(1) + v0(0) * v2(0) * v1(1) + v1(0) * v2(0) * v0(1)));
        Iyz += density * n(1) * 0.5 * (inv20 * (std::pow(v0(1), 2) * v0(2) + std::pow(v1(1), 2) * v1(2) + std::pow(v2(1), 2) * v2(2)) + inv60 * (std::pow(v0(1), 2) * v1(2) + 2 * v0(1) * v1(1) * v0(2) + std::pow(v0(1), 2) * v2(2) + 2 * v0(1) * v2(1) * v0(2) + std::pow(v1(1), 2) * v0(2) + 2 * v0(1) * v1(1) * v1(2) + std::pow(v1(1), 2) * v2(2) + 2 * v1(1) * v2(1) * v1(2) + std::pow(v2(1), 2) * v0(2) + 2 * v0(1) * v2(1) * v2(2) + std::pow(v2(1), 2) * v1(2) + 2 * v1(1) * v2(1) * v2(2) + v0(1) * v1(1) * v2(2) + v0(1) * v2(1) * v1(2) + v1(1) * v2(1) * v0(2)));
        Izx += density * n(2) * 0.5 * (inv20 * (std::pow(v0(2), 2) * v0(0) + std::pow(v1(2), 2) * v1(0) + std::pow(v2(2), 2) * v2(0)) + inv60 * (std::pow(v0(2), 2) * v1(0) + 2 * v0(2) * v1(2) * v0(0) + std::pow(v0(2), 2) * v2(0) + 2 * v0(2) * v2(2) * v0(0) + std::pow(v1(2), 2) * v0(0) + 2 * v0(2) * v1(2) * v1(0) + std::pow(v1(2), 2) * v2(0) + 2 * v1(2) * v2(2) * v1(0) + std::pow(v2(2), 2) * v0(0) + 2 * v0(2) * v2(2) * v2(0) + std::pow(v2(2), 2) * v1(0) + 2 * v1(2) * v2(2) * v2(0) + v0(2) * v1(2) * v2(0) + v0(2) * v2(2) * v1(0) + v1(2) * v2(2) * v0(0)));
    }
    Ixx -= mass * (std::pow(center(1), 2) + std::pow(center(2), 2));
    Iyy -= mass * (std::pow(center(0), 2) + std::pow(center(2), 2));
    Izz -= mass * (std::pow(center(0), 2) + std::pow(center(1), 2));
    Ixy -= mass * center(0) * center(1);
    Iyz -= mass * center(1) * center(2);
    Izx -= mass * center(2) * center(0);
    I << Ixx, -Ixy, -Izx, -Ixy, Iyy, -Iyz, -Izx, -Iyz, Izz;
}