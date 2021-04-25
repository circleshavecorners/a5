#include <rigid_body_jacobian.h>

void rigid_body_jacobian(Eigen::Matrix36d &J,
                         Eigen::Ref<const Eigen::Matrix3d> R, Eigen::Ref<const Eigen::Vector3d> p,
                         Eigen::Ref<const Eigen::Vector3d> x)
{
    J.setZero();
    Eigen::Vector3d X = x - p;
    Eigen::Matrix3d crossMatrix;
    crossMatrix << 0, -X.z(), X.y(), X.z(), 0, -X.x(), -X.y(), X.x(), 0;
    J.block(0, 0, 3, 3) = R * crossMatrix.transpose() * R.transpose();
    J.block(0, 3, 3, 3) = Eigen::Matrix3d::Identity();
}
