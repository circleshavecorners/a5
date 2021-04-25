#include <rodrigues.h>
#include <cmath>
#include <iostream>

void rodrigues(Eigen::Matrix3d &R, Eigen::Ref<const Eigen::Vector3d> omega, double dt)
{
    Eigen::Matrix3d omega_, omega_2;
    omega_ << 0, -omega.z(), omega.y(), omega.z(), 0, -omega.x(), -omega.y(), omega.x(), 0;
    omega_2 = omega * omega.transpose() - std::pow(omega.norm(), 2) * Eigen::Matrix3d::Identity();
    double omega_norm = omega.norm();

    Eigen::Matrix3d R_ = Eigen::Matrix3d::Identity() + omega_ * std::sin(omega.norm() * dt) / omega.norm() + omega_2 * (1 - std::cos(omega_norm * dt)) / std::pow(omega_norm, 2);
    R = (R_ * R).eval();
}