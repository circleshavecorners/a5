#include <Eigen/Dense>
#include <EigenTypes.h>
#include <rodrigues.h>
#include <iostream>

//Input:
//  q - 12n vector where n is the number of rigid bodies. Each rigid body is stored as 12 doubles.
//      The first 9 doubles are the columns of a 3x3 rotation matrix and the final 3 doubles are the world space position of the object's center of mass.
//  qdot - 6n vector of generalied velocities. The first 3 doubles of each body are the world space angular velocity and
//         the second 3 are the world space linear velocity.
//  dt - the integration time step
//  masses - a vector to mass matrices for each rigid body
//  forces - a 6n vector of generalized forces for n rigid bodies. The first 3 doubles of each rigid body are the torques acting on the object
//           while the second 3 doubles are the linear forces.
//Output:
//  q - updated generalized coordinates
//  qdot - updated generalized velocities
inline void exponential_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt,
                              std::vector<Eigen::Matrix66d> &masses, Eigen::Ref<const Eigen::VectorXd> forces)
{
    for (unsigned int irb = 0; irb < q.rows() / 12; ++irb)
    {
        Eigen::Matrix3d R = Eigen::Map<const Eigen::Matrix3d>(q.segment<9>(12 * irb).data());
        Eigen::Vector3d p = Eigen::Map<const Eigen::Vector3d>(q.segment<3>(12 * irb + 9).data());
        Eigen::Matrix3d Interia = masses[irb].block(0, 0, 3, 3);
        Eigen::Matrix3d RIRT = R * Interia * R.transpose();
        Eigen::Vector3d omega = qdot.segment<3>(6 * irb);
        qdot.segment<3>(6 * irb) = (RIRT).inverse() * (RIRT * omega + dt * omega.cross(RIRT * omega) + dt * forces.segment<3>(6 * irb));
        qdot.segment<3>(6 * irb + 3) += dt * forces.segment<3>(6 * irb + 3) / masses[irb](3, 3);
        rodrigues(R, qdot.segment<3>(6 * irb), dt);
        p += dt * qdot.segment<3>(6 * irb + 3);
        q.segment<9>(12 * irb) = Eigen::Map<const Eigen::Vector9d>(R.data());
        q.segment<3>(12 * irb + 9) = p;
    }
}