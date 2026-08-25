#pragma once

#include <array>

#include <Eigen/Core>

#include "common/POD_Event.hpp"

#include "KalmanFitter/KalmanFitterUtils.hxx"

namespace T2DS::KF {

// ## KF::Vertex ## //

struct Vertex {

    // Constructor //

    // Create a `KF::Vertex` from the event's reconstructed primary vertex.
    [[nodiscard]] static Vertex FromEvent(const POD::Event &e) {
        Vertex out;

        out.xyz(0) = static_cast<double>(e.PV_X);
        out.xyz(1) = static_cast<double>(e.PV_Y);
        out.xyz(2) = static_cast<double>(e.PV_Z);

        UnpackSym<3>(out.cov, e.PV_CovMatrix);

        return out;
    }

    // Getter //

    [[nodiscard]] std::array<double, 3> GetXYZ() const { return {xyz(0), xyz(1), xyz(2)}; }

    // Member Variables //

    Eigen::Matrix<double, 3, 3> cov{Eigen::Matrix<double, 3, 3>::Zero()};  // full symmetric
    Eigen::Vector<double, 3> xyz{Eigen::Vector<double, 3>::Zero()};
};

}  // namespace T2DS::KF
