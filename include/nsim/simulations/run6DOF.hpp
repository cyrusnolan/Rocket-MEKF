#pragma once
#include <Eigen/StdVector>
#include "Eigen/Dense"

#include "gnc/AttitudeController.hpp"
#include "gnc/PositionController.hpp"
#include "gnc/ThrustAllocator.hpp"
#include "nsim/GenerateForcesMoments.hpp"
#include "nsim/eom6DOF.hpp"
#include "nsim/step.hpp"
#include "nsim/SaveData.hpp"

/**
 * @brief Conducts an individual 6DOF simulation
 * 
 * @param z Initial state
 */
void run6DOF(Matrix<double, sim::num_states_6DOF, 1> z);