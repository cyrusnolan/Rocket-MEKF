#include "gnc/AttitudeController.hpp"
#include "gnc/PositionController.hpp"
#include "gnc/constants.hpp"
#include "gnc/utilities.hpp"
#include "nsim/SaveData.hpp"
#include "nsim/constants.hpp"
#include "nsim/eom.hpp"
#include "nsim/step.hpp"
#include <Eigen/StdVector>
#include <chrono>
#include <iostream>

int main() {

  // Initial Conditions
  Matrix<double, 13, 1> z;
  Vector3d M_des;
  Vector3d F_des;
  Vector3d q_int;
  Vector3d pos_int;
  std::vector<Matrix<double, sim::num_states_6DOF, 1>> data;

  q_int << 0, 0, 0;
  pos_int << 0, 0, 0;

  z << 10, 20, 300, 1, 0, 0, 0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;

  auto start = std::chrono::high_resolution_clock::now();
  while (z(2, 0) >= 0) {
    M_des =
        AttitudeController(z.block<4, 1>(3, 0), z.block<3, 1>(10, 0), q_int);
    F_des =
        PositionController(z.block<3, 1>(0, 0), z.block<3, 1>(7, 0), pos_int);
    z = step(&eom, z, sim::h, F_des, M_des);
    // Stores state to data vector
    data.push_back(z);
  }
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout << duration.count() << std::endl;
  SaveData(data);
}