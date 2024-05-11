#include <iomanip>  // Include this at the top of your file
#include "nsim/SaveData.hpp"
#include "gnc/constants.hpp"

void SaveData(const std::vector<Eigen::Matrix<double, 15, 1>> &data) {
  std::ofstream output_file("../analysis/fc_data.csv");
  output_file << std::setprecision(15);
  int j = 1;
  std::for_each(data.begin(), data.end(), [&](auto const &x) {
    output_file << FlightComputer::dt * j;
    output_file << ",";
    for (int i = 0; i < 15; i++) {
      output_file << x(i, 0);
      output_file << ",";
      if (i == 14) {
        output_file << std::endl;
      }
    }
    j++;
  });
}

void SaveData(const std::vector<Eigen::Matrix<double, 20, 1>> &data) {
  std::ofstream output_file("../analysis/sim_data.csv");
  output_file << std::setprecision(15);
  int j = 1;
  std::for_each(data.begin(), data.end(), [&](auto const &x) {
    output_file << sim::h * j;
    output_file << ",";
    for (int i = 0; i < 20; i++) {
      output_file << x(i, 0);
      output_file << ",";
      if (i == 19) {
        output_file << std::endl;
      }
    }
    j++;
  });
}
