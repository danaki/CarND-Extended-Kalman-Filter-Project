#include "FusionEKF.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;


FusionEKF::FusionEKF() {
  is_initialized_ = false;
  previous_timestamp_ = 0;
}

FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (! is_initialized_) {
    // first measurement
    cout << "EKF: " << endl;   
    
    ekf_.Init(measurement_pack);
    
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;    //dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.UpdateTimestamp(dt);

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);    
  } else {
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
