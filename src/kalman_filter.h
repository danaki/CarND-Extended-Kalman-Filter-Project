#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "Eigen/Dense"
#include "measurement_package.h"

class KalmanFilter {
private:
  // state transistion matrix
  Eigen::MatrixXd F_;

  // process covariance matrix
  Eigen::MatrixXd Q_;

  Eigen::MatrixXd H_;

  // measurement covariance matrix
  Eigen::MatrixXd R_laser_;
  Eigen::MatrixXd R_radar_;
  
  // measurement matrix  
  Eigen::MatrixXd H_laser_;
  Eigen::MatrixXd Hj_;
  
  float noise_ax;
  float noise_ay;  
  
public:
  // state vector
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;
  
  KalmanFilter();

  virtual ~KalmanFilter();

  void Init(const MeasurementPackage &measurement_pack);

  void UpdateTimestamp(float dt);
    
  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   * @param delta_T Time between k and k+1 in s
   */
  void Predict();

  /**
   * Updates the state by using standard Kalman Filter equations
   * @param z The measurement at k+1
   */
  void Update(const Eigen::VectorXd &z);

  /**
   * Updates the state by using Extended Kalman Filter equations
   * @param z The measurement at k+1
   */
  void UpdateEKF(const Eigen::VectorXd &z);

  void PrintState();
    
};

#endif /* KALMAN_FILTER_H_ */
