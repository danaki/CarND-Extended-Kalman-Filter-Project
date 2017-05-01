#include "kalman_filter.h"
#include "measurement_package.h"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {
  x_ = VectorXd(4);  
  
  P_ = MatrixXd(4, 4);
  F_ = MatrixXd(4, 4);
  Q_ = MatrixXd(4, 4);
  
  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);  
}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(const MeasurementPackage &measurement_pack) {
  P_ << 1, 0, 0,    0,
        0, 1, 0,    0,
        0, 0, 1000, 0,
        0, 0, 0,    1000;

  F_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1; 

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
            
  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0,      0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0,      0,
              0,    0.0009, 0,
              0,    0,      0.09;
              
  noise_ax = 9;
  noise_ay = 9;

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Convert radar from polar to cartesian coordinates and initialize state.
    float ro = measurement_pack.raw_measurements_[0];
    float phi = measurement_pack.raw_measurements_[1];
    float ro_dot = measurement_pack.raw_measurements_[2];
    
    x_ << ro * cos(phi), ro * sin(phi), ro_dot * cos(phi), ro_dot * sin(phi);
  }
  else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {    
    // Initialize state.
    x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
  }
}

void KalmanFilter::UpdateTimestamp(float dt) {
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  //Modify the F matrix so that the time is integrated
  F_(0, 2) = dt;
  F_(1, 3) = dt;

  //set the process covariance matrix Q
  Q_ << dt_4 / 4 * noise_ax, 0,                   dt_3 / 2 * noise_ax, 0,
        0,                   dt_4 / 4 * noise_ay, 0,                   dt_3 / 2 * noise_ay,
        dt_3 / 2 * noise_ax, 0,                   dt_2 * noise_ax,     0,
        0,                   dt_3 / 2 * noise_ay, 0,                   dt_2 * noise_ay;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd y = z - H_laser_ * x_;
  MatrixXd Ht = H_laser_.transpose();

  MatrixXd S = H_laser_ * P_ * Ht + R_laser_;  
  
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;
  
  x_ = x_ + (K * y);
  
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_laser_) * P_;    
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  if (z(0) == 0 && z(1) == 0 && z(2) == 0) {
      return;
  }
    
  Hj_ = Tools::CalculateJacobian(x_);    
   
  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];
  
  float h1 = sqrt(px * px + py * py);
  float h2 = atan2(py, px);
  float h3 = (px * vx + py * vy) / h1;  
  
  VectorXd hp = VectorXd(3);  
  hp << h1, h2, h3;
  
  VectorXd y = z - hp;
  MatrixXd Ht = Hj_.transpose();

  MatrixXd S = Hj_ * P_ * Ht + R_radar_;  
  
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;
  x_ = x_ + (K * y);
  
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj_) * P_;    
}
