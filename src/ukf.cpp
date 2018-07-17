#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

//https://github.com/ksakmann/Unscented-Kalman-Filter/blob/master/visualization.ipynb

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  outputLaser.open("../output/nisLaser.log", ios::out);
  if (!outputLaser.is_open())
  {
    cout << "Error opening nisLaser.log" << endl;
    exit(1);
  }

  outputRadar.open("../output/nisRadar.log", ios::out);
  if (!outputRadar.is_open())
  {
    cout << "Error opening nisRadar.log" << endl;
    exit(1);
  }

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.8;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;
  previous_timestamp_ = 0;
  n_x_ = 5;
  n_aug_ = n_x_ + 2;
  lambda_ = 3 - n_aug_;

  // initial weights vector
  weights_ = VectorXd(2 * n_aug_ + 1);

  // initial predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // initial aug state vector
  x_aug_ = VectorXd(n_aug_);

  // initial aug state covariance matrix
  P_aug_ = MatrixXd(n_aug_, n_aug_);

  // initial sigma point matrix
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // predicted measurement covariance matrices
  S_laser_ = MatrixXd(2,2);
  S_radar_ = MatrixXd(3,3);

  // Noise matrices
  R_laser_ = MatrixXd(2,2);
  R_radar_ = MatrixXd(3,3);

  // initial NIS
  NIS_radar_ = 0.0;
  NIS_laser_ = 0.0;
}

UKF::~UKF() {
  outputLaser.close();
  outputRadar.close();
}

void UKF::InitWeights(void) {
  weights_.fill(0.5 / (n_aug_ + lambda_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);
}

void UKF::InitAugMeanState(void) {
  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;
}

void UKF::InitAugCovarMatrix(void) {
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(n_x_, n_x_) = P_;
  P_aug_(5,5) = std_a_ * std_a_;
  P_aug_(6,6) = std_yawdd_ * std_yawdd_;
}

void UKF::InitSensorsNoiseCovarMatrix(void) {
  // measurement noise covariance matrix for laser
  R_laser_ << std_laspx_ * std_laspx_, 0,
              0, std_laspy_ * std_laspy_;

  // measurement noise covariance matrix for radar
  R_radar_ << std_radr_ * std_radr_, 0, 0,
              0, std_radphi_ * std_radphi_, 0,
              0, 0,std_radrd_ * std_radrd_;
}

void UKF::NormalizeAngle(VectorXd &in, int index) {
  while (in(index) > M_PI)
    in(index) -= 2.0 * M_PI;
  while (in(index) < -M_PI)
    in(index) += 2.0 * M_PI;
}
/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
  // initialize state vector
    x_ << 0, 0, 0, 0, 0;

    // initialize weights
    InitWeights();

    // initialize augmented state covariance matrix
    InitAugCovarMatrix();

    // initialize measurement noise covariance matrix for laser and radar sensors
    InitSensorsNoiseCovarMatrix();

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // intialize state covariance matrix
      double val = std_radr_ * std_radr_;
      P_ << val, 0, 0, 0, 0,
            0, val, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, std_radphi_, 0,
            0, 0, 0, 0, std_radphi_;

      // convert radar from polar to cartesian coordinates and initialize state.
      float rho = meas_package.raw_measurements_[0];     // range: radial distance from origin
      float phi = meas_package.raw_measurements_[1];     // bearing: angle between rho and x axis
      float rho_dot = meas_package.raw_measurements_[2]; // radial velocity: change of rho
      // Convert from polar to cartesian coordinates
      float x = rho * cos(phi);
      float y = rho * sin(phi);
      float vx = rho_dot * cos(phi);
      float vy = rho_dot * sin(phi);
      x_ << x, y, sqrt(vx * vx + vy * vy), 0, 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // intialize state covariance matrix
      double val = std_laspx_ * std_laspx_;
      P_ << val, 0, 0, 0, 0,
            0, val, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;

      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    // check for zero
    if (fabs(x_(0)) < 0.001)
      x_(0) = 0.001;
    if (fabs(x_(1)) < 0.001)
      x_(1) = 0.001;

    // Update augmented state with updated state vector
    InitAugMeanState();

    // Save the timestamp for dt calculation
    previous_timestamp_ = meas_package.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  // compute the time elapsed between the current and previous measurements
  float delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;

  while (delta_t > 0.1) {
    const double dt = 0.05;
    Prediction(dt);
    delta_t -= dt;
  }

  Prediction(delta_t);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    UpdateRadar(meas_package);
  else
    UpdateLidar(meas_package);
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  // initialize augmented mean state
  InitAugMeanState();

  // initialize augmented covariance matrix
  InitAugCovarMatrix();

  // create square root matrix
  MatrixXd L = P_aug_.llt().matrixL();

  //create augmented sigma points
  Xsig_aug_.col(0) = x_aug_;
  for (int i = 0; i < n_aug_; i++)
  {
    MatrixXd val = sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug_.col(i+1)        = x_aug_ + val;
    Xsig_aug_.col(i+1+n_aug_) = x_aug_ - val;
  }

  // predict sigma points
  for (int i = 0; i< 2 * n_aug_ + 1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * ( sin (yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw + yawd * delta_t) );
    }
    else {
      px_p = p_x + v*delta_t * cos(yaw);
      py_p = p_y + v*delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    // add noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    // write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // iterate over sigma points
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // iterate over sigma points
  // state difference
  VectorXd x_diff = Xsig_pred_.col(i) - x_;
  // angle normalization
  NormalizeAngle(x_diff, 3);

  P_ += weights_(i) * x_diff * x_diff.transpose() ;
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  // set measurement dimension, px, py
  int n_z_ = 2;

  // set the measurements
  VectorXd z_ = meas_package.raw_measurements_;

  // create matrix for sigma points in measurement space
  MatrixXd Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // 2n+1 simga points
    // sigma point predictions in process space
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);

    // sigma point predictions in measurement space
    Zsig_(0,i) = px;
    Zsig_(1,i) = py;
  }

  // initialize measurement prediction vector
  VectorXd z_pred_ = VectorXd(n_z_);
  z_pred_.fill(0.0);

  //mean predicted measurement
  z_pred_ = Zsig_ * weights_;

  // initialize matrix for predicted measurement covariance
  S_laser_.fill(0.0);

  // measurement covariance matrix S
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;
    S_laser_ += weights_(i) * z_diff * z_diff.transpose();
  }

  S_laser_ += R_laser_;

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  Tc.fill(0.0);

  //calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    // residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S_laser_.inverse();

  //residual
  VectorXd z_diff = z_ - z_pred_;

  // update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K * S_laser_ * K.transpose();

  // NIS Lidar Update
  NIS_laser_ = z_diff.transpose() * S_laser_.inverse() * z_diff;
  outputLaser  << NIS_laser_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  // set measurement dimension, radar can measure r, phi, and r_dot
  int n_z_ = 3;

  // set the measurements
  VectorXd z_ = meas_package.raw_measurements_;

  // create matrix for sigma points in measurement space
  MatrixXd Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // 2n+1 simga points
    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    // check for zero
    if (fabs(p_x) < 0.0001)
      p_x = 0.0001;
    if (fabs(p_y) < 0.0001)
      p_y = 0.0001;

    // measurement model
    Zsig_(0,i) = sqrt(p_x * p_x + p_y * p_y);                         //rho
    Zsig_(1,i) = atan2(p_y, p_x);                                     //phi
    Zsig_(2,i) = (p_x * v1 + p_y * v2 ) / sqrt(p_x * p_x + p_y * p_y);  //rho dot
  }

  //initialize measurement prediction vector
  VectorXd z_pred_ = VectorXd(n_z_);
  z_pred_.fill(0.0);

  // mean predicted measurement
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      z_pred_ = weights_(i) * Zsig_.col(i);
  }

  // initialize matrix for predicted measurement covariance
  S_radar_.fill(0.0);

  // measurement covariance matrix S
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;

    // angle normalization
    NormalizeAngle(z_diff, 1);

    S_radar_ += weights_(i) * z_diff * z_diff.transpose();
  }

  S_radar_ += R_radar_;

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  Tc.fill(0.0);

  //calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;
    // angle normalization
    NormalizeAngle(z_diff, 1);

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // angle normalization
    NormalizeAngle(x_diff, 3);

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S_radar_.inverse();

  // residual
  VectorXd z_diff = z_ - z_pred_;

  // angle normalization
  NormalizeAngle(z_diff, 1);

  // update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K * S_radar_ * K.transpose();

  // NIS Update
  NIS_radar_ = z_diff.transpose() * S_radar_.inverse() * z_diff;
  outputRadar << NIS_radar_ << endl;
}
