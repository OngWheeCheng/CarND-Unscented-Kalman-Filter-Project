#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
private:
  void InitWeights(void);
  void InitAugMeanState(void);
  void InitAugCovarMatrix(void);
  void InitSensorsNoiseCovarMatrix(void);
  void NormalizeAngle(VectorXd &in, int index);

public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* aug state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate nu psi] in SI units and rad
  VectorXd x_aug_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* aug state covariance matrix
  MatrixXd P_aug_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* aug sigma points matrix
  MatrixXd Xsig_aug_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* predicted measurement covariance matrix for laser
  MatrixXd S_laser_;

  ///* predicted measurement covariance matrix for radar
  MatrixXd S_radar_;

  ///* measurement noise covariance matrix for laser
  MatrixXd R_laser_;

  ///* measurement noise covariance matrix for radar
  MatrixXd R_radar_;

  ///* NIS for radar
  double NIS_radar_;

  ///* NIS for laser
  double NIS_laser_;

  // previous timestamp
  long long previous_timestamp_;

  // Output filestream for radar and laser NIS
  std::ofstream outputLaser;
  std::ofstream outputRadar;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);
};

#endif /* UKF_H */
