#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
   * End DO NOT MODIFY section for measurement noise values
   */

  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

  // state dimension
  int n_x_ = 5;

  // augmented state dimension
  int n_aug_ = 7;

  // spreading parameter for sigma points generation
  double lambda_ = 3 - n_aug_;

  // Initialise velocity, yaw angle, and yaw rate values for state vector x
  x_(2) = 0;
  x_(3) = 0;
  x_(4) = 0;

  // Initialise state covariance matrix P
  // TODO: change later to s.d. of measurement
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
   switch (meas_package.sensor_type_)
   {
      case MeasurementPackage::LASER:
        if (!is_initialized_)
        {
            x_.head(2) = meas_package.raw_measurements_;
            is_initialized_ = true;
        }
        break;
      case MeasurementPackage::RADAR:
        if (!is_initialized_)
        {
            double rho = meas_package.raw_measurements_[0];
            double phi = meas_package.raw_measurements_[1];
            x_(0) = rho*cos(phi);
            x_(1) = rho*sin(phi);
            is_initialized_ = true;
        }
        break;
   }
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location.
   * Modify the state vector, x_. Predict sigma points, the state,
   * and the state covariance matrix.
   */

  /**
   1. generate augmented sigma points
   */

  // create augmented state vector
  VectorXd x_aug = VectorXd(7);
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  // create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; ++i)
  {
      Xsig_aug.col(i+1)         = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
      Xsig_aug.col(i+1+n_aug_)  = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  /**
   2. predict sigma points
   */

  // create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // predict sigma points
  for (int i = 0; i < 2*n_aug_; ++i)
  {
      double p_x      = Xsig_aug(0,i);
      double p_y      = Xsig_aug(1,i);
      double v        = Xsig_aug(2,i);
      double yaw      = Xsig_aug(3,i);
      double yawd     = Xsig_aug(4,i);
      double nu_a     = Xsig_aug(5,i);
      double nu_yawdd = Xsig_aug(6,i);

      // predict state values
      double px_p, py_p;

      // avoid division by zero
      if (fabs(yawd) > 0.001)
      {
          px_p = p_x + v/yawd * ( sin(yaw + yawd*delta_t) - sin(yaw) );
          py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw + yawd * delta_t) );
      }
      else
      {
          px_p = p_x + v * delta_t * cos(yaw);
          py_p = p_y + v * delta_t * sin(yaw);
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

  /**
   3. predict state and state covariance
   */

  // create vector for weights
  VectorXd weights_ = VectorXd(2 * n_aug_ + 1);
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i = 1; i < 2 * n_aug_ + 1; ++i)
  {
      double weight = 0.5 / (n_aug_ + lambda_);
      weights_(i) = weight;
  }

  // predict state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
      x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  // predict state covariance
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
      // state difference
      VectorXd x_diff = Xsig_pred_.col(i) - x;
      // angle normalisation
      while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
      while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

      P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
  
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}
