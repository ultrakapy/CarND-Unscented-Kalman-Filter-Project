#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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
  std_a_ = 1;//30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI/4;//30;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  // set state dimension
  n_x_ = 5;

  //set augmented dimension
  n_aug_ = 7;

  //define spreading parameter
  lambda_ = 3 - n_aug_;
  
  //set vector for weights
  weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  x_ = VectorXd(n_x_);
  x_ << 1, 1, 1, 1, 0.1;
}

UKF::~UKF() {}

const double normalize_angle(const double angle) {
  double normalized_angle = angle;
  const double PI = std::atan(1.0)*4;
    
  while (std::fabs(normalized_angle) > PI) {
    if (normalized_angle > PI) {
      normalized_angle -= 2*PI;
    } else {
      normalized_angle += 2*PI;
    }
  }

  return normalized_angle;
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
    /**
    TODO:
      * Initialize the state ukf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "UKF: " << endl;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float ro = meas_package.raw_measurements_[0];
      float theta = meas_package.raw_measurements_[1];
      float normalized_theta = normalize_angle(theta);
      float px = ro*std::cos(normalized_theta);
      float py = ro*std::sin(normalized_theta);

      x_(0) = px;
      x_(1) = py;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_(0) = meas_package.raw_measurements_[0];
      x_(1) = meas_package.raw_measurements_[1];
    }

    //cout << "x_ INIT: " << x_ << "\n\n";
    
    P_ << std_laspx_*std_laspx_, 0, 0, 0, 0,
      0, std_laspy_*std_laspy_, 0, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 1, 0,
      0, 0, 0, 0, 1;
    
    //cout << "P_ INIT: " << P_ << "\n\n";

    previous_timestamp_ = meas_package.timestamp_;

    //use_laser_ = false;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  //compute the time elapsed between the current and previous measurements
  double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;
  //cout << "delta_t: " << dt  << "\n\n";
    
  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    if (use_radar_) {
      UpdateRadar(meas_package);
    }
  } else {
    if (use_laser_) {
      UpdateLidar(meas_package);
    }
  }

  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
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

  /* 1. Generate Sigma Points */
  
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  
  //create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  
  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  
  MatrixXd Q = MatrixXd(2, 2);
  Q.fill(0.0);
  Q(0,0) = std_a_*std_a_;
  Q(1,1) = std_yawdd_*std_yawdd_;
  
  P_aug.bottomRightCorner(2,2) = Q;
  
  //create square root matrix
  MatrixXd A_aug = P_aug.llt().matrixL();
  
  //create augmented sigma points
  double lambda_plus_naug = lambda_ + n_aug_;
  double sqrt_lambda_plus_naug = sqrt(lambda_plus_naug);
  MatrixXd sqrt_lambda_plus_naug_times_A = sqrt_lambda_plus_naug * A_aug;
  
  Xsig_aug.col(0) = x_aug;
  
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i+1) = x_aug + sqrt_lambda_plus_naug_times_A.col(i);
    Xsig_aug.col(i+n_aug_+1) = x_aug - sqrt_lambda_plus_naug_times_A.col(i);
  }
  
  /* 2. Predict Sigma Points */

  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
  Xsig_pred_.fill(0.0);

  double px, py, v, psi, psi_dot;
  double px_k, py_k, v_k, psi_k, psi_dot_k, v_a_noise_k, v_pdd_noise_k;
  double delta_t_sq = delta_t*delta_t;

  for (int i = 0; i < (2*n_aug_)+1; i++) {
      px_k = Xsig_aug(0,i);
      py_k = Xsig_aug(1,i);
      v_k = Xsig_aug(2,i);
      psi_k = Xsig_aug(3,i);
      psi_dot_k = Xsig_aug(4,i);
      v_a_noise_k = Xsig_aug(5,i);
      v_pdd_noise_k = Xsig_aug(6,i);

      if (fabs(psi_dot_k) > 0.001) {
        px = px_k + (v_k/psi_dot_k)*(sin(psi_k+(psi_dot_k*delta_t))-sin(psi_k)) + 0.5*delta_t_sq*cos(psi_k)*v_a_noise_k;
        py = py_k + (v_k/psi_dot_k)*(-cos(psi_k+(psi_dot_k*delta_t))+cos(psi_k)) + 0.5*delta_t_sq*sin(psi_k)*v_a_noise_k;
      } else {
        px = px_k + v_k*cos(psi_k)*delta_t + 0.5*delta_t_sq*cos(psi_k)*v_a_noise_k;
        py = py_k + v_k*sin(psi_k)*delta_t + 0.5*delta_t_sq*sin(psi_k)*v_a_noise_k;
      }
      
      v = v_k + delta_t*v_a_noise_k;
      psi = psi_k + psi_dot_k*delta_t + 0.5*delta_t_sq*v_pdd_noise_k;
      psi_dot = psi_dot_k + delta_t*v_pdd_noise_k;
      
      //Xsig_pred 
      Xsig_pred_.col(i) << px, py, v, psi, psi_dot;
  }

  /* 3. Predict Mean and Covariance */

  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);
  x.fill(0.0);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);
  P.fill(0.0);
  
  //predict state mean
  for (int i = 0; i < 2*n_aug_+1; i++) {
    x += weights_(i)*Xsig_pred_.col(i);
  }
  
  //predict state covariance matrix
  for (int i = 0; i < 2*n_aug_+1; i++) {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    x_diff(3) = normalize_angle(x_diff(3));
    
    P += weights_(i)*x_diff*x_diff.transpose();
  }

  x_ = x;
  P_ = P;
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
  // set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  
  //transform sigma points into measurement space
  for (int i = 0; i < 2*n_aug_+1; i++) {
    double px  = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    Zsig(0,i) = px;
    Zsig(1,i) = py;
  }
  
  //calculate mean predicted measurement
  for (int i = 0; i < 2*n_aug_+1; i++) {
    z_pred += weights_(i)*Zsig.col(i);
  }
  
  //calculate innovation covariance matrix S
  for (int i = 0; i < 2*n_aug_+1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    
    //angle normalization
    z_diff(1) = normalize_angle(z_diff(1));
    
    S += weights_(i)*z_diff*z_diff.transpose();
  }
  
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<  std_laspx_*std_laspx_, 0,
        0, std_laspy_*std_laspy_;
  S += R;
  
  /* Update State */
  //create vector for incoming lidar measurement
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0],   //px
    meas_package.raw_measurements_[1];   //py

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  
  //calculate cross correlation matrix
  for (int i = 0; i < 2*n_aug_+1; i++) {
    VectorXd xk_diff = Xsig_pred_.col(i) - x_;
    xk_diff(3) = normalize_angle(xk_diff(3));
    
    VectorXd zk_diff = Zsig.col(i) - z_pred;
    zk_diff(1) = normalize_angle(zk_diff(1));
    
    Tc += weights_(i)*xk_diff*zk_diff.transpose();
  }
  
  //calculate Kalman gain K;
  MatrixXd K = Tc*S.inverse();
  
  //update state mean and covariance matrix
  VectorXd z_diff = z - z_pred;
  z_diff(1) = normalize_angle(z_diff(1));
  
  x_ += K*z_diff;
  P_ -= K*S*K.transpose();
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
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  
  //transform sigma points into measurement space
  for (int i = 0; i < 2*n_aug_+1; i++) {
    double px  = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double psi = Xsig_pred_(3,i);

    double rho = sqrt(px*px+py*py);
    double phi = atan2(py,px);
    double rho_dot = (px*cos(psi)*v+py*sin(psi)*v)/rho;
    Zsig(0,i) = rho;
    Zsig(1,i) = phi;
    Zsig(2,i) = rho_dot;
  }
  
  //calculate mean predicted measurement
  for (int i = 0; i < 2*n_aug_+1; i++) {
    z_pred += weights_(i)*Zsig.col(i);
  }
  
  //calculate innovation covariance matrix S
  for (int i = 0; i < 2*n_aug_+1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    
    //angle normalization
    z_diff(1) = normalize_angle(z_diff(1));

    S += weights_(i)*z_diff*z_diff.transpose();
  }
  
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<  std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0, std_radrd_*std_radrd_;
  S += R;
  
  /* Update State */
  //create vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z <<
      meas_package.raw_measurements_[0],   //rho in m
      meas_package.raw_measurements_[1],   //phi in rad
      meas_package.raw_measurements_[2];   //rho_dot in m/s

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  
  //calculate cross correlation matrix
  for (int i = 0; i < 2*n_aug_+1; i++) {
    VectorXd xk_diff = Xsig_pred_.col(i) - x_;
    xk_diff(3) = normalize_angle(xk_diff(3));
    
    VectorXd zk_diff = Zsig.col(i) - z_pred;
    zk_diff(1) = normalize_angle(zk_diff(1));
    
    Tc += weights_(i)*xk_diff*zk_diff.transpose();
  }
  
  //calculate Kalman gain K;
  MatrixXd K = Tc*S.inverse();
  
  //update state mean and covariance matrix
  VectorXd z_diff = z - z_pred;
  z_diff(1) = normalize_angle(z_diff(1));
  
  x_ += K*z_diff;
  P_ -= K*S*K.transpose();
}
