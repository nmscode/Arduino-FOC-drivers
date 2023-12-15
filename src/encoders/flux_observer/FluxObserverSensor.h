#ifndef FLUX_OBSERVER_SENSOR_H
#define FLUX_OBSERVER_SENSOR_H

#include "Arduino.h"
#include "common/base_classes/FOCMotor.h"
#include "common/base_classes/Sensor.h"
#include "common/multi_filter.h"
#include "BLDCMotor.h"
/**
  
*/

class FluxObserverSensor : public Sensor
{
  public:
    /**
    FluxObserverSensor class constructor
    @param m  Motor that the FluxObserverSensor will be linked to
    */
    FluxObserverSensor(BLDCMotor* m);
    void update() override;
 
    void init() override;

    // Abstract functions of the Sensor class implementation
    /** get current angle (rad) */
    float getSensorAngle() override;

    BLDCMotor* _motor;
    // For sensors with slow communication, use these to poll less often
    unsigned int sensor_downsample = 0; // parameter defining the ratio of downsampling for sensor update
    unsigned int sensor_cnt = 0; // counting variable for downsampling
    float flux_alpha = 0; // Flux Alpha 
    float flux_beta = 0; // Flux Beta
    float flux_linkage = 0; // Flux linkage, calculated based on KV and pole number
    float i_alpha_prev = 0; // Previous Alpha current
    float i_beta_prev = 0; // Previous Beta current
    float electrical_angle_prev = 0; // Previous electrical angle
    float electrical_angle;
    float angle_track = 0; // Total Electrical angle
    float bemf_threshold = 1000; // Bemf voltage amplitude when the flux observer should start tracking
    int8_t first = 1; // To skip angle difference calculation the first time
    float i_alpha, i_beta, i_qh, i_qh_prev, i_dh, i_dh_prev; //Stores the band passed currents and previous difference values
    float Ts, e, theta_out, wrotor,accel, second_integral_input, second_integral_input_prev, third_integral_input, third_integral_input_prev, input, input_prev, ka, kw, ktheta, sigma; //Position Observer values
    MultiFilter filter_calc_q, q_lp, filter_calc_d, d_lp, theta_lpf_sin,theta_lpf_cos; //Filters for HFI
    bool hfi_calculated=false;
    PhaseCurrent_s current;
    float smooth_theta_cos, smooth_theta_sin;
    float delta_i_qh, delta_i_dh ,atan_test;
    float prev_pll_time;
    float hfi_frequency=2000;
    float flux_observer_theta_out;
  
};

#endif
