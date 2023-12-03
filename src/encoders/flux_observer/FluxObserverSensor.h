#ifndef FLUX_OBSERVER_SENSOR_H
#define FLUX_OBSERVER_SENSOR_H

#include "Arduino.h"
#include "common/base_classes/FOCMotor.h"
#include "common/base_classes/Sensor.h"
#include "common/multi_filter.h"
/**
  
*/

class FluxObserverSensor : public Sensor
{
  public:
    /**
    FluxObserverSensor class constructor
    @param m  Motor that the FluxObserverSensor will be linked to
    */
    FluxObserverSensor(const FOCMotor& m);
    void update() override;
 
    void init() override;

    // Abstract functions of the Sensor class implementation
    /** get current angle (rad) */
    float getSensorAngle() override;

    
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
    float i_alpha, i_beta, i_ah, i_bh, i_ah_prev, i_bh_prev; //Stores the band passed currents and previous difference values
    float Ts, e, e_in_prev, theta_in, theta_out, theta_out_prev, wrotor, wrotor_prev, kp, ki, ke; //PLL values
    MultiFilter filter_calc_a, filter_calc_b, a_lpf, b_lpf, e_lpf, db_lpf; //Filters for HFI
    bool hfi_calculated=false;
    bool hfi_converged=false;
    float convergence_threshold;
    float grad_db, prev_db;
    float db;
    long pll_samp_time_prev;
    PhaseCurrent_s current;

  protected:    
    const FOCMotor& _motor;

};

#endif
