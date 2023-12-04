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
    float i_alpha,i_beta,i_qh, i_qh_prev; //Stores the band passed currents and previous difference values
    float Ts, e, e_in_prev, theta_out, theta_out_prev, wrotor, wrotor_prev, kp, ki, ke; //PLL values
    MultiFilter filter_calc_q,q_lp; //Filters for HFI
    bool hfi_calculated=false;
    PhaseCurrent_s current;


};

#endif
