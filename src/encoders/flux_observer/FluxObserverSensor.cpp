#include "FluxObserverSensor.h"
#include "common/foc_utils.h"
#include "common/time_utils.h"
#include "common/multi_filter.h"

FluxObserverSensor::FluxObserverSensor(const FOCMotor& m) : _motor(m)
{
  // Derive Flux linkage from KV_rating and pole_pairs
  if (_isset(_motor.pole_pairs) && _isset(_motor.KV_rating)){
    flux_linkage = 60 / ( _sqrt(3) * _PI * (_motor.KV_rating) * (_motor.pole_pairs * 2));
  }
  filter_calc_a = MultiFilter(1.0f/1500.0f);
  filter_calc_b = MultiFilter(1.0f/1500.0f);
  a_lpf=MultiFilter(1.0f/(200.0f));
  b_lpf=MultiFilter(1.0f/(200.0f));
  e_lpf=MultiFilter(1.0f/(200.0f));
  theta_in_lpf=MultiFilter(1.0f/(200.0f));
  db_lpf=MultiFilter(1.0f/(200.0f));
  grad_db_lpf=MultiFilter(1.0f/(200.0f));
  e_in_prev=0; //n-1 e into PLL
  theta_out=0;
  theta_out_prev=0;
  wrotor=0; //PLL speed output
  wrotor_prev=0; //n-1 PLL speed output
  kp=0.001/(0.5/1500);//PI value set based on desired dampening/settling time
  ki=0.01/(0.5/1500);//PI value set based on desired dampening/settling time
  ke=0.3;
  convergence_threshold=0.01;
  pll_samp_time_prev=micros();
  prev_db=0;
  grad_db=-0.1;
  e=0;
}


void FluxObserverSensor::update() {
  // Current sense is required for the observer
  if (!_motor.current_sense) return;
  
  // Exit if one of the parameter needed for the flux observer is 0
  if ((_motor.phase_inductance == 0) ||
      (_motor.phase_resistance == 0) ||
      (flux_linkage == 0)) return;

  // Update sensor, with optional downsampling of update rate
  if (sensor_cnt++ < sensor_downsample) return;

  // Close to zero speed the flux observer can resonate
  // Estimate the BEMF and use HFI if it's below the threshold and HFI is enabled
  kp=0.01;//0.1/(0.5/_motor.hfi_frequency);//PI value set based on desired dampening/settling time
  ki=0.01;//0.1/(0.5/_motor.hfi_frequency);//PI value set based on desired dampening/settling time
  float bemf = _motor.voltage.q - _motor.phase_resistance * _motor.current.q;
  if (abs(bemf < bemf_threshold)){
    if(_motor.hfi_enabled){
        if(!_motor.hfi_injection_started){
          return;
        }
        sensor_cnt = 0;
        // read current phase currents
        current = _motor.current_sense->getPhaseCurrents();

        // calculate clarke transform
        if(!current.c){
            // if only two measured currents
            i_alpha = current.a;  
            i_beta = _1_SQRT3 * current.a + _2_SQRT3 * current.b;
        }else if(!current.a){
            // if only two measured currents
            float a = -current.c - current.b;
            i_alpha = a;  
            i_beta = _1_SQRT3 * a + _2_SQRT3 * current.b;
        }else if(!current.b){
            // if only two measured currents
            float b = -current.a - current.c;
            i_alpha = current.a;  
            i_beta = _1_SQRT3 * current.a + _2_SQRT3 * b;
        } else {
            // signal filtering using identity a + b + c = 0. Assumes measurement error is normally distributed.
            float mid = (1.f/3) * (current.a + current.b + current.c);
            float a = current.a - mid;
            float b = current.b - mid;
            i_alpha = a;
            i_beta = _1_SQRT3 * a + _2_SQRT3 * b;
        }
        
        i_ah=filter_calc_a.getBp(i_alpha);
        i_bh=filter_calc_b.getBp(i_beta);
        
        theta_in = theta_in_lpf.getLp(_atan2(b_lpf.getLp(_motor.hfi_state*i_bh),a_lpf.getLp(_motor.hfi_state*i_ah)));
        db=db_lpf.getLp(_motor.hfi_state*((i_bh-i_bh_prev)));
        if(!hfi_converged){
          if(fabs(db)<convergence_threshold && i_ah>0.2){
            hfi_converged=true;
            theta_out=theta_in;
            e=e_lpf.getLp(theta_in-theta_out); //If converged switch over to regular error
            first=0;//Reset angle when converged
           }
          else{
            if(grad_db_lpf.getLp(db-prev_db)>0){
              grad_db*=-1;
            }
            theta_out+=grad_db;
            theta_out=_normalizeAngle(theta_out);
            electrical_angle=theta_out;
            //e=e_lpf.getLp(theta_in-theta_out);
          }
        }
        else{
          e=e_lpf.getLp(theta_in-theta_out); //If converged use regular error
        }
        i_ah_prev=i_ah;
        i_bh_prev=i_bh;
        if(hfi_converged){
        //PLL
        float curr_pll_time=micros();
        Ts=_motor.hfi_dt/1000000.0; //Sample time can be dynamically calculated
        pll_samp_time_prev=curr_pll_time;
        wrotor = ((2*kp+ki*Ts)*e + (ki*Ts-2*kp)*e_in_prev + 2 * (wrotor_prev))/2; //bilinear transform based difference equation of transfer function kp+ki/s
        theta_out = ((Ts/2)*(wrotor+wrotor_prev)+theta_out_prev); //#1/s transfer function. just integration

        //Shift values over
        wrotor_prev=wrotor;
        e_in_prev=e;
        theta_out_prev=theta_out;

        //Set angle
        electrical_angle=theta_out;
        }
        hfi_calculated=true;
    }
    else{
      return;
    }
  }
  float now;
  if(!hfi_calculated){
    sensor_cnt = 0;

  // read current phase currents
    current = _motor.current_sense->getPhaseCurrents();

    // calculate clarke transform
    if(!current.c){
        // if only two measured currents
        i_alpha = current.a;  
        i_beta = _1_SQRT3 * current.a + _2_SQRT3 * current.b;
    }else if(!current.a){
        // if only two measured currents
        float a = -current.c - current.b;
        i_alpha = a;  
        i_beta = _1_SQRT3 * a + _2_SQRT3 * current.b;
    }else if(!current.b){
        // if only two measured currents
        float b = -current.a - current.c;
        i_alpha = current.a;  
        i_beta = _1_SQRT3 * current.a + _2_SQRT3 * b;
    } else {
        // signal filtering using identity a + b + c = 0. Assumes measurement error is normally distributed.
        float mid = (1.f/3) * (current.a + current.b + current.c);
        float a = current.a - mid;
        float b = current.b - mid;
        i_alpha = a;
        i_beta = _1_SQRT3 * a + _2_SQRT3 * b;
    }

      // This work deviates slightly from the BSD 3 clause licence.
      // The work here is entirely original to the MESC FOC project, and not based
      // on any appnotes, or borrowed from another project. This work is free to
      // use, as granted in BSD 3 clause, with the exception that this note must
      // be included in where this code is implemented/modified to use your
      // variable names, structures containing variables or other minor
      // rearrangements in place of the original names I have chosen, and credit
      // to David Molony as the original author must be noted.

      // Flux linkage observer    
      now = _micros();
      float Ts = ( now - angle_prev_ts) * 1e-6f; 
      flux_alpha = _constrain( flux_alpha + (_motor.Ualpha - _motor.phase_resistance * i_alpha) * Ts -
            _motor.phase_inductance * (i_alpha - i_alpha_prev),-flux_linkage, flux_linkage);
      flux_beta  = _constrain( flux_beta  + (_motor.Ubeta  - _motor.phase_resistance * i_beta)  * Ts -
            _motor.phase_inductance * (i_beta  - i_beta_prev) ,-flux_linkage, flux_linkage);
      
      // Calculate angle
        electrical_angle = _normalizeAngle(_atan2(flux_beta,flux_alpha));
  }
  
  // Electrical angle difference
  float d_electrical_angle = 0;
  if (first){
    // Skip angle difference calculation the first time
    first = 0;
    d_electrical_angle = electrical_angle;
  }else{
    d_electrical_angle = electrical_angle - electrical_angle_prev;
    if(abs(d_electrical_angle) > _2PI * 0.8 ){ //change the  factor based on sample rate can also just use _PI for simplicity 
      if (d_electrical_angle > 0){
        d_electrical_angle -= _2PI;
      }else{
        d_electrical_angle += _2PI;
      }
    }
  }
  angle_track += d_electrical_angle;

  // Mechanical angle and full_rotations
  if(abs(angle_track) > _2PI * _motor.pole_pairs){
    if (angle_track>0){
      full_rotations += 1;
      angle_track -= _2PI * _motor.pole_pairs;
    }else{
      full_rotations -= 1;
      angle_track += _2PI * _motor.pole_pairs;
    }
  }
  angle_prev = angle_track /_motor.pole_pairs;
  
  // Store Previous values
  i_alpha_prev = i_alpha;
  i_beta_prev = i_beta;
  angle_prev_ts = now;
  electrical_angle_prev = electrical_angle;

}

void FluxObserverSensor::init(){
  this->Sensor::init(); // call base class
}

/*
	Shaft angle calculation
*/
float FluxObserverSensor::getSensorAngle(){
  return 0;
}