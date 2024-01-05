#include "FluxObserverSensor.h"
#include "common/foc_utils.h"
#include "common/time_utils.h"
#include "common/multi_filter.h"
#include "BLDCmotor.h"

FluxObserverSensor::FluxObserverSensor(BLDCMotor* m)
{
  // Derive Flux linkage from KV_rating and pole_pairs
  if (_isset(_motor->pole_pairs) && _isset(_motor->KV_rating)){
    flux_linkage = 60 / ( _sqrt(3) * _PI * (_motor->KV_rating) * (_motor->pole_pairs * 2));
  }
  float hpf_freq=1600.0f;
  float q_val=0.7;
  q_hp = MultiFilter(1.0f/hpf_freq,q_val);
  q_hp2 = MultiFilter(1.0f/hpf_freq,q_val);
  q_hp3 = MultiFilter(1.0f/hpf_freq,q_val);
  q_hp4 = MultiFilter(1.0f/hpf_freq,q_val);

  float lpf_freq=20.0f;
  q_lp=MultiFilter(1.0f/lpf_freq,q_val);
  q_lp2=MultiFilter(1.0f/lpf_freq,q_val);
  q_lp3=MultiFilter(1.0f/lpf_freq,q_val);
  q_lp4=MultiFilter(1.0f/lpf_freq,q_val);

 
  theta_lpf_sin=MultiFilter(1.0f/100.0f);
  theta_lpf_cos=MultiFilter(1.0f/100.0f);

  // filter_calc_d = MultiFilter(1.0f/1500.0f);
  //d_lp=MultiFilter(1.0f/200.0f);

  theta_out=0.0f;
  _motor=m;
  accel=0.0f;
  theta_out_prev=0.0f;
  wrotor=0.0f;
  wrotor_prev=0.0f;
  input=0.0f;
  input_prev=0.0f;
  second_integral_input=0.0f;
  second_integral_input_prev=0.0f;
  prev_pll_time=micros();
  sigma=0.0f;
  ka=0.0f;
  kw=10.0f;
  ktheta=0.0f;
  electrical_angle=0;
}


void FluxObserverSensor::update() {
  // Current sense is required for the observer
  if (!_motor->current_sense) return;
  
  // Exit if one of the parameter needed for the flux observer is 0
  if ((_motor->phase_inductance == 0) ||
      (_motor->phase_resistance == 0) ||
      (flux_linkage == 0)) return;

  // Update sensor, with optional downsampling of update rate
  if (sensor_cnt++ < sensor_downsample) return;

  // Close to zero speed the flux observer can resonate
  // Estimate the BEMF and use HFI if it's below the threshold and HFI is enabled
  //kp=1.0f;//0.1/(0.5/_motor->hfi_frequency);//PI value set based on desired dampening/settling time
  //ki=10.0f;//0.1/(0.5/_motor->hfi_frequency);//PI value set based on desired dampening/settling time

  // read current phase currents
  
  current = current_sense->getPhaseCurrents();//_motor->current;
  unsigned long heterodyne_time=micros();
  unsigned long curr_pll_time=micros();
  Ts=(curr_pll_time-prev_pll_time)/1000000.0f; //Sample time can be dynamically calculated
  prev_pll_time=curr_pll_time;
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

  float bemf = _motor->voltage.q - _motor->phase_resistance * _motor->current.q;
  if (fabs(bemf < bemf_threshold)){
    if(_motor->hfi_enabled){
        if(!_motor->hfi_injection_started){
          prev_pll_time=micros();
          return;
        }
        sensor_cnt = 0;
        
        float ct;
        float st;
        _sincos(_normalizeAngle(electrical_angle-_PI/4.0f), &st, &ct);
        i_qh=(i_beta * ct - i_alpha * st);//q_hp4.getHp(q_hp3.getHp(q_hp2.getHp(q_hp.getHp((current.q)-_motor->current_sp))));
        i_dh=(i_alpha * ct + i_beta * st);;
        //i_dh=filter_calc_d.getBp(i_alpha * ct + i_beta * st);

        

        delta_i_qh=(i_qh-i_qh_prev);
        delta_i_dh=(i_dh-i_dh_prev);
        
        //atan_test=_atan2(i_qh-i_qh_prev,i_dh-i_dh_prev);
        e=_motor->hfi_state*(delta_i_qh-delta_i_dh) //(q_lp.getLp((i_qh)*_cos(_normalizeAngle(_motor->hfi_dt*_2PI/((1.0f/hfi_frequency)*1000000.0f)))));//ke*delta_i_qh;

        
        //Position Observer
        if(e>0.0f){
          sigma=1.0f;
        }
        else if (e<0.0f){
          sigma=-1.0f;
        }

        else{
          sigma=e;
        }
        accel+=ka*sigma*Ts;
        input=kw*sigma+accel;
        wrotor = _constrain(0.1f*input+wrotor_prev,-10,10);
        second_integral_input=wrotor+ktheta*sigma;
        theta_out = Ts*second_integral_input+theta_out_prev; //#1/s transfer function. just integration
        if(theta_out<-_PI){
          theta_out+=_2PI;
        }
        if(theta_out>_PI){
          theta_out-=_2PI;
        }
        //i_qh_prev=i_qh;
        //i_dh_prev=i_dh;
        //Shift values over
        wrotor_prev=wrotor;
        second_integral_input_prev=second_integral_input;
        input_prev=input;

        theta_out_prev=theta_out;
        //smooth_theta_cos=theta_lpf_cos.getLp(_cos(theta_out));
        //smooth_theta_sin=theta_lpf_sin.getLp(_sin(theta_out));
        //Set angle
        electrical_angle=_normalizeAngle(theta_out);//_normalizeAngle(_atan2(smooth_theta_sin,smooth_theta_cos));
        
        //angle_prev = electrical_angle /_motor->pole_pairs;
        hfi_calculated=true;
        //angle_prev_ts=micros();
        //return;
    }
    else{
    return;
    }
  }
  
  float now=micros();
  if(!hfi_calculated){
    sensor_cnt = 0;
      // This work deviates slightly from the BSD 3 clause licence.
      // The work here is entirely original to the MESC FOC project, and not based
      // on any appnotes, or borrowed from another project. This work is free to
      // use, as granted in BSD 3 clause, with the exception that this note must
      // be included in where this code is implemented/modified to use your
      // variable names, structures containing variables or other minor
      // rearrangements in place of the original names I have chosen, and credit
      // to David Molony as the original author must be noted.

      // Flux linkage observer    
      float Ts = ( now - angle_prev_ts) * 1e-6f; 
      flux_alpha = _constrain( flux_alpha + (_motor->Ualpha - _motor->phase_resistance * i_alpha) * Ts -
            _motor->phase_inductance * (i_alpha - i_alpha_prev),-flux_linkage, flux_linkage);
      flux_beta  = _constrain( flux_beta  + (_motor->Ubeta  - _motor->phase_resistance * i_beta)  * Ts -
            _motor->phase_inductance * (i_beta  - i_beta_prev) ,-flux_linkage, flux_linkage);
      
      // Calculate angle
        electrical_angle = _normalizeAngle(_atan2(flux_beta,flux_alpha));
  }
  
  // Electrical angle difference
  float d_electrical_angle = 0;
  if (first){
    // Skip angle difference calculation the first time
    first = 0;
    d_electrical_angle = electrical_angle;
    angle_track=0;
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
  if(abs(angle_track) > _2PI * _motor->pole_pairs){
    if (angle_track>0){
      full_rotations += 1;
      angle_track -= _2PI * _motor->pole_pairs;
    }else{
      full_rotations -= 1;
      angle_track += _2PI * _motor->pole_pairs;
    }
  }
  angle_prev = angle_track /_motor->pole_pairs;
  
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
  return angle_prev;
}