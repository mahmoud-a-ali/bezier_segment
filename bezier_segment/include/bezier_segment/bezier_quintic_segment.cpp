#include "bezier_quintic_segment.h"

bezier_quintic_segment::bezier_quintic_segment()
{

}

void bezier_quintic_segment::init(double start_time, double duration, std::vector<double> start_state, std::vector<double> end_state){
    if (start_time < 0)
      throw(std::invalid_argument("segment can't be constructed: start_time should be positive "));

    if ( duration < 0)
      throw(std::invalid_argument("segment can't be constructed: duration should be positive."));

    if (start_state.empty() || end_state.empty())
      throw(std::invalid_argument("segment can't be constructed: start & end points  can't be empty."));

    if (start_state.size() != 3)
      throw(std::invalid_argument("segment can't be constructed: startpoint should have 3 values (pos, vel, acc)"));

    if (end_state.size() != 3)
      throw(std::invalid_argument("segment can't be constructed: endpoint should have 3 values (pos, vel, acc)"));

    t_start_ = start_time;
    set_duration_( duration);
    start_state_ = start_state;
    end_state_ = end_state;
    coef_.resize(6);
    iter_=0;
}


// set maximum values of position, velocit, acceleration and jerk
void bezier_quintic_segment::set_absolute_limits(double max_pos, double max_vel, double max_acc, double max_jrk){
    if(max_pos < 0 || max_vel < 0 || max_acc < 0|| max_jrk < 0)
        throw(std::invalid_argument("segment can't be constructed: absolute limits should be positive"));
    max_pos_ = max_pos;
    max_vel_ = max_vel;
    max_acc_ = max_acc;
    max_jrk_ = max_jrk;
}


// compute bezir curve coefficients suing the given initial and final states of the segment 
void bezier_quintic_segment::compute_coef(){
    coef_.resize(6);
    double t1=t_start_, t2=t_end_;
    double p0=start_state_[0], v0=start_state_[1], a0=start_state_[2], pn=end_state_[0], vn=end_state_[1], an=end_state_[2];
    coef_[0] = -p0/pow((t1 - t2), 5) ;
    coef_[1] = -(5*p0 - t1*v0 + t2*v0)/(5*pow((t1 - t2),5));
    coef_[2] = -(a0*t1*t1 - 2*a0*t1*t2 - 8*v0*t1 + a0*t2*t2 + 8*v0*t2 + 20*p0)/(20*pow((t1 - t2), 5));
    coef_[3] = -(an*t1*t1 - 2*an*t1*t2 + 8*vn*t1 + an*t2*t2 - 8*vn*t2 + 20*pn)/(20*pow((t1 - t2), 5));
    coef_[4] = -(5*pn + t1*vn - t2*vn)/(5*pow((t1 - t2), 5));
    coef_[5] = -pn/pow((t1 - t2), 5);
}


// find time of inflection points where a curve (vel/acc/jrk) has its maximum and minimum values 
void bezier_quintic_segment::compute_maxmin_times(){
    double A=0, B=0, C=0, D=0, a=0, b=0, c=0, d=0, disc=0;
    double t1=t_start_, t2=t_end_;
    double P0=coef_[0], P1=coef_[1], P2=coef_[2], P3=coef_[3], P4=coef_[4], P5=coef_[5];
    t_maxjrk_ = (120*P0*t2 - 120*P1*t1 - 480*P1*t2 + 480*P2*t1 + 720*P2*t2 - 720*P3*t1 - 480*P3*t2 + 480*P4*t1 + 120*P4*t2 - 120*P5*t1)/(120*P0 - 600*P1 + 1200*P2 - 1200*P3 + 600*P4 - 120*P5);
    A=  180*P1 -  60*P0 - 180*P2 + 60*P3;
    B= -60*P2  + 180*P3 - 180*P4 + 60*P5;
    C=  120*P1 - 360*P2 + 360*P3 - 120*P4;
    a= A+B+C;
    b= - (2*t2*A + 2*t1*B + t1*C + t2*C);
    c= (A*t2*t2 + B*t1*t1 + C*t1*t2);
    disc= b*b - 4*a*c;
    t_minacc_ = (-b -sqrt(disc) ) / (2*a);
    t_maxacc_ = (-b +sqrt(disc) ) / (2*a);
}


// check if the maximum jerk value is voilated or not 
bool bezier_quintic_segment::check_jerk_limit(){
    double P0=coef_[0], P1=coef_[1], P2=coef_[2], P3=coef_[3], P4=coef_[4], P5=coef_[5];
    double t1=t_start_, t2=t_end_;
    double tm[]= {t1, t2, t_maxjrk_, t_maxacc_, t_minacc_};
    double jrk =0;
     for(auto t : tm){
        jrk = 180*P1*pow( (t-t2), 2) - 60*P0*pow( (t-t2), 2) - 60*P2*pow( (t-t1), 2) - 180*P2*pow( (t-t2), 2) + 180*P3*pow( (t-t1), 2) + 60*P3*pow( (t-t2), 2) - 180*P4*pow( (t-t1), 2) + 60*P5*pow( (t-t1), 2) + 60*P1*(2*t - 2*t2)*(t - t1) - 90*P2*(2*t - 2*t1)*(2*t - 2*t2) + 90*P3*(2*t - 2*t1)*(2*t - 2*t2) - 12*P4*(2*t - 2*t1)*(5*t - 5*t2);
        if(fabs(jrk) > max_jrk_)
           return false;
     }
     return true;
}


// check if the maximum acceleration value is voilated or not 
bool bezier_quintic_segment::check_acc_limit(){
    double P0=coef_[0], P1=coef_[1], P2=coef_[2], P3=coef_[3], P4=coef_[4], P5=coef_[5];
    double t1=t_start_, t2=t_end_;
    double tm[]= {t1, t2, t_maxjrk_, t_maxacc_, t_minacc_};
    double acc =0;
     for(auto t : tm){
         acc = 40*P1*pow( (t-t2), 3) - 20*P0*pow( (t-t2), 3) - 20*P2*pow( (t-t2), 3) + 20*P3*pow( (t-t1), 3) - 40*P4*pow( (t-t1), 3) + 20*P5*pow( (t-t1), 3) + 60*P1*(t - t1)*pow( (t-t2), 2) - 60*P2*(2*t - 2*t1)*pow( (t-t2), 2) - 30*P2*(2*t - 2*t2)*pow( (t-t1), 2) + 30*P3*(2*t - 2*t1)*pow( (t-t2), 2) + 60*P3*(2*t - 2*t2)*pow( (t-t1), 2) - 12*P4*(5*t - 5*t2)*pow( (t-t1), 2);
         if(fabs(acc) > max_acc_)
           return false;
     }
     return true;
}


// check if the maximum velocity value is voilated or not 
bool bezier_quintic_segment::check_vel_limit(){
    double P0=coef_[0], P1=coef_[1], P2=coef_[2], P3=coef_[3], P4=coef_[4], P5=coef_[5];
    double t1=t_start_, t2=t_end_;
    double tm[]= {t1, t2, t_maxjrk_, t_maxacc_, t_minacc_};
    double vel =0;
    for(auto t: tm){
         vel = 5*P1*pow( (t-t2), 4) - 5*P0*pow( (t-t2), 4) - 5*P4*pow( (t-t1), 4) + 5*P5*pow( (t-t1), 4) + 20*P1*(t - t1)*pow( (t-t2), 3) - 10*P2*(2*t - 2*t1)*pow( (t-t2), 3) + 10*P3*(2*t - 2*t2)*pow( (t-t1), 3) - 4*P4*(5*t - 5*t2)*pow( (t-t1), 3) - 30*P2*pow( (t-t1), 2)*pow( (t-t2), 2) + 30*P3*pow( (t-t1), 2)*pow( (t-t2), 2);
        if(fabs(vel) > max_vel_)
            return false;
     }
     return true;
}


// check if the final segment is monotonic or not
bool bezier_quintic_segment::check_monotonic(){
        compute_coef();
        compute_maxmin_times();
        double P0=coef_[0], P1=coef_[1], P2=coef_[2], P3=coef_[3], P4=coef_[4], P5=coef_[5];
        double t1=t_start_, t2=t_end_;
        double t= t_maxjrk_;
        double reached_vel = 5*P1*pow( (t-t2), 4) - 5*P0*pow( (t-t2), 4) - 5*P4*pow( (t-t1), 4) + 5*P5*pow( (t-t1), 4) + 20*P1*(t - t1)*pow( (t-t2), 3) - 10*P2*(2*t - 2*t1)*pow( (t-t2), 3) + 10*P3*(2*t - 2*t2)*pow( (t-t1), 3) - 4*P4*(5*t - 5*t2)*pow( (t-t1), 3) - 30*P2*pow( (t-t1), 2)*pow( (t-t2), 2) + 30*P3*pow( (t-t1), 2)*pow( (t-t2), 2);
        if( (reached_vel >=0 && start_state_[1]>=0  &&  end_state_[1]>=0) || (reached_vel <=0 && start_state_[1]<=0  &&  end_state_[1]<=0)  )
            return true;
        else
            return false;
}



// updated segment duration whenever any maximum value of pos/vel/acc/jrk is voilated 
double  bezier_quintic_segment::update_duration( ){
    if( max_pos_>0 && max_vel_>0 && max_acc_>0 && max_jrk_>0)
        return update_duration(max_pos_, max_vel_, max_acc_, max_jrk_);
    else
        throw( std::invalid_argument("absolute limits (pos, vel, acc, jrk) should be positive and non-zero values") );
}
double  bezier_quintic_segment::update_duration(double max_pos, double max_vel, double max_acc, double max_jrk){
    set_absolute_limits( max_pos,  max_vel,  max_acc,  max_jrk);
    compute_coef();
    compute_maxmin_times();
    if( check_jerk_limit() && check_acc_limit() && check_vel_limit() ){
        std::cout<<"requisted time already grantee jrk, pos, vel limits" <<std::endl;
        return duration_;
    }
    else{
        while (!check_jerk_limit() || !check_acc_limit() || !check_vel_limit()  ){
            set_duration_( duration_ + 0.01);
            compute_coef();
            compute_maxmin_times();
            if( check_jerk_limit() && check_acc_limit() && check_vel_limit() ){
                return duration_;
            }
            iter_ ++;
        }
    }
}


// calculate segment minimum duration
double bezier_quintic_segment::compute_min_duration( ){
    if( max_pos_>0 && max_vel_>0 && max_acc_>0 && max_jrk_>0)
        return compute_min_duration(max_pos_, max_vel_, max_acc_, max_jrk_);
    else
        throw( std::invalid_argument("absolute limits (pos, vel, acc, jrk) should be positive and non-zero values") );
}
double bezier_quintic_segment::compute_min_duration(double max_pos, double max_vel, double max_acc, double max_jrk){
    duration_ = 0.01; //starting with one time sample
    return update_duration();
}


// calculate segment maximum duration
double bezier_quintic_segment::compute_max_duration(){
    if( max_pos_>0 && max_vel_>0 && max_acc_>0 && max_jrk_>0)
        return compute_max_duration(max_pos_, max_vel_, max_acc_, max_jrk_);
    else
        throw( std::invalid_argument("absolute limits (pos, vel, acc, jrk) should be positive and non-zero values") );
}
double bezier_quintic_segment::compute_max_duration(double max_pos, double max_vel, double max_acc, double max_jrk){
    if( !check_generic_transition()){
        std::cout<<"segment transition is from zero state(vel, acc) to zero state(vel , acc), so it does not have max duration (consider 10e5)"<<std::endl;
        return 100000;
    }
    set_absolute_limits(max_pos, max_vel, max_acc, max_jrk);
    duration_ = compute_min_duration();
    compute_coef();
    compute_maxmin_times();
    while(check_monotonic()){
        t_end_ += 0.01;
        compute_coef();
        compute_maxmin_times();
    }
    duration_ = t_end_ - t_start_;
    return t_end_ - t_start_;
}


// check if it is 0 to 0 state
bool bezier_quintic_segment::check_generic_transition(){
    if( fabs(start_state_[1])<0.001 && fabs(start_state_[2])<0.001 && fabs(end_state_[1])<0.001 && fabs(end_state_[2])<0.001 )
        return false;
    else
        return true;
}


// sample segment (find pos/vel/acc/jrk values) at time 
void bezier_quintic_segment::sample_segment(double t, std::vector<double> &state) {
    double t1=t_start_, t2=t_end_;
    double P0=coef_[0], P1=coef_[1], P2=coef_[2], P3=coef_[3], P4=coef_[4], P5=coef_[5];
    state.resize(4); //pos, vel, acc, jrk
    state[0] = P5*pow( (t-t1), 5) - P0*pow( (t-t2), 5) + 5*P1*(t - t1)*pow( (t-t2), 4) - P4*(5*t - 5*t2)*pow( (t-t1), 4) - 10*P2*pow( (t-t1), 2)*pow( (t-t2), 3) + 10*P3*pow( (t-t1), 3)*pow( (t-t2), 2);
    state[1] = 5*P1*pow( (t-t2), 4) - 5*P0*pow( (t-t2), 4) - 5*P4*pow( (t-t1), 4) + 5*P5*pow( (t-t1), 4) + 20*P1*(t - t1)*pow( (t-t2), 3) - 10*P2*(2*t - 2*t1)*pow( (t-t2), 3) + 10*P3*(2*t - 2*t2)*pow( (t-t1), 3) - 4*P4*(5*t - 5*t2)*pow( (t-t1), 3) - 30*P2*pow( (t-t1), 2)*pow( (t-t2), 2) + 30*P3*pow( (t-t1), 2)*pow( (t-t2), 2);
    state[2] = 40*P1*pow( (t-t2), 3) - 20*P0*pow( (t-t2), 3) - 20*P2*pow( (t-t2), 3) + 20*P3*pow( (t-t1), 3) - 40*P4*pow( (t-t1), 3) + 20*P5*pow( (t-t1), 3) + 60*P1*(t - t1)*pow( (t-t2), 2) - 60*P2*(2*t - 2*t1)*pow( (t-t2), 2) - 30*P2*(2*t - 2*t2)*pow( (t-t1), 2) + 30*P3*(2*t - 2*t1)*pow( (t-t2), 2) + 60*P3*(2*t - 2*t2)*pow( (t-t1), 2) - 12*P4*(5*t - 5*t2)*pow( (t-t1), 2);
    state[3] = 180*P1*pow( (t-t2), 2) - 60*P0*pow( (t-t2), 2) - 60*P2*pow( (t-t1), 2) - 180*P2*pow( (t-t2), 2) + 180*P3*pow( (t-t1), 2) + 60*P3*pow( (t-t2), 2) - 180*P4*pow( (t-t1), 2) + 60*P5*pow( (t-t1), 2) + 60*P1*(2*t - 2*t2)*(t - t1) - 90*P2*(2*t - 2*t1)*(2*t - 2*t2) + 90*P3*(2*t - 2*t1)*(2*t - 2*t2) - 12*P4*(2*t - 2*t1)*(5*t - 5*t2);
}

// caclulate all the states of the segment at all sampling point of time 
void bezier_quintic_segment::segment_states(const double frq, std::vector<double> &T_vec, std::vector<double> &POS, std::vector<double> &VEL, std::vector<double> &ACC, std::vector<double> &JRK ){
    double t=t_start_;
    std::vector<double> state;
    while(t < t_end_){
        sample_segment(t, state);
        T_vec.push_back(t);
        POS.push_back(state[0]);
        VEL.push_back(state[1]);
        ACC.push_back(state[2]);
        JRK.push_back(state[3]);
        t+=1/frq;
    }
}

// print attributes values for the segment 
void bezier_quintic_segment::print_attributes(){
    std::cout<<"======================== segment attributes ========================="<<std::endl;
    std::cout<<"### coef_: \n" <<"P0-P5: " <<coef_[0] <<",  " <<coef_[1]<<",  " <<coef_[2]<<",  " <<coef_[3]<<",  " <<coef_[4]<<",  " <<coef_[5]<<std::endl;
    std::cout<<"### limits: \n" <<"pos: " <<max_pos_ <<",  vel: " <<max_vel_<<", acc: " <<max_acc_<<",  jrk: " <<max_jrk_<<std::endl;
    std::cout<<"### times: \n" <<"t_start: " <<t_start_ <<",  t_end: " <<t_end_<<", duration: " <<duration_<<std::endl;
    std::cout<<"### start_state: \n" <<"pos0: " <<start_state_[0] <<",  vel0: " <<start_state_[1]<<", acc0: " <<start_state_[2]<<std::endl;
    std::cout<<"### end_state: \n" <<"pos1: " <<end_state_[0] <<",  vel1: " <<end_state_[1]<<", acc1: " <<end_state_[2]<<std::endl;
    std::cout<<"### iterations: "<<iter_<<"     ### monotonic: "<<check_monotonic()<<std::endl;
    std::cout<<"=========================================================================="<<std::endl;

}

