#ifndef BEZIER_QUINTIC_SEGMENT_H
#define BEZIER_QUINTIC_SEGMENT_H

#include<vector>
#include<math.h>
#include<iostream>

class bezier_quintic_segment
{
private:
    double max_pos_=0, max_vel_=0, max_acc_=0, max_jrk_=0;
    double t_start_=0, duration_=0, t_end_= t_start_+duration_, iter_=0;
    std::vector<double> start_state_, end_state_, coef_;
    double t_maxjrk_=0, t_minacc_=0, t_maxacc_=0;

public:
    bezier_quintic_segment();
    bezier_quintic_segment( double start_time, double duration, const std::vector<double> &start_state, const std::vector<double> &end_state ) {
        init(start_time, duration, start_state, end_state);
    }

    void init(double start_time, double duration, std::vector<double> start_state, std::vector<double> end_state);
    void set_absolute_limits(double max_pos, double max_vel, double max_acc, double max_jrk);


    // compute bezir curve coefficients using initial & final values 
    void compute_coef();    
    void compute_maxmin_times(); // calculate times of inflection points
    
    
    // check if the maximum values (pos, vel, acc, jrk) are voilated or not 
    bool check_jerk_limit();
    bool check_acc_limit();
    bool check_vel_limit();
    bool check_monotonic();
    bool check_generic_transition();


    // updated segment duration whenever any maximum value of pos/vel/acc/jrk is voilated 
    double update_duration();
    double update_duration(double max_pos, double max_vel, double max_acc, double max_jrk);
    void set_duration_(double duration ){ duration_ = duration; t_end_ = t_start_ + duration;  }
    

    // check the minimum and maximum duration a segment can has without voilating maximum values 
    double compute_min_duration();
    double compute_min_duration(double max_pos, double max_vel, double max_acc, double max_jrk);
    double compute_max_duration();
    double compute_max_duration(double max_pos, double max_vel, double max_acc, double max_jrk);


    // segment sampling for plotting purpose 
    void sample_segment(double t, std::vector<double> & state);
    void segment_states(const double frq, std::vector<double> &T_vec, std::vector<double> &POS,
                        std::vector<double> &VEL, std::vector<double> &ACC,
                        std::vector<double> &JRK );


    // print segment attributes in case of debugging
    void print_attributes();
};

#endif // BEZIER_QUINTIC_SEGMENT_H
