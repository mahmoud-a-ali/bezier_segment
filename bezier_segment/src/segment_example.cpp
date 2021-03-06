#include <ros/ros.h>
#include "bezier_segment/bezier_quintic_segment.h"
#include <yaml_trajectory_1.h>

#include <python2.7/Python.h>
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
using namespace std;



const double max_pos=180, max_vel= 130, max_acc=250, max_jrk=1000, frq=125; // pos, vel, acc, jrk max limits


int main(int argc, char **argv)
{
    ros::init(argc, argv, "bezier_segment_example");
    ros::NodeHandle nh("~");
    ros::console::set_logger_level(ROSCONSOLE_DEFAULT_NAME, ros::console::levels::Debug);
    ROS_INFO("initializing path .... ");

    // load trajectory from header file
    trajectory_msgs::JointTrajectory  traj;
    traj = generate_traj();
    int n_jts = traj.joint_names.size();
    int n_pts = traj.points.size();

    // extract positions, velocities ,accelerations
    std::vector< std::vector<double> > P_jt_wpt, V_jt_wpt, A_jt_wpt,T_jt_wpt;
    P_jt_wpt.resize(n_jts);
    V_jt_wpt.resize(n_jts);
    A_jt_wpt.resize(n_jts);
    T_jt_wpt.resize(n_jts);
    for(int jt=0; jt<n_jts; jt++){
        for(int pt=0; pt<n_pts; pt++){
            P_jt_wpt[jt].push_back( traj.points[pt].positions[jt] );
            V_jt_wpt[jt].push_back( traj.points[pt].velocities[jt] );
            A_jt_wpt[jt].push_back( traj.points[pt].accelerations[jt] );
            T_jt_wpt[jt].push_back(0);
        }
    }


    //choose a segment from the loaded trajectory
    int jt=4, seg=6;
    nh.getParam("segment", seg);

    std::vector<double> start_state, end_state;
    start_state.resize(3);
    end_state.resize(3);

    ROS_INFO_STREAM("=========================");
    ROS_INFO_STREAM("       seg: "<< seg);
    ROS_INFO_STREAM("=========================");

    // assign start state 
    start_state[0]= P_jt_wpt[jt][seg] ;
    start_state[1]=V_jt_wpt[jt][seg] ;
    start_state[2]= A_jt_wpt[jt][seg] ;
    // assign end state
    end_state[0]= P_jt_wpt[jt][seg+1] ;
    end_state[1]= V_jt_wpt[jt][seg+1] ;
    end_state[2]= A_jt_wpt[jt][seg+1] ;
    // assign starting time and initial_duration
    double old_duration=0.01, t_start = T_jt_wpt[jt][seg];
    
    //create bezier segment object
    bezier_quintic_segment traj_seg( t_start, old_duration, start_state, end_state );
    double new_duration = traj_seg.update_duration(max_pos, max_vel, max_acc, max_jrk);
    T_jt_wpt[jt][seg+1] = t_start + new_duration;
    traj_seg.print_attributes();

    // new vectors for ploting purpose 
    std::vector<double> T_vec, POS, VEL,ACC, JRK;
    traj_seg.segment_states( frq, T_vec, POS, VEL,ACC, JRK );


    // plot segment
    plt::figure(1);
    plt::subplot(2, 2, 1);
    plt::named_plot( "pos",T_vec, POS); 
    plt::title("pos"); plt::grid(true); plt::title("pos"); //plt::legend();
    plt::subplot(2, 2, 2);
    plt::named_plot( "vel",T_vec, VEL); plt::grid(true); plt::title("vel");
    plt::subplot(2, 2, 3);
    plt::named_plot( "acc",T_vec, ACC); plt::grid(true); plt::title("acc");
    plt::subplot(2, 2, 4);
    plt::named_plot( "jrk",T_vec, JRK); plt::grid(true); plt::title("jrk");
    plt::show();


}
