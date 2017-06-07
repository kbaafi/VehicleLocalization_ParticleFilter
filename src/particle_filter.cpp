/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"
#include "helper_functions.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles = 200;

	// Normal distributions for x position, y position and heading theta
	default_random_engine gen;

	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y,std[1]);
	normal_distribution<double> dist_theta(theta,std[2]);

	for(int i=0;i<num_particles;i++){
		Particle p;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1.0;
		p.id = i;

		particles.push_back(p);
	}

	is_initialized = true;

	return;	
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;

	for (auto& i:particles){
		double theta = i.theta;
		double x = i.x;
		double y = i.y;
		
		/*
		if(fabs(yaw_rate)>0.01){
			std::cout<<"\n"<<"big fabs yaw_rate"<<std::endl;
			x+=((velocity/yaw_rate)*(sin(theta+(yaw_rate*delta_t))-sin(theta)));
			y+=	((velocity/yaw_rate)*(cos(theta)-cos(theta+yaw_rate*delta_t)));
			theta+=(yaw_rate*delta_t);
		}
		else{
			x+=velocity*delta_t*cos(theta);
			y+=velocity*delta_t*sin(theta);
		}
		*/

		x+=velocity * delta_t * cos(theta);
		y+=velocity * delta_t * sin(theta);
		theta+=(yaw_rate*delta_t);


		// adding requisite white noise for predicted x,y, and theta
		normal_distribution<double> dist_x(x, std_pos[0]);
		normal_distribution<double> dist_y(y,std_pos[1]);
		normal_distribution<double> dist_theta(theta,std_pos[2]);

		i.x = dist_x(gen);
		i.y = dist_y(gen);
		i.theta = dist_theta(gen);
	}
	return;
}

/*
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	return;
}
*/

std::vector<LandmarkAssociation> ParticleFilter::dataAssociation(std::vector<LandmarkObs> observations,Map mapLandmarks){

	std::vector<LandmarkAssociation> assoc;

	for(auto& i:observations)  {

		double predicted_x = i.x;
		double predicted_y = i.y;
		double distance =0.00;
		double curr_distance = 999999;//dist(predicted_x,predicted_y,mapLandmarks.landmark_list.at(0).x_f,mapLandmarks.landmark_list.at(0).y_f);

		LandmarkAssociation assignment;
		assignment.predicted_x = predicted_x;
		assignment.predicted_y = predicted_y;

		//redundant
		int selected_map_landmark_id; //= mapLandmarks.landmark_list.at(0).id_i;
		//assignment.map_landmark_id = selected_map_landmark_id;

		for(auto& lm:mapLandmarks.landmark_list){	
			double map_x = lm.x_f;
			double map_y = lm.y_f;

			distance = dist(predicted_x,predicted_y,map_x,map_y);

			if(distance<curr_distance){
				assignment.map_x = map_x;
				assignment.map_y = map_y;
				curr_distance = distance;
				selected_map_landmark_id = lm.id_i;
				assignment.map_landmark_id = selected_map_landmark_id;
			}
		}
		
		assoc.push_back(assignment);
	}
	return assoc;
}



void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	for(auto& particle:particles){
		std::vector<LandmarkObs> obsInRange;

		double px = particle.x;
		double py = particle.y;
		double ptheta = particle.theta;

		for(auto& i:observations){
			//transform the observations
			LandmarkObs tObs;
			
			double lidar_x = i.x;
			double lidar_y = i.y;
			
			tObs = particle_transform_lidar_to_map(px,py,ptheta,lidar_x,lidar_y);

			double d = dist(px,py,tObs.x,tObs.y);

			if(d<=sensor_range){			
				obsInRange.push_back(tObs);
			}
		}
		
		// perform data association
		std::vector<LandmarkAssociation> landmarkAssocs;
		landmarkAssocs = ParticleFilter::dataAssociation(obsInRange,map_landmarks);

				// saving information for simulator
		particle.associations.clear();
		particle.sense_x.clear();
		particle.sense_y.clear();

		for (auto& assoc:landmarkAssocs){
			particle.associations.push_back(assoc.map_landmark_id);
			particle.sense_x.push_back(assoc.map_x);
			particle.sense_y.push_back(assoc.map_y);
		}

		// calculate multivariate gaussian
		particle.weight = particle_multivariate_location_probability(landmarkAssocs, std_landmark);	
	}	
	
	return;
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	std::random_device rd;
	std::mt19937 gen(rd());

	std::vector<double> weights;
	for(auto& el:particles){
		weights.push_back(el.weight);
	}

	std::discrete_distribution<> d(weights.begin(),weights.end());

	std::vector<Particle> new_particle_set;

	for (int i=0;i<num_particles;++i){
		Particle p = particles.at(d(gen));
		p.id = i;
        new_particle_set.push_back(p);
	
  	}

	// update with optimized particles
	particles = new_particle_set;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
