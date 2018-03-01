/*
 * particle_filter.cpp
 *
 *  Created on: March 1, 2018
 *      Author: Hiu Chan
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

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
  //   x, y, theta and their uncertainties from GPS) and all weights to 1. 
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).

  // I tried a few different number for the number of particles and even 5 will pass the project (3 will fail, I did not try 4.). However, look like the larger it is
  // The smaller the errors. For 10, I got x = 0.154 and y = 0.147, which is not bad. If I put 70, my error was down to about 0.11. Not sure it is worth it for the
  // this kind of improvement with more computation resource waste.
  num_particles = 10;
  is_initialized = true;
  
  particles.clear();

  default_random_engine generator;
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  for (int i=0; i<num_particles; ++i){
    Particle p;
    p.id = i;
    p.x = dist_x(generator);
    p.y = dist_y(generator);
    p.theta = dist_theta(generator);
    p.weight = 1;
    particles.push_back(p);
  }
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // TODO: Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/

  default_random_engine generator;

  double distance = velocity*delta_t;
  for(std::vector<Particle>::iterator it = particles.begin(); it != particles.end(); ++it) {
    double new_theta;
    double new_x;
    double new_y;

    if (fabs(yaw_rate) != 0){
      new_theta = it->theta + yaw_rate * delta_t;
      new_x = it->x + (velocity/yaw_rate) * (sin(new_theta) - sin(it->theta));
      new_y = it->y + (velocity/yaw_rate) * (cos(it->theta) - cos(new_theta));
    }
    else{ // if yaw_rate is zero, just do the regular tri-geometry since the car is moving at a straight line
      new_theta = it->theta;
      new_x = it->x + (velocity*delta_t) * cos(it->theta);
      new_y = it->y + (velocity*delta_t) * sin(it->theta);
    }

    normal_distribution<double> dist_x(new_x, std_pos[0]);
    normal_distribution<double> dist_y(new_y, std_pos[1]);
    normal_distribution<double> dist_theta(new_theta, std_pos[2]);
    
    it->x = dist_x(generator);
    it->y = dist_y(generator);
    it->theta = dist_theta(generator);
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.


  // Don't really think I need this function but keep it here anyway.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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

  double gau_front = 2 * M_PI * std_landmark[0] * std_landmark[1] ;
  gau_front = 1/gau_front;

  weights.clear();

  // loop through each particle first
  for(std::vector<Particle>::iterator it = particles.begin(); it != particles.end(); ++it) {
    
    it->associations.clear();
    it->sense_x.clear();
    it->sense_y.clear();
    it->weight = 1;

    // for each particle, we check the observed landmarks and convert to the map coordiate.
    for (std::vector<LandmarkObs>::const_iterator lm_it = observations.begin(); lm_it != observations.end(); ++lm_it) {      
      double sin_ang = sin(it->theta);
      double cos_ang = cos(it->theta);

      double map_x = lm_it->x * cos_ang - lm_it->y * sin_ang + it->x ;
      double map_y = lm_it->x * sin_ang + lm_it->y * cos_ang + it->y ;

      int asso = 0;
      double distance = 999999999999999;
      double mu_x , mu_y;

      // For each observation, we find the association of the landmark on the map base on the smallest distance.
      for (std::vector<Map::single_landmark_s>::const_iterator map_it = map_landmarks.landmark_list.begin(); map_it != map_landmarks.landmark_list.end(); ++map_it){
	double d = dist(map_x, map_y, map_it->x_f, map_it->y_f);
	if (d  < distance){
	  distance = d;
	  asso = map_it->id_i;
	  mu_x = map_it->x_f;
	  mu_y = map_it->y_f;		
	}
      }
      double exponent = ((mu_x - map_x)*(mu_x - map_x))/(2*std_landmark[0]*std_landmark[0]) 
			 + ((mu_y - map_y)*(mu_y - map_y))/(2*std_landmark[1]*std_landmark[1]);
      exponent = -exponent;
      it->weight *= gau_front * exp(exponent);

      it->associations.push_back(asso);
      it->sense_x.push_back(map_x);
      it->sense_y.push_back(map_y);
    }
    weights.push_back(it->weight);
  }  
}

void ParticleFilter::resample() {
  // TODO: Resample particles with replacement with probability proportional to their weight. 
  // NOTE: You may find std::discrete_distribution helpful here.
  //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  vector<Particle> resampled_particles;

  default_random_engine gen;
  discrete_distribution<int> weight_distribution(weights.begin(), weights.end()); 

  for (int i = 0; i < num_particles; i++) {
    resampled_particles.push_back(particles[weight_distribution(gen)]);
  }

  particles = resampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
