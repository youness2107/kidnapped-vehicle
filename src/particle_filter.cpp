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

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 100;

        std::default_random_engine gen;

        std::normal_distribution<double> N_x(x, std[0]);
        std::normal_distribution<double> N_y(y, std[1]);
        std::normal_distribution<double> N_theta(theta, std[2]);

        for (int i = 0; i < num_particles; i++) {

		Particle particle;
		particle.id = i;
		particle.x = N_x(gen);
		particle.y = N_y(gen);
                particle.theta = N_theta(gen);
                particle.weight = 1;

                particles.push_back(particle);
                weights.push_back(1);
         }

	is_initialized = true;



}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;

	for (int i = 0; i < num_particles; i++) {
            
		double p_x = particles[i].x;
		double p_y = particles[i].y;
                double p_theta = particles[i].theta;

           	double new_x;
	   	double new_y;
		double new_theta;

		if (yaw_rate == 0) {
			new_x = p_x + velocity*delta_t*cos(p_theta);
                        new_y = p_y + velocity*delta_t*sin(p_theta);
                        new_theta = p_theta;
               }
               else {
                       new_x = p_x + (velocity/yaw_rate)*(sin(p_theta + yaw_rate*delta_t) - sin(p_theta));
                       new_y = p_y + (velocity/yaw_rate)*(cos(p_theta) - cos(p_theta + (yaw_rate * delta_t)));
                       new_theta = p_theta + yaw_rate*delta_t;
               }
		normal_distribution<double> N_x(new_x, std_pos[0]);
		normal_distribution<double> N_y(new_y, std_pos[1]);
		normal_distribution<double> N_theta(new_theta, std_pos[2]);
	
                particles[i].x = N_x(gen);
		particles[i].y = N_y(gen);
		particles[i].theta = N_theta(gen);
       }


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations, double sensor_range) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.


	int i, j;
	for (i = 0; i < observations.size(); i++) {
	//Maximum distance can be square root of 2 times the range of sensor.
		double lowest_distance = sensor_range;
		double observation_x = observations[i].x;
		double observation_y = observations[i].y;
		int closest_landmark_id = -1;

		for (j = 0; j < predicted.size(); j++) {
		  double predicted_x = predicted[j].x;
		  double predicted_y = predicted[j].y;
		  int predicted_id = predicted[j].id;
		  double current_distance = dist(observation_x, observation_y, predicted_x, predicted_y);

		  if (lowest_distance > current_distance) {
		    lowest_distance = current_distance;
		    closest_landmark_id = predicted_id;
		  }
		}
		observations[i].id = closest_landmark_id;
	}


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



  // normalize weight with this
  double normalizer = 0.0;
  int i, j;


  for (i = 0; i < num_particles; i++) {
    double particle_x = particles[i].x;
    double particle_y = particles[i].y;
    double particle_theta = particles[i].theta;

    // Vehicle coordinates to Map coordiantes transformation
    vector<LandmarkObs> transformed_observations;

    for (j = 0; j < observations.size(); j++) {
      LandmarkObs transformed_obs;
      transformed_obs.id = j;
      transformed_obs.x = particle_x + (cos(particle_theta) * observations[j].x) - (sin(particle_theta) * observations[j].y);
      transformed_obs.y = particle_y + (sin(particle_theta) * observations[j].x) + (cos(particle_theta) * observations[j].y);
      transformed_observations.push_back(transformed_obs);
    }

    // keep only ones in sensor_range 
    vector<LandmarkObs> in_range_landmarks;

    for (int j = 0;  j < map_landmarks.landmark_list.size(); j++) {

      const int mid = map_landmarks.landmark_list[j].id_i;
      const double mx = map_landmarks.landmark_list[j].x_f;
      const double my = map_landmarks.landmark_list[j].y_f;

      const double dx = mx - particle_x;
      const double dy = my - particle_y;
      const double error = sqrt(dx * dx + dy * dy);

      if (error < sensor_range) {

        LandmarkObs in_range_landmark = {
          mid,
          mx,
          my
         };

        in_range_landmarks.push_back(in_range_landmark);
      }
    }

    //Associate observations
    dataAssociation(in_range_landmarks, transformed_observations, sensor_range);

    // Weight Calculation
  
    particles[i].weight = 1.0;

	
    double sigma_x = std_landmark[0];
    double sigma_y = std_landmark[1];
    double sigma_x_2 = pow(sigma_x, 2);
    double sigma_y_2 = pow(sigma_y, 2);
	
    int k, l;

    for (k = 0; k < transformed_observations.size(); k++) {
     
      double prob = 1.0;

      for (l = 0; l < in_range_landmarks.size(); l++) {
        double pred_landmark_x = in_range_landmarks[l].x;
        double pred_landmark_y = in_range_landmarks[l].y;
        double pred_landmark_id = in_range_landmarks[l].id;

        if (transformed_observations[k].id == pred_landmark_id) {
          prob = (1.0/(2.0 * M_PI * sigma_x * sigma_y)) * exp(-1.0 * ((pow((transformed_observations[k].x - pred_landmark_x), 2)/(2.0 * sigma_x_2)) + (pow((transformed_observations[k].y - pred_landmark_y), 2)/(2.0 * sigma_y_2))));
          particles[i].weight *= prob;
        }
      }
    }

  normalizer += particles[i].weight;
  }

  // Normalize the weights
  for (int i = 0; i < particles.size(); i++) {
    double normalized_weight = particles[i].weight/normalizer;
    particles[i].weight = normalized_weight;
    weights[i] = normalized_weight;
  }


}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;
        discrete_distribution<int> distribution(weights.begin(), weights.end());

        vector<Particle> resample_particles;

        for (int i = 0; i < num_particles; i++) {
            resample_particles.push_back(particles[distribution(gen)]);

        }

	particles = resample_particles;

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
