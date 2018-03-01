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
static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	//@TODO Set the number of particles
	num_particles = 100;
	// Standard deviations for x, y, and theta
	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];

	//Create a normal (Gaussian) distribution for x
	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);
	
	double sample_x, sample_y, sample_theta;

	// Add random Gaussian noise to each particle.
	for (int i = 0; i < num_particles; i++) {
	    sample_x = dist_x(gen);
		sample_y = dist_y(gen);
		sample_theta = dist_theta(gen);
		Particle p;
		p.x = sample_x;
		p.y = sample_y;
		p.theta = sample_theta;
		p.weight = 1.0;
		particles.push_back(p);
		weights.push_back(p.weight);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// Standard deviations for x, y, and theta
	double std_x = std_pos[0];
	double std_y = std_pos[1];
	double std_theta = std_pos[2];
	
	double sample_x, sample_y, sample_theta;

	for (int i = 0; i < num_particles; i++) {
		Particle particle = particles[i];
		double x, y, theta;
		if (fabs(yaw_rate) < 0.00001) {
			x = particle.x + velocity * delta_t * cos(particles[i].theta);
			y = particle.y + velocity * delta_t * sin(particles[i].theta);
			theta =particle.theta + yaw_rate * delta_t;
		} else {
			x = particle.x + (velocity / yaw_rate) * (sin( particle.theta + yaw_rate * delta_t) - sin(particle.theta));
			y = particle.y + (velocity / yaw_rate) * (cos(particle.theta) - cos(particle.theta + yaw_rate * delta_t));
			theta =particle.theta + yaw_rate * delta_t;
		}
		//Create a normal (Gaussian) distribution for x
		normal_distribution<double> dist_x(x, std_x);
		normal_distribution<double> dist_y(y, std_y);
		normal_distribution<double> dist_theta(theta, std_theta);
		// add random Gaussian noise.
		sample_x = dist_x(gen);
		sample_y = dist_y(gen);
		sample_theta = dist_theta(gen);
		particles[i].x = sample_x;
		particles[i].y = sample_y;
		particles[i].theta = sample_theta;
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	for (int j = 0; j < observations.size(); j++){
		int match_id = -1;
		float min_diff = -1;
		for (int i =0; i < predicted.size(); i++) {
			float diff = dist(observations[j].x, observations[j].y, predicted[i].x, predicted[i].y);
			if (min_diff == -1 || diff < min_diff) {
				min_diff = diff;
				match_id = predicted[i].id;
				//std::cout << "[DEBUG] match_id:" << match_id << std::endl;
			}
			//std::cout << "[DEBUG] match_id2:" << match_id << std::endl;
		}
		if (match_id == -1) {
			std::cout << "[WARNING] data association doesn't match" << std::endl;
		} else {
			observations[j].id = match_id;
		} 	
		//std::cout << "[DEBUG] start SetAssociations" << std::endl;
		//SetAssociations(particle, associations, sense_x, sense_y);
		//std::cout << "[DEBUG] start SetAssociations success" << std::endl;
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
	// filter particles by sensor_range
	double std_x = std_landmark[0];
	double std_y = std_landmark[1];
 	for (int i = 0; i <  particles.size(); i++) {
		double p_x = particles[i].x;
    	double p_y = particles[i].y;
    	double p_theta = particles[i].theta;
		std::vector<LandmarkObs> predictions;
		for (int j = 0 ; j< map_landmarks.landmark_list.size(); j++) {
			Map::single_landmark_s map_landmark = map_landmarks.landmark_list[j];			
			double distant = dist(particles[i].x, particles[i].y, map_landmark.x_f, map_landmark.y_f);
			if (fabs(map_landmark.x_f - particles[i].x) <= sensor_range && fabs(map_landmark.y_f - particles[i].y) <= sensor_range) {
			
				predictions.push_back(LandmarkObs {map_landmark.id_i, map_landmark.x_f, map_landmark.y_f});
				//std::cout << "j: " << j << "map_landmark.x_f: " << map_landmark.x_f
				//<< "map_landmark.id_i: "<< map_landmark.id_i << std::endl;
			}		
		}

		std::vector<LandmarkObs> transformed_os;
		for (int j = 0; j < observations.size(); j++){
			std::vector<double> map_coordinate = transform_coordinate(observations[j].x,
			observations[j].y,  particles[i].theta, particles[i].x, particles[i].y);
			LandmarkObs landmark_obs;
			landmark_obs.x = map_coordinate[0];
			landmark_obs.y = map_coordinate[1];
			landmark_obs.id = observations[j].id;
			//std::cout << "j: " << j << "observations[j].x: " << observations[j].x
			//	<< "landmark_obs.x: "<< landmark_obs.x << std::endl;
			transformed_os.push_back(landmark_obs);
		}
		dataAssociation(predictions, transformed_os);
    	//std::cout << "[DEBUG] start dataAssociation success" << std::endl;

		double p = 1.0;
		for (int j =0; j < transformed_os.size(); j++){
			double x_obs = transformed_os[j].x;
			double y_obs = transformed_os[j].y;
			double mu_x = 0;
			double mu_y = 0;
			bool is_find = false;
			for (int k = 0; k < predictions.size(); k++) {
				LandmarkObs map_landmark = predictions[k];
				//std::cout << "obs[j].id:" << obs[j].id << "map_landmark.id:" << map_landmark.id << std::endl;

				if (map_landmark.id == transformed_os[j].id){
					mu_x = map_landmark.x;
					mu_y = map_landmark.y;
					is_find = true;
				}
			}
			if(!is_find){
				std::cout << "error" << std::endl;
			}
			// calculate exponent
			double exponent= (x_obs - mu_x)*(x_obs - mu_x) /(2 * std_x * std_x) + (y_obs - mu_y) * (y_obs - mu_y)/(2 *std_y * std_y);

			// calculate weight using normalization terms and exponent
			double gauss_norm= (1/(2 * M_PI * std_x * std_y));
			double weight = gauss_norm * exp(-exponent);
			//std::cout  << " j: " << j << " obs[j].id: " << predictions[j].id << "x_obs:" << x_obs 
			//		<< "mu_x:" << mu_x << " exponent: " << exponent 
			//		<< " weights:" << weight << std::endl;
		    p = p * weight ; 
		}

		weights[i] = p;

		particles[i].weight = p;

	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
  vector<Particle> new_particles;

  // get all of the current weights

  // generate random starting index for resampling wheel
  uniform_int_distribution<int> uniintdist(0, num_particles-1);
  auto index = uniintdist(gen);

  // get max weight
  double max_weight = *max_element(weights.begin(), weights.end());

  // uniform random distribution [0.0, max_weight)
  uniform_real_distribution<double> unirealdist(0.0, max_weight);

  double beta = 0.0;

  // spin the resample wheel!
  for (int i = 0; i < num_particles; i++) {
    beta += unirealdist(gen) * 2.0;
    while (beta > weights[index]) {
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    new_particles.push_back(particles[index]);
  }

  particles = new_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
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
