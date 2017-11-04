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
    Particle sample;
    default_random_engine gen;
    int particle_num = 300;
	// This line creates a normal (Gaussian) distribution for x
	normal_distribution<double> dist_x(x, std[0]);
	
	// TODO: Create normal distributions for y and theta
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	
	for (int i = 0; i < particle_num; ++i) {
		// TODO: Sample  and from these normal distrubtions like this: 
		//	 sample_x = dist_x(gen);
		//	 where "gen" is the random engine initialized earlier.
	    sample.id = i;	
		sample.x = dist_x(gen);
		sample.y = dist_y(gen);
		sample.theta = dist_theta(gen);
        sample.weight = 1.0;
        particles.push_back(sample);
        weights.push_back(1.0);
	}
    num_particles = particles.size();
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    double x, y, theta;
    for (int i=0; i<num_particles; i++) {
        x = particles[i].x;
        y = particles[i].y;
        theta = particles[i].theta;
        if (!yaw_rate) {
           x +=  velocity * sin(theta) * delta_t;
           y +=  velocity * cos(theta) * delta_t;
        }
        else {
           x += velocity/yaw_rate * (sin(theta + yaw_rate * delta_t) - sin(theta));
           y += velocity/yaw_rate * (cos(theta) - cos(theta + yaw_rate * delta_t));
           theta +=  yaw_rate * delta_t;
        }
	
    	normal_distribution<double> dist_x(x, std_pos[0]);
    	normal_distribution<double> dist_y(y, std_pos[1]);
    	normal_distribution<double> dist_theta(theta, std_pos[2]);
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);
        particles[i].id = i;
    }
    return;
}
/*
double ParticleFilter::CalculateLandmarkWeight(double sigma_x_2sq, double sigma_y_2sq, double outer_term, LandmarkObs observed_lm, LandmarkObs predicted_lm) const
{
	auto x_term = pow(predicted_lm.x - observed_lm.x, 2) / sigma_x_2sq;
	auto y_term = pow(predicted_lm.y - observed_lm.y, 2) / sigma_y_2sq;
	auto pow_term = -(x_term + y_term);
	return outer_term * exp(pow_term);
}

bool ParticleFilter::CheckLandmarkRange(LandmarkObs& observed_lm, LandmarkObs*& predicted_lm, double& min_distance, LandmarkObs current_landmark) const
{
	auto distance = pow(current_landmark.x - observed_lm.x, 2) + pow(current_landmark.y - observed_lm.y, 2);
	if (distance < min_distance || min_distance == -1)
	{
		min_distance = distance;
		observed_lm.id = current_landmark.id;
		predicted_lm = &current_landmark;

		// If the particle is closer than a certain threshold, end the loop early for faster process.
		if (min_distance <= 1.0) return true;
	}
	return false;
}
*/
std::vector<LandmarkObs> ParticleFilter::observationtoMap(std::vector<LandmarkObs> &observations,
        double x, double y, double theta)
{
    std::vector<LandmarkObs> observations_map;
    LandmarkObs observation_m;
    for (int k = 0; k < observations.size(); k++) {
        observation_m.x = observations[k].x * cos(theta) - observations[k].y * sin(theta) + x;
        observation_m.y = observations[k].x * sin(theta) + observations[k].y * cos(theta) + y;
        observation_m.id = observations[k].id;
        observations_map.push_back(observation_m);
    }
    return observations_map;
}

double ParticleFilter::CalculateParticleWeight(double sensor_range, vector<LandmarkObs> observations, Map map_landmarks, double sigma_x_2sq, double sigma_y_2sq, double outer_term, Particle& particle)
{
	particle.weight = 1;
    double x = particle.x;
    double y = particle.y;
    double theta = particle.theta;
    std::vector<LandmarkObs> predictes;
    observations = observationtoMap(observations,x,y,theta);
    landmarkinrange(predictes, sensor_range, map_landmarks, particle);
    dataAssociation(predictes, observations);
    for (int k = 0; k < observations.size(); k++) {
        double x_m = observations[k].x;
        double y_m = observations[k].y;
        double id = observations[k].id;
        double x_l = map_landmarks.landmark_list[id].x_f;
        double y_l = map_landmarks.landmark_list[id].y_f;
        
        
	    double x_term = pow(x_l - x_m, 2) / sigma_x_2sq;
	    double y_term = pow(y_l - y_m, 2) / sigma_y_2sq;
	    double pow_term = -(x_term + y_term);
        particle.weight*=outer_term*exp(-((pow((x_l-x_m),2)/sigma_x_2sq)+(pow(y_l-y_m,2)/sigma_y_2sq)));
    }
	return particle.weight;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   vector<LandmarkObs> &observations, const Map &map_landmarks)
{
	auto sigma_x = std_landmark[0];
	auto sigma_y = std_landmark[1];
	auto sigma_x_2sq = 2 * sigma_x * sigma_x;
	auto sigma_y_2sq = 2 * sigma_y * sigma_y;
	auto outer_term = 1 / (2 * M_PI * sigma_x * sigma_y);

	weights.clear();
	for (int i = 0, m = particles.size(); i < m; i++)
        weights.push_back(CalculateParticleWeight(sensor_range, observations, map_landmarks, sigma_x_2sq, sigma_y_2sq,outer_term, particles[i]));
    return;
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs> &observations) {
	double min_dist = 100000;
    double dist1 = 0;
    for ( int j = 0; j < observations.size(); j++) {
        min_dist = 100000;
        for ( int i = 0; i < predicted.size(); i++) {
            dist1 = dist(predicted[i].x, predicted[i].y, observations[j].x,observations[j].y);
            if (min_dist > dist1&&dist1<=1.0) {
                min_dist = dist1;
                observations[j].id = predicted[i].id;
                break;
            }
        }
    }
    return ;
}



void ParticleFilter::landmarkinrange(std::vector<LandmarkObs> &predictes, double sansor_range,
        const Map &map_landmarks,Particle particle )
{
    double l_x,l_y,p_x,p_y;
    double dist;
    p_x = particle.x;
    p_y = particle.y;
    LandmarkObs predict_tmp;
    auto num_landmarks = map_landmarks.landmark_list.size();
    for ( int i = 0; i < num_landmarks; i++) {
        l_x = map_landmarks.landmark_list[i].x_f;
        l_y = map_landmarks.landmark_list[i].y_f;
        dist = sqrt(pow((p_x-l_x),2)+pow((p_y-l_y),2));
        if (dist < 1.5*sansor_range) {
            predict_tmp.x = l_x;
            predict_tmp.y = l_y;
            predict_tmp.id = i; 
            predictes.push_back(predict_tmp);
        }
    }
    return ;
}

void ParticleFilter::resample() {
    discrete_distribution<> distribution(weights.begin(),weights.end());
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd());
    std::vector<Particle> particles_temp;
    int index = (int)distribution(gen);
    double bate = 0.0;
    double mw = 0;
    for (int i = 0; i < weights.size(); i++)
        if(mw < weights[i])
            mw = weights[i];
    for(int i = 0; i < num_particles; i++) {
       bate += (double)distribution(gen)/10.0 *2.0 * mw;
       while (bate>weights[index]){
          bate -=weights[index];
          index = (index+1) % num_particles;
       }
       particles_temp.push_back(particles[index]);
    }
    particles.clear();
    for( int i = 0; i< num_particles; i++) {
        particles.push_back(particles_temp[i]);
    }
    return;
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
