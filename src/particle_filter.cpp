/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <iostream>
#include <sstream>

#include "particle_filter.h"


using namespace std;

default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
    //   x, y, theta and their uncertainties from GPS) and all weights to 1.
    // Add random Gaussian noise to each particle.
    // NOTE: Consult particle_filter.h for more information about this method (and others in this file).


    //	initialize the numner of particles .
    num_particles = 50;

    //	applying sensor noise
    // This line creates a normal (Gaussian) distribution for x
    normal_distribution<double> normal_dist_x(0, std[0]);
    normal_distribution<double> normal_dist_y(0, std[1]);
    normal_distribution<double> normal_dist_psi(0, std[2]);


    for (int i = 0; i < num_particles; ++i) {
        // intiallizing particle object and its values
        // all the normal distribution represent the noise here.
        Particle particle;
        particle.id = i;
        particle.x = x + normal_dist_x(gen);
        particle.y = y + normal_dist_y(gen);
        particle.theta = theta + normal_dist_psi(gen);
        particle.weight = 1.0;

        //adding additional noise
        particle.x = particle.x + normal_dist_x(gen);
        particle.y = particle.y + normal_dist_y(gen);
        particle.theta = particle.theta + normal_dist_psi(gen);

        particles.push_back(particle);
    }

    is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // TODO: Add measurements to each particle and add random Gaussian noise.
    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/
    normal_distribution<double> normal_dist_x(0, std_pos[0]);
    normal_distribution<double> normal_dist_y(0, std_pos[1]);
    normal_distribution<double> normal_dist_psi(0, std_pos[2]);

    for (int i = 0; i < num_particles; i++) {
        Particle particle = particles[i];

        if (fabs(yaw_rate) < 0.00001) {
            particle.x += velocity * delta_t * cos(particle.theta);
            particle.y += velocity * delta_t * sin(particle.theta);
        }
        else {
            particle.x += velocity / yaw_rate * (sin(particle.theta + yaw_rate*delta_t) - sin(particle.theta));
            particle.y += velocity / yaw_rate * (cos(particle.theta) - cos(particle.theta + yaw_rate * delta_t));
            particle.theta += yaw_rate * delta_t;
        }

        //adding additional noise
        particle.x = particle.x + normal_dist_x(gen);
        particle.y = particle.y + normal_dist_y(gen);
        particle.theta = particle.theta + normal_dist_psi(gen);
    }


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs> &observations) {
    // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
    //   observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
    //   implement this method and use it as a helper during the updateWeights phase.

//    https://stackoverflow.com/questions/409348/iteration-over-stdvector-unsigned-vs-signed-index-variable |iterate over a vector
    for (unsigned i=0; i < observations.size(); i++) { // itrating over each landmark
        LandmarkObs landmarkObsObservation = observations[i];


        int tempId = -1;
        //initially setting to max
        double minDistance = numeric_limits<double>::max() ;//http://en.cppreference.com/w/cpp/types/numeric_limits/max

        for (unsigned j = 0; j < predicted.size(); j++) { // itrating over each prediction over each landmark
            LandmarkObs landmarkObsPredicted = predicted[j];

            double euclidianDistance = dist(landmarkObsObservation.x,landmarkObsObservation.y,landmarkObsPredicted.x,landmarkObsPredicted.y);
            if (euclidianDistance < minDistance){
                minDistance = euclidianDistance;
                tempId = landmarkObsPredicted.id;

            }
        }
        landmarkObsObservation.id = tempId;
    }
};

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

    // for each particle...
    for (int i = 0; i < num_particles; i++) {

        Particle particle = particles[i];

        // create a vector to hold the map landmark locations predicted to be within sensor range of the particle
        vector<LandmarkObs> listOfPredictions;

        // itrating each landmark...
        for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {

            // getting landmark details of  id and x,y coordinates
            float singleLandmarkX = map_landmarks.landmark_list[j].x_f;
            float singleLandmarkY = map_landmarks.landmark_list[j].y_f;
            int singleLandmarkId = map_landmarks.landmark_list[j].id_i;

            // only consider landmarks within sensor range of the particle (rather than using the "dist" method considering a circular
            // region around the particle, this considers a rectangular region but is computationally faster)
            if (fabs(singleLandmarkX - particle.x) <= sensor_range && fabs(singleLandmarkY - particle.y) <= sensor_range) {

                // add prediction to vector
                listOfPredictions.push_back(LandmarkObs{ singleLandmarkId, singleLandmarkX, singleLandmarkY });
            }
        }

        //converted vehicle Coordinates to map coordinates
        vector<LandmarkObs> transformedGroundtruthToMapCoordinates;

        // performing conversion for each car/ground truth coordinates
        for (unsigned int j = 0; j < observations.size(); j++) {

            double t_x = cos(particle.theta)*observations[j].x - sin(particle.theta)*observations[j].y + particle.x;
            double t_y = sin(particle.theta)*observations[j].x + cos(particle.theta)*observations[j].y + particle.y;
            transformedGroundtruthToMapCoordinates.push_back(LandmarkObs{ observations[j].id, t_x, t_y });
        }

        dataAssociation(listOfPredictions, transformedGroundtruthToMapCoordinates);

        particle.weight = 1.0;

        for (unsigned int j = 0; j < transformedGroundtruthToMapCoordinates.size(); j++) {

            // placeholders for observation and associated prediction coordinates
            double observedX, observedY, predictedX, predictedY;
            observedX = transformedGroundtruthToMapCoordinates[j].x;
            observedY = transformedGroundtruthToMapCoordinates[j].y;

            int associatedPrediction = transformedGroundtruthToMapCoordinates[j].id;


            for (unsigned int k = 0; k < listOfPredictions.size(); k++) {
                if (listOfPredictions[k].id == associatedPrediction) {
                    predictedX = listOfPredictions[k].x;
                    predictedY = listOfPredictions[k].y;
                }
            }

            // calculate weight for this observation with multivariate Gaussian
            double stdLandmarkX = std_landmark[0];
            double stdLandmarkY = std_landmark[1];
            double observedWeight = ( 1 / (2 * M_PI * stdLandmarkX * stdLandmarkY)) * exp( -( pow(predictedX - observedX , 2)/( 2 * pow(stdLandmarkX , 2)) + (pow(predictedY - observedY , 2)/(2 * pow(stdLandmarkY, 2))) ) );

            // product of this obersvation weight with total observations weight
            particle.weight *= observedWeight;
        }
    }


}

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
//Reference  python code
//    def resample(myrobot, p):
//    Z = myrobot.sense()
//    w = []
//    for i in range(N):
//          w.append(p[i].measurement_prob(Z))
//
//    p3 = []
//    index = int(random.random() * N)
//    beta = 0.0
//    mw = max(w)
//    for i in range(N):
//      beta += random.random() * 2.0 * mw
//      while beta > w[index]:
//          beta -= w[index]
//          index = (index + 1) % N
//      p3.append(p[index])
//    p = p3


    // get all of the current weights
    vector<double> listWeights;

    for (int i = 0; i < num_particles; i++) {
        listWeights.push_back(particles[i].weight);
    }

    // generate random starting index for resampling wheel
    uniform_int_distribution<int> uniintdist(0, num_particles-1);

    auto index = uniintdist(gen);

    // get max weight  in the
    double max_weight = *max_element(listWeights.begin(), listWeights.end());


    vector<Particle> tempListParticles;

    // uniform random distribution instead of random distribution.
    uniform_real_distribution<double> unirealdist(0.0, max_weight);

    double beta = 0.0;

    // spin the resample wheel!
    for (int i = 0; i < num_particles; i++) {
        beta += unirealdist(gen) * 2.0;
        while (beta > listWeights[index]) {
            beta -= listWeights[index];
            index = (index + 1) % num_particles;
        }
        tempListParticles.push_back(particles[index]);
    }

    particles = tempListParticles;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x,
                                         std::vector<double> sense_y) {
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    //Clear the previous associations
    particle.associations.clear();
    particle.sense_x.clear();
    particle.sense_y.clear();

    particle.associations = associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;

    return particle;
}

string ParticleFilter::getAssociations(Particle best) {
    vector<int> v = best.associations;
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length() - 1);  // get rid of the trailing space
    return s;
}

string ParticleFilter::getSenseX(Particle best) {
    vector<double> v = best.sense_x;
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length() - 1);  // get rid of the trailing space
    return s;
}

string ParticleFilter::getSenseY(Particle best) {
    vector<double> v = best.sense_y;
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length() - 1);  // get rid of the trailing space
    return s;
}
