#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <omp.h>
#include <sstream>
#include <algorithm>
#include <limits>
#include "nqueue/queue.hpp"

class VisitedCity {
	public: 
		std::vector<int> _tour;
		float _cost;
		float _lowerBound;
		int _length;
		int _node;

        virtual ~VisitedCity() = default;
		VisitedCity() = default;

        struct CompareCityByLowerBound {
            bool operator()(const VisitedCity& a, const VisitedCity& b) {
                return a._lowerBound > b._lowerBound;
            }
        };

		VisitedCity(std::vector<int> tour, float cost, float lowerBound, int length, int node) {
			_tour = tour;
			_cost = cost;
			_lowerBound = lowerBound;
			_length = length;
			_node = node;
		}

		std::vector<int> getTour() {
			return _tour;
		}
		float getCost() {
			return _cost;
		}
		float getLB() {
			return _lowerBound;
		}
		int getLength() {
			return _length;
		}
		int getNode() {
			return _node;
		}
};

class TSP {
    public:
        int _numberOfCities;
        int _numberOfRoads;
        float _maxValue;
        float _bestTourCost = 999999999;
        std::vector<int> _bestTour;
        std::vector<std::vector<float>> _roadsCost;
        std::vector<std::vector<int>> _citiesNeighbors;
        std::vector<std::pair<float,float>> _minPairs;

        TSP() = default;
        virtual ~TSP() = default;

        int parse_inputs(int argc, char *argv[]) {

            if (argc != 3) {
                std::cerr << "Usage: tsp <cities file> <max value>" << std::endl;
                return 1;
            }

            // Read input from stdin
            std::string cities_file = argv[1];
            _maxValue = std::stof(argv[2]);

            // Initiate structures
            std::ifstream file(cities_file);
            std::string line;

            // Read cities file storing its information on <vector> cities
            bool firstLine = true;
            while (std::getline(file, line)) {
                std::stringstream line_ss(line);

                if (firstLine) {
                    line_ss >> _numberOfCities >> _numberOfRoads;
                    _roadsCost.resize(_numberOfCities, std::vector<float>(_numberOfCities, std::numeric_limits<float>::infinity()));
                    _citiesNeighbors.resize(_numberOfCities);
                    firstLine = false;
                    continue;
                }

                int cityO;
                int cityD;
                int cost;
                line_ss >> cityO >> cityD >> cost;
                _roadsCost[cityO][cityD] = cost;
                _roadsCost[cityD][cityO] = cost;
                _citiesNeighbors[cityO].push_back(cityD);
                _citiesNeighbors[cityD].push_back(cityO);
            }
            file.close();

            // Show the cost matrix
            std::cout << "----- INIT PHASE ------" << std::endl;
            std::cout << "Showing roads cost between each city" << std::endl;
            for (std::size_t i = 0; i < _roadsCost.size(); i++) {
                for (std::size_t j = 0; j < _roadsCost[i].size(); j++) {
                    std::cout << _roadsCost[i][j] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << "Showing cities neighbours list" << std::endl;
            for (std::size_t i = 0; i < _citiesNeighbors.size(); i++) {
                for (std::size_t j = 0; j < _citiesNeighbors[i].size(); j++) {
                    std::cout << _citiesNeighbors[i][j] << " ";
                }
                std::cout << std::endl;
            }
            return 0;
        }

        void print_result() {
            if (_bestTourCost > _maxValue) {
                std::cout << "NO SOLUTION" << std::endl;
            }
            else {
                std::cout << "Best tour cost: " << _bestTourCost << std::endl;
                std::cout << "Best tour:";
                for (const auto& city: _bestTour) {
                    std::cout << " " << city;
                }
                std::cout << std::endl;
            }
        }

        std::pair<float, float> findMinPairs(const int origin) {
            float min_one = std::numeric_limits<float>::max();
            float min_two = std::numeric_limits<float>::max();
            for (const int destiny: _citiesNeighbors[origin]) {
                if (_roadsCost[origin][destiny] < min_one) {
                    min_two = min_one;
                    min_one = _roadsCost[origin][destiny];
                } else if (_roadsCost[origin][destiny] < min_two && _roadsCost[origin][destiny] != min_one) {
                    min_two = _roadsCost[origin][destiny];
                }
            }
            return std::make_pair(min_one, min_two);
        }

        float computeInitialLowerBound() {
            float sum = 0;
            for (int origin = 0; origin < _numberOfCities; origin++) {
                std::pair<float, float> minPair = findMinPairs(origin);
                _minPairs.push_back(minPair);
                sum += minPair.first + minPair.second; 
            }
            std::cout << "Cities min1 and min2 info:" << std::endl;
            for (int i = 0; i < _numberOfCities; i++) {
                std::cout << " -> City " << i << " | min1: " << _minPairs[i].first << "; min2: " << _minPairs[i].second << std::endl;
            }
            return sum/2;
        }

        float newLowerBound(int origin, int destiny, float LB, float cost) {
            float cf = cost >= _minPairs[origin].second ? _minPairs[origin].second : _minPairs[origin].first;
            float ct = cost >= _minPairs[destiny].second ? _minPairs[destiny].second : _minPairs[destiny].first;
            //std::cout << "LB: " << LB << " cost: " << cost << " " << "cf: " << cf << " ct: " << ct << std::endl;
            return LB + cost - (cf + ct)/2;
        }

        void showTour(std::vector<int> tour) {
            for (const int city: tour) {
                std::cout << city << " ";
            }
        }

        void findSolution() {
            std::cout << "------ BRANCH AND BOUND ------" << std::endl;
            // Initialize Branch and Bound
            std::vector<int> tourOrigin;
            tourOrigin.push_back(0);
            float initial_lb = computeInitialLowerBound();
            VisitedCity newVisitedCity(tourOrigin, 0, initial_lb, 1, 0);
            PriorityQueue<VisitedCity, VisitedCity::CompareCityByLowerBound> queue;
            queue.push(newVisitedCity);
            
            VisitedCity city;
            std::vector<int> tour, newTour;
            float tourCost, bound, costUntilEnd, cost, newBound, newTourCost;
            int length, node;
            // Branch and Bound main loop
            while (!queue.empty()) {
                city = queue.pop();
                //std::cout << ".. Visiting new city in the context of tour: ";
                //showTour(city.getTour());
                //std::cout << "| Node: " << city.getNode() << "; Cost: " << city.getCost() <<
                //    "; Lower-bound: " << city.getLB() << "; Length: " << city.getLength() << std::endl;
                
                tour = city.getTour();
                tourCost = city.getCost();
                bound = city.getLB();
                length = city.getLength();
                node = city.getNode();

                if (bound >= _bestTourCost) {
                    return;
                }
                if (length+1 > _numberOfCities) {
                    costUntilEnd = tourCost + _roadsCost[node][0];
                    if (costUntilEnd < _bestTourCost) {
                        tour.push_back(0);
                        _bestTour = tour;
                        _bestTourCost = costUntilEnd;
                        std::cout << "!! Tour complete: ";
                        showTour(_bestTour);
                        std::cout << " -- Tour cost: " << _bestTourCost << std::endl;
                    }
                } else {
                    for (const int& destiny: _citiesNeighbors[node]) {
                        // If this destiny isn't yet part of the tour (std::find returns tour.end())
                        if (std::find(tour.begin(), tour.end(), destiny) == tour.end()) {
                            cost = _roadsCost[node][destiny];
                            newBound = newLowerBound(node, destiny, initial_lb, cost);
                            if (newBound > _bestTourCost) continue;
                            newTour = tour;
                            newTour.push_back(destiny);
                            newTourCost = tourCost + cost;
                            newVisitedCity = VisitedCity(newTour, newTourCost, newBound, length+1, destiny);
                            queue.push(newVisitedCity);
                        }
                    }
                }
            }
        }
};

int main(int argc, char *argv[]) {

    double exec_time;
    TSP tsp = TSP();
    if (tsp.parse_inputs(argc, argv)) return 1;
    exec_time = -omp_get_wtime();

    tsp.findSolution();
    
    exec_time += omp_get_wtime();

    std::cout << "----- RESULTS ------" << std::endl;
    fprintf(stderr, "Time of execution: %.6fs\n", exec_time);
    tsp.print_result();
    return 0;
}