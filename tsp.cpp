#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <omp.h>
#include <sstream>
#include <algorithm>
#include <limits>
#include <memory>
#include <stack>
#include "queue.hpp"

struct tour_t {
    int node;
    std::shared_ptr<tour_t> prev;
    int64_t visitedCities;
};

class VisitedCity {
public:
    std::shared_ptr<tour_t> _tour;
    double _cost;
    double _lowerBound;
    int _length;

    virtual ~VisitedCity() = default;
    VisitedCity() = default;

    struct CompareCityByLowerBound {
        bool operator()(const std::shared_ptr<VisitedCity>& a, const std::shared_ptr<VisitedCity>& b) {
            if (a->_lowerBound == b->_lowerBound) {
                return a->_tour->node > b->_tour->node;
            }
            return a->_lowerBound > b->_lowerBound;
        }
    };

    VisitedCity(std::shared_ptr<tour_t> tour, double cost, double lowerBound, int length) {
        _tour = tour;
        _cost = cost;
        _lowerBound = lowerBound;
        _length = length;
    }

    std::shared_ptr<tour_t> getTour() {
        return _tour;
    }
    double getCost() {
        return _cost;
    }
    double getLB() {
        return _lowerBound;
    }
    int getLength() {
        return _length;
    }
};

class TSP {
public:
    int _numberOfCities;
    int _numberOfRoads;
    double _maxValue;
    double _bestTourCost = std::numeric_limits<double>::infinity();
    std::shared_ptr<tour_t> _bestTour;
    std::vector<std::vector<double>> _roadsCost;
    std::vector<std::vector<int>> _citiesNeighbors;
    std::vector<std::pair<double, double>> _minPairs;

    std::vector<PriorityQueue<std::shared_ptr<VisitedCity>, VisitedCity::CompareCityByLowerBound>> _queues;
    std::vector<omp_lock_t*> _queueLocks;
    omp_lock_t* _bestTourLock = new omp_lock_t;
    int _numThreads = omp_get_num_threads();

    TSP() {
        for (int tid = 0; tid < _numThreads; tid++) {
            _queues.push_back(PriorityQueue<std::shared_ptr<VisitedCity>, VisitedCity::CompareCityByLowerBound>());
            _queueLocks.push_back(new omp_lock_t);
            omp_init_lock(_queueLocks[tid]);
        }
        omp_init_lock(_bestTourLock);
    }

    ~TSP() {
        for (int tid = 0; tid < _numThreads; tid++) {
            omp_destroy_lock(_queueLocks[tid]);
        }
        omp_destroy_lock(_bestTourLock);
    }

    int parse_inputs(int argc, char* argv[]) {

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
                _roadsCost.resize(_numberOfCities, std::vector<double>(_numberOfCities, std::numeric_limits<double>::infinity()));
                _citiesNeighbors.resize(_numberOfCities);
                firstLine = false;
                continue;
            }

            int cityO;
            int cityD;
            double cost;
            line_ss >> cityO >> cityD >> cost;
            _roadsCost[cityO][cityD] = cost;
            _roadsCost[cityD][cityO] = cost;
            _citiesNeighbors[cityO].push_back(cityD);
            _citiesNeighbors[cityD].push_back(cityO);
        }
        file.close();

        // Show the cost matrix
        //std::cout << "----- INIT PHASE ------" << std::endl;
        //std::cout << "Showing roads cost between each city" << std::endl;
        //for (std::size_t i = 0; i < _roadsCost.size(); i++) {
        //    for (std::size_t j = 0; j < _roadsCost[i].size(); j++) {
        //        std::cout << _roadsCost[i][j] << " ";
        //    }
        //    std::cout << std::endl;
        //}
        //std::cout << "Showing cities neighbours list" << std::endl;
        //for (std::size_t i = 0; i < _citiesNeighbors.size(); i++) {
        //    for (std::size_t j = 0; j < _citiesNeighbors[i].size(); j++) {
        //        std::cout << _citiesNeighbors[i][j] << " ";
        //    }
        //    std::cout << std::endl;
        //}
        return 0;
    }

    void print_result() {
        // When bestourcost == maxvalue, since they are doubles,
        // they might not be exactly equal, thus we add 0.01 
        // (since the cost goes only up to 1 decimal place)
        if (_bestTourCost > _maxValue + 0.01) {
            std::cout << "NO SOLUTION" << std::endl;
        }
        else {
            bool first = true;
            std::stack<int> bestTour;
            while (_bestTour != nullptr) {
                bestTour.push(_bestTour->node);
                _bestTour = _bestTour->prev;
            }
            std::cout << std::fixed << std::setprecision(1) << _bestTourCost << std::endl;
            while (bestTour.size() > 0) {
                if (first) first = false;
                else std::cout << " ";
                std::cout << bestTour.top();
                bestTour.pop();
            }
            std::cout << std::endl;
        }
    }

    std::pair<double, double> findMinPairs(const int origin) {
        double min_one = std::numeric_limits<double>::max();
        double min_two = std::numeric_limits<double>::max();
        for (const int destiny : _citiesNeighbors[origin]) {
            if (_roadsCost[origin][destiny] < min_one) {
                min_two = min_one;
                min_one = _roadsCost[origin][destiny];
            }
            else if (_roadsCost[origin][destiny] < min_two) {
                min_two = _roadsCost[origin][destiny];
            }
        }
        return std::make_pair(min_one, min_two);
    }

    double computeInitialLowerBound() {
        double sum = 0;
        for (int origin = 0; origin < _numberOfCities; origin++) {
            std::pair<double, double> minPair = findMinPairs(origin);
            _minPairs.push_back(minPair);
            sum += minPair.first + minPair.second;
        }
        //std::cout << "Cities min1 and min2 info:" << std::endl;
        //for (int i = 0; i < _numberOfCities; i++) {
        //    std::cout << " -> City " << i << " | min1: " << _minPairs[i].first << "; min2: " << _minPairs[i].second << std::endl;
        //}
        return sum / 2;
    }

    double newLowerBound(int origin, int destiny, double LB, double cost) {
        double cf = cost >= _minPairs[origin].second ? _minPairs[origin].second : _minPairs[origin].first;
        double ct = cost >= _minPairs[destiny].second ? _minPairs[destiny].second : _minPairs[destiny].first;
        //std::cout << "LB: " << LB << " cost: " << cost << " " << "cf: " << cf << " ct: " << ct << std::endl;
        return LB + cost - (cf + ct) / 2;
    }

    void showTour(std::shared_ptr<tour_t> tour) {
        std::stack<int> tourStack;
        while (tour != nullptr) {
            tourStack.push(tour->node);
            tour = tour->prev;
        }
        while (tourStack.size() > 0) {
            std::cout << " " << tourStack.top();
            tourStack.pop();
        }
    }

    bool tourHasNode(std::shared_ptr<tour_t> tour, int node) {
        while (tour != nullptr) {
            if (tour->node == node) {
                return true;
            }
            tour = tour->prev;
        }
        return false;
    }

    int64_t visit(int64_t visitedCities, int destiny) {
        return visitedCities |= static_cast<int64_t>(1 << destiny);
    }

    bool hasNotVisitedCity(int64_t visitedCities, int city) {
        return (visitedCities &= static_cast<int64_t>(1 << city)) == 0;
    }

    std::shared_ptr<tour_t> extendTour(std::shared_ptr<tour_t> tour, int destiny) {
        std::shared_ptr<tour_t> newTour = std::make_shared<tour_t>();
        newTour->prev = tour;
        newTour->node = destiny;
        newTour->visitedCities = visit(tour == nullptr ? 0 : tour->visitedCities, destiny);
        return newTour;
    }

    void printBits(int num) {
        // Determine the number of bits needed to represent the integer
        int numBits = sizeof(num) * 8;

        // Create a mask to extract each bit from the integer
        unsigned int mask = 1u << (numBits - 1);

        // Loop through each bit and print it as 0 or 1
        for (int i = 0; i < numBits; i++) {
            int bit = (num & mask) ? 1 : 0;
            std::cout << bit;
            mask >>= 1;
        }

        std::cout << std::endl;
    }

    int existElementsInQueues(std::vector<PriorityQueue<std::shared_ptr<VisitedCity>, VisitedCity::CompareCityByLowerBound>> queues) {
        for (int thread_id = omp_get_thread_num(); thread_id < omp_get_thread_num() + _numThreads; thread_id++) {
            if (queues.size() != 0) return thread_id % _numThreads;
        }
        return -1;
    }


    void findSolution() {
        //std::cout << "------ BRANCH AND BOUND ------" << std::endl;
        // Initialize Branch and Bound
# pragma omp master 
        {
            std::shared_ptr<tour_t> tour = extendTour(nullptr, 0);
        }
        double initial_lb = computeInitialLowerBound();
# pragma omp barrier 

        int idx;
        // Branch and Bound main loop
        while ((idx = existElementsInQueues(_queues)) != -1) {

            omp_set_lock(_queueLocks[idx]);
            std::shared_ptr<VisitedCity> city = _queues[idx].pop();
            omp_unset_lock(_queueLocks[idx]);
            std::shared_ptr<tour_t> currentTour = city->getTour();
            double tourCost = city->getCost();
            double bound = city->getLB();
            int length = city->getLength();
            int node = currentTour->node;

            //std::cout << ".. Visiting new city in the context of tour: (";
            //showTour(city->getTour());
            //std::cout << " ) | Node: " << node << "; Cost: " << tourCost <<
            //    "; Lower-bound: " << bound << "; Length: " << length << "; Visited cities: ";
            //printBits(currentTour->visitedCities);
            //std::cout << std::endl;

            if (bound >= _bestTourCost) {
                return;
            }

            if (length == _numberOfCities) {
                double costUntilEnd = tourCost + _roadsCost[node][0];

                omp_set_lock(_bestTourLock);
                if (costUntilEnd < _bestTourCost) {
                    if (!_bestTour) {
                        _bestTour = std::make_shared<tour_t>();
                        _bestTour->node = 0;
                    }
                    _bestTour->prev = currentTour;
                    _bestTourCost = costUntilEnd;
                    //std::cout << "!! Tour complete: ";
                    //showTour(_bestTour);
                    //std::cout << " -- Tour cost: " << _bestTourCost << std::endl;
                }
                omp_unset_lock(_bestTourLock);
            }
            else {
                for (const int& destiny : _citiesNeighbors[node]) {
                    if (hasNotVisitedCity(currentTour->visitedCities, destiny)) {
                        double cost = _roadsCost[node][destiny];
                        double newBound = newLowerBound(node, destiny, bound, cost);
                        if (newBound > _bestTourCost) continue;
                        std::shared_ptr<tour_t> newTour = extendTour(currentTour, destiny);
                        omp_set_lock(_queueLocks[omp_get_thread_num()]);
                        _queues[omp_get_thread_num()].push(std::make_shared<VisitedCity>(newTour, tourCost + cost, newBound, length + 1));
                        omp_unset_lock(_queueLocks[omp_get_thread_num()]);
                    }
                }
            }

        }
    }
};

int main(int argc, char* argv[]) {

    double exec_time;
    TSP tsp = TSP();
    if (tsp.parse_inputs(argc, argv)) return 1;
    exec_time = -omp_get_wtime();

#pragma omp parallel num_threads(4) 
    {
        tsp.findSolution();
    }

    exec_time += omp_get_wtime();

    //std::cout << "----- RESULTS ------" << std::endl;
    //fprintf(stderr, "Time of execution: %.6fs\n", exec_time);
    tsp.print_result();
    return 0;
}
