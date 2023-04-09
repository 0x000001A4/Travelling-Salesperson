#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <mpi.h>
#include <sstream>
#include <algorithm>
#include <limits>
#include <memory>
#include <stack>
#include "queue.hpp"

// MPI comm tags
#define NODEREQ_T 420
// process state flags
#define TERMINATED 1
#define NOT_TERMINATED_YET 0

#define MINIMUM_NODE_THREASHOLD_REQ 1

struct tour_t {
    int node;
    std::shared_ptr<tour_t> prev;
    double visitedCities;
};

class VisitedCity {
public:
    std::shared_ptr<tour_t> _tour;
    double _cost;
    double _lowerBound;
    int _length;

    VisitedCity() = default;

    virtual ~VisitedCity() = default;

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

    TSP() = default;
    virtual ~TSP() = default;

    // Communication
    typedef struct {
        short terminated;
        double visitedCities;
        short tour[1];
    } node_request;

    std::vector<short> processes_state;
    std::vector<node_request*> sendReqsBuffer, recvReqsBuffer;
    std::vector<MPI_Request> send_req, recv_req;
    int NODE_REQUEST_SIZE;
    int pid;

    void initMPICommunication(int pid, int number_of_processes) {
        processes_state.resize(number_of_processes, NOT_TERMINATED_YET);
        NODE_REQUEST_SIZE = sizeof(node_request) + (_numberOfCities - 1) * sizeof(int);
        sendReqsBuffer.resize(number_of_processes, (node_request*) new char[NODE_REQUEST_SIZE]);
        recvReqsBuffer.resize(number_of_processes, (node_request*) new char[NODE_REQUEST_SIZE]);
        send_req.resize(number_of_processes);
        recv_req.resize(number_of_processes);
        for (int i = 0; i < number_of_processes; i++) {
            MPI_Send_init(sendReqsBuffer[i], NODE_REQUEST_SIZE, MPI_BYTE, i, NODEREQ_T, MPI_COMM_WORLD, &send_req[i]);
            MPI_Recv_init(recvReqsBuffer[i], NODE_REQUEST_SIZE, MPI_BYTE, i, NODEREQ_T, MPI_COMM_WORLD, &recv_req[i]);
        }
    }

    void destroyMPICommunication(int number_of_processes) {
        for (int i = 0; i < number_of_processes; i++) {
            MPI_Request_free(&recv_req[i]);
            MPI_Request_free(&send_req[i]);
        }
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

    int visit(int visitedCities, int destiny) {
        return visitedCities |= (1 << destiny);
    }

    bool hasNotVisitedCity(int visitedCities, int city) {
        return (visitedCities &= (1 << city)) == 0;
    }

    std::shared_ptr<tour_t> extendTour(std::shared_ptr<tour_t> tour, int destiny) {
        std::shared_ptr<tour_t> newTour = std::make_shared<tour_t>();
        newTour->prev = tour;
        newTour->node = destiny;
        newTour->visitedCities = visit(tour == nullptr ? 0 : tour->visitedCities, destiny);
        return newTour;
    }

    void findSolution() {
        //std::cout << "------ BRANCH AND BOUND ------" << std::endl;
        // Initialize Branch and Bound
        std::shared_ptr<tour_t> tour = extendTour(nullptr, 0);
        double initial_lb = computeInitialLowerBound();
        PriorityQueue<std::shared_ptr<VisitedCity>, VisitedCity::CompareCityByLowerBound> queue;
        queue.push(std::make_shared<VisitedCity>(tour, 0, initial_lb, 1));

        std::shared_ptr<VisitedCity> city;
        double tourCost, bound, costUntilEnd, cost, newBound, newTourCost;
        int length, node;
        // Branch and Bound main loop
        while (std::any_of(processes_state.begin(), processes_state.end(), [](short state){ 
            return state == NOT_TERMINATED_YET; // termination condition
        })){
            if (!queue.empty()) {
                if (queue.size() == MINIMUM_NODE_THREASHOLD_REQ) {
                    // perform request nodes
                }
                city = queue.pop();
                tour = city->getTour();
                tourCost = city->getCost();
                bound = city->getLB();
                length = city->getLength();
                node = tour->node;

                //std::cout << ".. Visiting new city in the context of tour: (";
                //showTour(city->getTour());
                //std::cout << " ) | Node: " << node << "; Cost: " << tourCost <<
                //    "; Lower-bound: " << bound << "; Length: " << length << "; Visited cities: ";
                //printBits(tour->visitedCities);
                //std::cout << std::endl;

                if (bound >= _bestTourCost) {
                    return;
                }

                
                if (length == _numberOfCities) {
                    costUntilEnd = tourCost + _roadsCost[node][0];
                    if (costUntilEnd < _bestTourCost) {
                        if (!_bestTour) {
                            _bestTour = std::make_shared<tour_t>();
                            _bestTour->node = 0;
                        }
                        _bestTour->prev = tour;
                        _bestTourCost = costUntilEnd;
                        //std::cout << "!! Tour complete: ";
                        //showTour(_bestTour);
                        //std::cout << " -- Tour cost: " << _bestTourCost << std::endl;
                    }
                }
                else {
                    for (const int& destiny : _citiesNeighbors[node]) {
                        if (hasNotVisitedCity(tour->visitedCities, destiny)) {
                            cost = _roadsCost[node][destiny];
                            newBound = newLowerBound(node, destiny, bound, cost);
                            if (newBound > _bestTourCost) continue;
                            queue.push(std::make_shared<VisitedCity>(extendTour(tour, destiny), tourCost + cost, newBound, length + 1));
                        }
                    }
                }
                
            }
            else {
                mark_as_terminated(pid);
                // wait for requests
            }
        }
    }

    void mark_as_terminated(int id) {
        processes_state[id] = TERMINATED;
    }

    void revoke_terminated(int id) {
        processes_state[id] = NOT_TERMINATED_YET;
    }
};


int main(int argc, char* argv[]) {

    double elapsed_time;

    int pid, number_of_processes;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    TSP tsp = TSP();
    if (tsp.parse_inputs(argc, argv)) return 1;
    tsp.initMPICommunication(pid, number_of_processes);

    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();
    tsp.findSolution();
    elapsed_time += MPI_Wtime();

    tsp.destroyMPICommunication(number_of_processes);

    std::cout << "----- RESULTS ------" << std::endl;
    fprintf(stderr, "Time of execution: %.6fs\n", elapsed_time);
    tsp.print_result();

    // Attention! please 
    MPI_Finalize();

    return 0;
}
