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
#include <list>
#include <unistd.h>
#include <cstring>

// MPI comm tags
#define NODEREQ_T 1
#define BESTTOURCOST_T 2
#define ALL_GATHER_T 3
// process state flags
#define NO_NODES_TO_OFFER 0
#define OFFERING_NODE 1
// other flags
#define NO_NODE_REQUESTS -1

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
    double _bestTourCost = std::numeric_limits<double>::infinity(); // only updated when node finds a complete path
    double _bestLowerBound = std::numeric_limits<double>::infinity(); // updated every time it knows a lower cost
    std::shared_ptr<tour_t> _bestTour;
    std::vector<std::vector<double>> _roadsCost;
    std::vector<std::vector<int>> _citiesNeighbors;
    std::vector<std::pair<double, double>> _minPairs;
    PriorityQueue<std::shared_ptr<VisitedCity>, VisitedCity::CompareCityByLowerBound> queue;

    TSP() = default;
    virtual ~TSP() = default;

    int _numberOfProcesses, _pid;
    // Communication
    typedef struct {
        short destPid;
        short offer;
        double visitedCities;
        double cost;        
        double lowerBound;  
        int length;         
        short tour[1];      // MUST BE DECLARED AT LAST POSITION IN THE STRUCT!!
    } node_request;

    node_request* finalBestTourInfo;
    int NODE_REQUEST_SIZE;
    bool waiting = false;

    std::vector<int> recvdRequests;

    // node request
    std::vector<MPI_Request> sendNodeRequest, recvNodeRequest;
    int sendNodeRequest_buffer, *recvNodeRequest_buffer;

    // best tour cost
    std::vector<MPI_Request> sendBestTourCost, recvBestTourCost;
    std::vector<double> recvBestTourCost_buffer;
    double sendBestTourCost_buffer;

    // all gather
    char *recvAllGather_buffer;
    char *sendAllGather_buffer;

    void initMPICommunication(int pid, int number_of_processes) {
        _numberOfProcesses = number_of_processes;
        _pid = pid;
        NODE_REQUEST_SIZE = sizeof(node_request) + (_numberOfCities) * sizeof(short); // tour has size n_cities + 1

        // node requests
        sendNodeRequest.resize(number_of_processes);
        recvNodeRequest.resize(number_of_processes);
        recvNodeRequest_buffer = new int[number_of_processes];

        // best tour cost
        sendBestTourCost.resize(number_of_processes);
        recvBestTourCost.resize(number_of_processes);
        recvBestTourCost_buffer.resize(number_of_processes);
        finalBestTourInfo = (node_request*) new char[NODE_REQUEST_SIZE];
        finalBestTourInfo->cost = std::numeric_limits<double>::infinity();

        // all gather
        recvAllGather_buffer = new char[NODE_REQUEST_SIZE * number_of_processes * number_of_processes];
        sendAllGather_buffer = new char[NODE_REQUEST_SIZE * number_of_processes];
        recvdRequests.resize(number_of_processes, NO_NODE_REQUESTS);

        for (int i = 0; i < number_of_processes; i++) {
            // node request / return
            MPI_Send_init(&sendNodeRequest_buffer   , 1, MPI_INT, i, NODEREQ_T, MPI_COMM_WORLD, &sendNodeRequest[i]);
            MPI_Recv_init(&recvNodeRequest_buffer[i], 1, MPI_INT, i, NODEREQ_T, MPI_COMM_WORLD, &recvNodeRequest[i]);
            MPI_Start(&recvNodeRequest[i]);
            
            // best tour cost
            MPI_Send_init(&sendBestTourCost_buffer   , 1, MPI_DOUBLE, i, BESTTOURCOST_T, MPI_COMM_WORLD, &sendBestTourCost[i]);
            MPI_Recv_init(&recvBestTourCost_buffer[i], 1, MPI_DOUBLE, i, BESTTOURCOST_T, MPI_COMM_WORLD, &recvBestTourCost[i]);
            MPI_Start(&recvBestTourCost[i]);
        }
    }

    void destroyMPICommunication(int number_of_processes) {
        for (int i = 0; i < number_of_processes; i++) { 
            MPI_Request_free(&sendNodeRequest[i]);
            MPI_Request_free(&recvNodeRequest[i]);
            MPI_Request_free(&sendBestTourCost[i]);
            MPI_Request_free(&recvBestTourCost[i]);
        }
    }
 
    node_request* getNodeReqFromByteArray(char* array, int index) {
        return (node_request*)(array + index * NODE_REQUEST_SIZE);
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
        return 0;
    }

    void print_result() {
        std::cout << _pid << " printing results" << std::endl;
        // When bestourcost == maxvalue, since they are doubles,
        // they might not be exactly equal, thus we add 0.01 
        // (since the cost goes only up to 1 decimal place)
        if (finalBestTourInfo->cost > _maxValue + 0.01) {
            std::cout << "NO SOLUTION" << std::endl;
        }
        else {
            bool first = true;
            std::cout << std::fixed << std::setprecision(1) << finalBestTourInfo->cost << std::endl;
            for (int i = 0; i < _numberOfCities + 1; i++) {
                if (first) first = false;
                else std::cout << " ";
                std::cout << finalBestTourInfo->tour[i];
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
        return sum / 2;
    }

    double newLowerBound(int origin, int destiny, double LB, double cost) {
        double cf = cost >= _minPairs[origin].second ? _minPairs[origin].second : _minPairs[origin].first;
        double ct = cost >= _minPairs[destiny].second ? _minPairs[destiny].second : _minPairs[destiny].first;
        return LB + cost - (cf + ct) / 2;
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

    void addNodeToQueue(node_request * node_request_object) {
        std::shared_ptr<tour_t> newTour = nullptr;
        for (int node = 0; node < _numberOfCities; node++) {
            if (node_request_object->tour[node] == -1) break;
            newTour = extendTour(newTour, node_request_object->tour[node]);
        }
        queue.push(std::make_shared<VisitedCity>(newTour, node_request_object->cost, node_request_object->lowerBound, node_request_object->length));
    }

    void mpi_testRecvBestTourCost() {
        MPI_Status status;
        int finished;
        for (int i = 0; i < _numberOfProcesses; i++) {
            MPI_Test(&recvBestTourCost[i], &finished, &status);
            if (finished) {
                int senderPid = status.MPI_SOURCE;
                double newBestLowerBound = recvBestTourCost_buffer[senderPid];
                if (newBestLowerBound < _bestLowerBound) _bestLowerBound = newBestLowerBound;
                MPI_Start(&recvBestTourCost[senderPid]);
            }
        }
    }


    void mpi_testRecvNodeReq() {
        int didntReceiveAnything = NO_NODE_REQUESTS;
        MPI_Status status;
        int finished;
        for (int i = 0; i < _numberOfProcesses; i++) {
            MPI_Test(&recvNodeRequest[i], &finished, &status);
            if (finished) {
                int senderPid = status.MPI_SOURCE;
                recvdRequests[i] = senderPid;
                MPI_Start(&recvNodeRequest[senderPid]);
            }
        }
    }

    void copyTourIntoNodeRequest(node_request* node_request_object, std::shared_ptr<tour_t> tour, int length) {
        std::fill(node_request_object->tour, node_request_object->tour + _numberOfCities, (short) -1);
        for (int i = 0; i <= length - 1; i++) {
            node_request_object->tour[length - 1 - i] = static_cast<short>(tour->node);
            tour = tour->prev;
        }
    }

    void broadcastNodeReq() {
        std::vector<MPI_Status> statuses(_numberOfProcesses);
        for (int i = 0; i < _numberOfProcesses; i++) {
            //if (i == _pid) continue;
            MPI_Start(&sendNodeRequest[i]);
        }
        MPI_Waitall(_numberOfProcesses, sendNodeRequest.data(), statuses.data());
    }

    void broadcastBestTourCost() {
        std::vector<MPI_Status> statuses(_numberOfProcesses);
        for (int i = 0; i < _numberOfProcesses; i++) {
            if (i == _pid) continue;
            MPI_Start(&sendBestTourCost[i]);
        }
        MPI_Waitall(_numberOfProcesses, sendBestTourCost.data(), statuses.data());
    }


    bool mpi_shareStateAndNodes() {
        bool allAgreeOnTerminate = true;

        for (int i = 0 ; i < _numberOfProcesses ; i++) {
            // set nodes to send
            node_request* node_to_send = getNodeReqFromByteArray(sendAllGather_buffer, i);
            node_to_send->offer = NO_NODES_TO_OFFER;
            if (recvdRequests[i] != NO_NODE_REQUESTS && queue.size() > 3) {
                std::cout << _pid << " sending node to " << recvdRequests[i] << std::endl;
                node_to_send->destPid = recvdRequests[i];
                std::shared_ptr<VisitedCity> city = queue.pop();
                node_to_send->offer = OFFERING_NODE;
                node_to_send->cost = city->getCost();
                node_to_send->lowerBound = city->getLB();
                node_to_send->length = city->getLength();
                copyTourIntoNodeRequest(node_to_send, city->getTour(), city->getLength());
                recvdRequests[i] = NO_NODE_REQUESTS;
            }
        }

        // synchronize will all other processes
        MPI_Allgather(sendAllGather_buffer, NODE_REQUEST_SIZE * _numberOfProcesses, MPI_BYTE, 
        recvAllGather_buffer, NODE_REQUEST_SIZE * _numberOfProcesses , MPI_BYTE, MPI_COMM_WORLD);

        // check proposals if any. Else return a terminate flag
        for (int i = 0 ; i < _numberOfProcesses ; i++) { // all responses
            for (int j = 0 ; j < _numberOfProcesses ; j++) { // all nodes from each response
                node_request* recvdNode = getNodeReqFromByteArray(recvAllGather_buffer, _numberOfProcesses*i + j);
                if (recvdNode->offer == OFFERING_NODE) {
                    allAgreeOnTerminate = false;
                    if (recvdNode->destPid == _pid) {
                        std::cout << _pid << " received node from " << i << std::endl;
                        addNodeToQueue(recvdNode);
                        revoke_waiting();
                    }
                    recvdRequests[recvdNode->destPid] = NO_NODE_REQUESTS;
                }

            }
        }

        return allAgreeOnTerminate;
    }


    void mpi_shareBestToursAndDecideBest() {
        char* recvBuffer = new char[NODE_REQUEST_SIZE * _numberOfProcesses];
        // create send buffer
        finalBestTourInfo->cost = _bestTourCost;

        if (_bestTourCost < std::numeric_limits<double>::infinity()) { // might not have encoutered a full route
            copyTourIntoNodeRequest(finalBestTourInfo, _bestTour, _numberOfCities + 1); // best tour starts and ends at 0, thus n_cities + 1
        }
        // gather
        MPI_Gather((char*) finalBestTourInfo, NODE_REQUEST_SIZE, MPI_BYTE, recvBuffer,
            NODE_REQUEST_SIZE , MPI_BYTE, 0, MPI_COMM_WORLD);
        // choose the best from all received
        if (_pid == 0) {
            for (int i = 0; i < _numberOfProcesses; i++) {
                node_request* request = getNodeReqFromByteArray(recvBuffer, i);
                if ( request->cost < finalBestTourInfo->cost) {
                    memcpy(finalBestTourInfo, request, NODE_REQUEST_SIZE);
                }
            }
        }
    }

    void mark_as_waiting() {
        waiting = true;
    }

    void revoke_waiting() {
        waiting = false;
    }

    bool is_waiting() {
        return waiting;
    }

    void findSolution() {
        std::shared_ptr<tour_t> tour;
        if (_pid == 0) {
            queue.push(std::make_shared<VisitedCity>(extendTour(nullptr, 0), 0, computeInitialLowerBound(), 1));
        }
        else computeInitialLowerBound();
        // Branch and Bound main loop
        while (true) {
            mpi_testRecvBestTourCost();

            if (queue.empty()) {
                if (!is_waiting()) {
                    std::cout <<_pid << " issuing node req broadcast " << std::endl;
                    mark_as_waiting();
                    broadcastNodeReq();
                }
            } 
            else {
                std::shared_ptr<VisitedCity> city = queue.pop();
                tour = city->getTour();
                double tourCost = city->getCost();
                double bound = city->getLB();
                int length = city->getLength();
                int node = tour->node;

                if (bound >= _bestLowerBound) {
                    queue = PriorityQueue<std::shared_ptr<VisitedCity>, VisitedCity::CompareCityByLowerBound>();
                    continue;
                }
                
                if (length == _numberOfCities) {
                    double costUntilEnd = tourCost + _roadsCost[node][0];
                    if (costUntilEnd < _bestTourCost) {
                        if (!_bestTour) {
                            _bestTour = std::make_shared<tour_t>();
                            _bestTour->node = 0;
                        }
                        _bestTour->prev = tour;
                        _bestTourCost = costUntilEnd;
                        if (_bestTourCost < _bestLowerBound) {
                            _bestLowerBound = _bestTourCost;
                            sendBestTourCost_buffer = _bestTourCost;
                            broadcastBestTourCost();
                        }
                    }
                }
                else {  
                    for (const int& destiny : _citiesNeighbors[node]) {
                        if (hasNotVisitedCity(tour->visitedCities, destiny)) {
                            double cost = _roadsCost[node][destiny];
                            double newBound = newLowerBound(node, destiny, bound, cost);
                            if (newBound > _bestLowerBound) {
                                continue;
                            }
                            queue.push(std::make_shared<VisitedCity>(extendTour(tour, destiny), tourCost + cost, newBound, length + 1));
                        }
                    }
                }
            }

            mpi_testRecvNodeReq(); // get the ids of a processes that want nodes
            if (std::any_of(recvdRequests.begin(), recvdRequests.end(), [](int request) { return request != NO_NODE_REQUESTS; })) {
                std::cout << _pid << " got requests  "  << std::endl;
                bool allFinished = mpi_shareStateAndNodes();
                if (allFinished) break;
            }
        }
        mpi_shareBestToursAndDecideBest();
    };
};


int main(int argc, char* argv[]) {
    // sleep(10);
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
    MPI_Finalize();

    if (pid == 0) {
        std::cout << "----- RESULTS ------" << std::endl;
        fprintf(stderr, "Time of execution: %.6fs\n", elapsed_time);
        tsp.print_result();
    }

    return 0;
}
