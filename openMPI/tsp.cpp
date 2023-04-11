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

// MPI comm tags
#define NODEREQ_T 1
# define NODERET_T 2
# define BESTTOURCOST_T 3
// process state flags
#define NOT_TERMINATED 0
#define TERMINATED 1

#define NODE_LIMIT 3

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
    PriorityQueue<std::shared_ptr<VisitedCity>, VisitedCity::CompareCityByLowerBound> queue;

    TSP() = default;
    virtual ~TSP() = default;

    int _numberOfProcesses, _pid;
    // Communication
    typedef struct {
        short terminated;
        double visitedCities;
        double cost;        
        double lowerBound;  // needed ?
        int length;         // needed ?
        short tour[1];      // MUST BE DECLARED AT LAST POSITION IN THE STRUCT!!
    } node_request;

    node_request* finalBestTourInfo;
    std::vector<short> processes_state;
    int NODE_REQUEST_SIZE;
    

    // Mpi Request structures
    std::vector<MPI_Request> sendNodeRequest, sendNodeReturn, sendBestTourCost;
    std::vector<MPI_Request> recvNodeRequest, recvNodeReturn, recvBestTourCost;
    // Buffers
    std::vector<node_request*> sendNodeBuffer, recvNodeBuffer;
    std::vector<double> recv_bestTourCostBuffer;
    double bestTourCostBuffer;
    char* sendNodeRequestBuffer = new char[1], *recvNodeRequestBuffer = new char[1];


    void initMPICommunication(int pid, int number_of_processes) {
        _numberOfProcesses = number_of_processes;
        _pid = pid;
        processes_state.resize(number_of_processes, NOT_TERMINATED);
        NODE_REQUEST_SIZE = sizeof(node_request) + (_numberOfCities) * sizeof(short); // tour has size n_cities + 1
        finalBestTourInfo = (node_request*) new char[NODE_REQUEST_SIZE]; // test
        finalBestTourInfo->cost = std::numeric_limits<double>::infinity();
        sendNodeBuffer.resize(number_of_processes, (node_request*) new char[NODE_REQUEST_SIZE]);
        recvNodeBuffer.resize(number_of_processes, (node_request*) new char[NODE_REQUEST_SIZE]);
        sendNodeRequest.resize(number_of_processes);
        sendNodeReturn.resize(number_of_processes);
        sendBestTourCost.resize(number_of_processes);
        recvNodeRequest.resize(number_of_processes);
        recvNodeReturn.resize(number_of_processes);
        recvBestTourCost.resize(number_of_processes);
        recv_bestTourCostBuffer.resize(number_of_processes);
        for (int i = 0; i < number_of_processes; i++) {
            MPI_Send_init(sendNodeRequestBuffer, 1, MPI_BYTE, i, NODEREQ_T, MPI_COMM_WORLD, &sendNodeRequest[i]);

            MPI_Send_init(sendNodeBuffer[i], NODE_REQUEST_SIZE, MPI_BYTE, i, NODERET_T, MPI_COMM_WORLD, &sendNodeReturn[i]);
            MPI_Send_init(&bestTourCostBuffer, 1, MPI_DOUBLE, i, BESTTOURCOST_T, MPI_COMM_WORLD, &sendBestTourCost[i]);
            
            MPI_Recv_init(recvNodeRequestBuffer, 1, MPI_BYTE, i, NODEREQ_T, MPI_COMM_WORLD, &recvNodeRequest[i]);
            MPI_Start(&recvNodeRequest[i]);
            
            MPI_Recv_init(recvNodeBuffer[i], NODE_REQUEST_SIZE, MPI_BYTE, i, NODERET_T, MPI_COMM_WORLD, &recvNodeReturn[i]);
            MPI_Start(&recvNodeReturn[i]);
            
            MPI_Recv_init(&recv_bestTourCostBuffer[i], 1, MPI_DOUBLE, i, BESTTOURCOST_T, MPI_COMM_WORLD, &recvBestTourCost[i]);
            MPI_Start(&recvBestTourCost[i]);
            
        }
    }

    void destroyMPICommunication(int number_of_processes) {
        for (int i = 0; i < number_of_processes; i++) { // maybe after each request?
            MPI_Request_free(&sendNodeReturn[i]);
            MPI_Request_free(&sendNodeRequest[i]);
            MPI_Request_free(&sendBestTourCost[i]);
            MPI_Request_free(&recvNodeRequest[i]);
            MPI_Request_free(&recvNodeReturn[i]);
            MPI_Request_free(&recvBestTourCost[i]);
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
        std::cout << _pid << " printing results" << std::endl;
        // When bestourcost == maxvalue, since they are doubles,
        // they might not be exactly equal, thus we add 0.01 
        // (since the cost goes only up to 1 decimal place)
        if (_bestTourCost > _maxValue + 0.01) {
            std::cout << "NO SOLUTION" << std::endl;
        }
        else {
            bool first = true;
            std::cout << std::fixed << std::setprecision(1) << _bestTourCost << std::endl;
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
            //std::cout << " " << tourStack.top();
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

    void mark_as_terminated(int id) {
        processes_state[id] = TERMINATED;
    }

    void revoke_terminated(int id) {
        processes_state[id] = NOT_TERMINATED;
    }

    int get_process_state() {
        return processes_state[_pid];
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
        revoke_terminated(_pid);
    }

    void mpi_testRecvNodeReturn() {
        std::cout << _pid << " test Recv Nodes " << std::endl;
        MPI_Status status;
        int finished;
        for (int i = 0; i < _numberOfProcesses; i++) {
            MPI_Test(&recvNodeReturn[i], &finished, &status);
            if (finished) {
                std::cout << _pid << " received request from " << i << " with MPI-SOURCE = " << status.MPI_SOURCE  << std::endl;
                int senderPid = status.MPI_SOURCE;
                node_request* node_request_object = recvNodeBuffer[senderPid];
                if (node_request_object->terminated) {
                    mark_as_terminated(senderPid);
                }
                else {
                    addNodeToQueue(node_request_object);
                    revoke_terminated(senderPid);
                }
                MPI_Start(&recvNodeReturn[senderPid]);
            }
        }

    }


    void mpi_testRecvBestTourCost() {
        std::cout << _pid << " test Recv bestTourCost " << std::endl;
        MPI_Status status;
        int finished;
        for (int i = 0; i < _numberOfProcesses; i++) {
            MPI_Test(&recvBestTourCost[i], &finished, &status);
            if (finished) {
                int senderPid = status.MPI_SOURCE;
                double newBestTourCost = recv_bestTourCostBuffer[senderPid];
                if (_bestTourCost < newBestTourCost) {
                    _bestTourCost = newBestTourCost;
                }
                std::cout <<_pid << " Receive Cost " << newBestTourCost << " From " << senderPid  << " " << std::endl;
                MPI_Start(&recvBestTourCost[senderPid]);
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

    void handleRecvNodeRequest(int destPid) {
        MPI_Status status;
        int queueSize = queue.size();
        if (queueSize > 0 && queueSize < NODE_LIMIT) {
            std::cout << _pid << " dont have enough nodes to send to " << destPid << std::endl;
            return;
        }
        std::shared_ptr<node_request> node_request_object = std::make_shared<node_request>();
        if (queueSize == 0) node_request_object->terminated = TERMINATED;
        else {
            std::shared_ptr<VisitedCity> city = queue.pop();
            node_request_object->terminated = NOT_TERMINATED;
            node_request_object->cost = city->getCost();
            std::cout << _pid << " city cost = " << city->getCost() << std::endl;
            std::cout << _pid << " node cost = " << node_request_object->cost << std::endl;
            node_request_object->lowerBound = city->getLB();
            copyTourIntoNodeRequest(node_request_object.get(), city->getTour(), city->getLength());
        }
        sendNodeBuffer[destPid] = node_request_object.get();
        node_request_object = nullptr;
        MPI_Wait(&sendNodeReturn[destPid], &status); // must wait before issuing next send
        MPI_Start(&sendNodeReturn[destPid]);
    }

    void mpi_testRecvNodeRequest() {
        std::cout << _pid << " testing recv node requests" << std::endl;
        MPI_Status status;
        int finished;
        for (int i = 0; i < _numberOfProcesses; i++) {
            MPI_Test(&recvNodeRequest[i], &finished, &status);
            if (finished) {
                std::cout << _pid << " received a node request from " << i << " MPI_SOURCE= " << status.MPI_SOURCE << std::endl;
                int senderPid = status.MPI_SOURCE;
                handleRecvNodeRequest(senderPid);
                MPI_Start(&recvNodeRequest[senderPid]);
            }
        }
    }

    void broadcast(std::vector<MPI_Request>& requestObj) {
        //std::cout << _pid << " broadcast start " << std::endl;
        std::vector<MPI_Status> statuses(_numberOfProcesses);
        for (int i = 0; i < _numberOfProcesses; i++) {
            // attention lack of & before requestObj on start because of & in parameters
            MPI_Start(&requestObj[i]);
            std::cout << _pid << " sending to " << i << std::endl;

        }
        int err_code = MPI_Waitall(_numberOfProcesses, requestObj.data(), statuses.data());
        if (err_code != MPI_SUCCESS) {
            std::cout << _pid << " MPI_Waitall error " << std::endl;
        }
        //std::cout << _pid << " broadcast end" << std::endl;
    }

    void mpi_shareBestToursAndDecideBest() {
        std::cout << _pid << " START gather results" << std::endl;
        char* recvBuffer = (char*) malloc(NODE_REQUEST_SIZE * _numberOfProcesses);
        // create send buffer
        finalBestTourInfo->cost = _bestTourCost;
        std::cout << _pid << " cost = " << finalBestTourInfo->cost << std::endl;

        if (_bestTourCost < std::numeric_limits<double>::infinity()) { // might not have encoutered a full route
            std::cout << _pid << " copying tour" << std::endl;
            copyTourIntoNodeRequest(finalBestTourInfo, _bestTour, _numberOfCities + 1); // best tour starts and ends at 0, thus n_cities + 1
        }
        // gather
        std::cout << _pid << " gather results" << std::endl;
        MPI_Gather(finalBestTourInfo, NODE_REQUEST_SIZE, MPI_BYTE, recvBuffer,
            NODE_REQUEST_SIZE, MPI_BYTE, 0, MPI_COMM_WORLD);
        std::cout << _pid << " after gather results" << std::endl;
        // choose the best from all received
        for (int i = 0; i < _numberOfProcesses; i++) {
            node_request* request = (node_request*)(recvBuffer + i * NODE_REQUEST_SIZE);
            std::cout << _pid << " received cost of " << request->cost << " from " << i << std::endl;
            if ( request->cost < finalBestTourInfo->cost) {
                std::cout << _pid << " choose from " << i << std::endl;
                finalBestTourInfo = request;
            }
        }
        std::cout << _pid << " gather results end" << std::endl;
    }


    void findSolution() {
        std::shared_ptr<tour_t> tour;
        if (_pid == 0) {
            queue.push(std::make_shared<VisitedCity>(extendTour(nullptr, 0), 0, computeInitialLowerBound(), 1));
        }

        // Branch and Bound main loop
        while (std::any_of(processes_state.begin(), processes_state.end(), [](short state) { return state == NOT_TERMINATED; })) {
            if (queue.empty()) {
                std::cout << _pid << " has empty queue" << std::endl;
                mark_as_terminated(_pid);
                std::cout << _pid << " requesting nodes" << std::endl;
                mpi_testRecvNodeRequest();
                broadcast(sendNodeRequest);
                mpi_testRecvNodeReturn();
                mpi_testRecvBestTourCost();

                // Continue (while receiving messages and updating process_state)
            } else {

                mpi_testRecvNodeRequest();
                mpi_testRecvBestTourCost();

                std::shared_ptr<VisitedCity> city = queue.pop();
                tour = city->getTour();
                double tourCost = city->getCost();
                double bound = city->getLB();
                int length = city->getLength();
                int node = tour->node;
                std::cout << _pid << " has " << queue.size() << " nodes in queue" << std::endl;


                //std::cout << ".. Visiting new city in the context of tour: (";
                //showTour(city->getTour());
                //std::cout << " ) | Node: " << node << "; Cost: " << tourCost <<
                //    "; Lower-bound: " << bound << "; Length: " << length << "; Visited cities: ";
                //printBits(tour->visitedCities);
                //std::cout << std::endl;
                if (bound >= _bestTourCost) {
                    queue = PriorityQueue<std::shared_ptr<VisitedCity>, VisitedCity::CompareCityByLowerBound>();
                    // no need to set buffer (nullptr)
                    std::cout << _pid << " bound >= bestTourCost" << std::endl;
                    broadcast(sendNodeRequest);
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
                        bestTourCostBuffer = _bestTourCost;
                        std::cout << _pid << " sending bestTourCost = " << bestTourCostBuffer << std::endl;
                        broadcast(sendBestTourCost);
                        //showTour(_bestTour);
                        //std::cout << " -- Tour cost: " << _bestTourCost << std::endl;
                    }
                }
                else {
                    for (const int& destiny : _citiesNeighbors[node]) {
                        if (hasNotVisitedCity(tour->visitedCities, destiny)) {
                            double cost = _roadsCost[node][destiny];
                            double newBound = newLowerBound(node, destiny, bound, cost);
                            if (newBound > _bestTourCost) continue;
                            queue.push(std::make_shared<VisitedCity>(extendTour(tour, destiny), tourCost + cost, newBound, length + 1));
                        }
                    }
                }
                std::cout << _pid << " testing node req on ELSE" << std::endl;
            }
        }
        std::cout << _pid << " finished findSolution" << std::endl;
        mpi_shareBestToursAndDecideBest();
    };
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
    std::cout << pid << " exited findSol" << std::endl;
    elapsed_time += MPI_Wtime();

    tsp.destroyMPICommunication(number_of_processes);

    if (pid == 0) {
        std::cout << "----- RESULTS ------" << std::endl;
        fprintf(stderr, "Time of execution: %.6fs\n", elapsed_time);
        tsp.print_result();
    }

    // Attention! please 
    MPI_Finalize();

    return 0;
}
