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
#define NODERET_T 2
#define NODEREQRET_T 3
#define BESTTOURCOST_T 4
#define TERMINATION_T 5
#define NODERET_WITHNODE_T 6
#define NODERET_WITHOUTNODE_T 7
// process state flags
#define NOT_TERMINATED 0
#define TERMINATED 1
#define READY_TO_TERMINATE 2

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
        short type;
        int timestamp;
        short terminated;
        double visitedCities;
        double cost;        
        double lowerBound;  // needed ?
        int length;         // needed ?
        short tour[1];      // MUST BE DECLARED AT LAST POSITION IN THE STRUCT!!
    } node_request;

    node_request* finalBestTourInfo;
    std::vector<short> processes_state;
    std::vector<int> recvTimestamps;
    std::vector<int> sendTimestamps;
    // keep track of all sent requests to wait for all received node messages
    std::vector<int> requestCount;
    std::vector<int> responseCount;
    int NODE_REQUEST_SIZE;
    

    // Mpi Request structures
    std::vector<MPI_Request> sendNodeRequest, sendNodeReturn, sendBestTourCost;
    std::vector<MPI_Request> recvNodeRequestReturn, recvBestTourCost;

    std::vector<MPI_Request> sendTerminationReq;
    std::vector<MPI_Request> sendTerminationReturn, recvTerminationReqReturn;

    // Buffers
    char *sendNodeBuffer_m, *recvNodeReqRetBuffer;
    node_request* sendNodeRequestBuffer;

    int* recvTerminationReqReturn_buffer;
    int  sendTerminationReq_buffer, sendTerminationReturn_buffer, recvTerminationReq_buffer;

    std::vector<double> recv_bestTourCostBuffer;
    double bestTourCostBuffer;

    node_request* getNodeReqFromByteArray(char* array, int index) {
        return (node_request*)(array + index * NODE_REQUEST_SIZE);
    }

    void initMPICommunication(int pid, int number_of_processes) {
        _numberOfProcesses = number_of_processes;
        _pid = pid;
        processes_state.resize(number_of_processes, NOT_TERMINATED);
        NODE_REQUEST_SIZE = sizeof(node_request) + (_numberOfCities) * sizeof(short); // tour has size n_cities + 1

        sendTimestamps.resize(number_of_processes, 0);
        recvTimestamps.resize(number_of_processes, 0);

        // keep track of expected responses
        requestCount.resize(number_of_processes, 0);
        responseCount.resize(number_of_processes, 0);

        // node requests / returns
        sendNodeRequest.resize(number_of_processes);
        sendNodeReturn.resize(number_of_processes);
        recvNodeRequestReturn.resize(number_of_processes);
        sendNodeBuffer_m = new char[NODE_REQUEST_SIZE * _numberOfProcesses];
        recvNodeReqRetBuffer = new char[NODE_REQUEST_SIZE * _numberOfProcesses];
        sendNodeRequestBuffer = (node_request*) new char[NODE_REQUEST_SIZE];

        // best tour cost
        sendBestTourCost.resize(number_of_processes);
        recvBestTourCost.resize(number_of_processes);
        recv_bestTourCostBuffer.resize(number_of_processes);
        finalBestTourInfo = (node_request*) new char[NODE_REQUEST_SIZE]; // test
        finalBestTourInfo->cost = std::numeric_limits<double>::infinity();

        // termination
        sendTerminationReq.resize(number_of_processes);
        recvTerminationReqReturn.resize(number_of_processes);
        sendTerminationReturn.resize(number_of_processes);
        recvTerminationReqReturn_buffer = (int*) calloc(number_of_processes, sizeof(int));

        for (int i = 0; i < number_of_processes; i++) {
            // node request / return
            MPI_Send_init(sendNodeRequestBuffer                           , NODE_REQUEST_SIZE, MPI_BYTE, i, NODEREQRET_T, MPI_COMM_WORLD, &sendNodeRequest[i]);
            MPI_Send_init(getNodeReqFromByteArray(sendNodeBuffer_m, i)    , NODE_REQUEST_SIZE, MPI_BYTE, i, NODEREQRET_T, MPI_COMM_WORLD, &sendNodeReturn[i]);
            MPI_Recv_init(getNodeReqFromByteArray(recvNodeReqRetBuffer, i), NODE_REQUEST_SIZE, MPI_BYTE, i, NODEREQRET_T, MPI_COMM_WORLD, &recvNodeRequestReturn[i]);
            MPI_Start(&recvNodeRequestReturn[i]);
            
            // best tour cost
            MPI_Send_init(&bestTourCostBuffer        , 1, MPI_DOUBLE, i, BESTTOURCOST_T, MPI_COMM_WORLD, &sendBestTourCost[i]);
            MPI_Recv_init(&recv_bestTourCostBuffer[i], 1, MPI_DOUBLE, i, BESTTOURCOST_T, MPI_COMM_WORLD, &recvBestTourCost[i]);
            MPI_Start(&recvBestTourCost[i]);

            // termination
            MPI_Send_init(&sendTerminationReq_buffer         , 1, MPI_INT, i, TERMINATION_T, MPI_COMM_WORLD, &sendTerminationReq[i]);
            MPI_Send_init(&sendTerminationReturn_buffer      , 1, MPI_INT, i, TERMINATION_T, MPI_COMM_WORLD, &sendTerminationReturn[i]);
            MPI_Recv_init(&recvTerminationReqReturn_buffer[i], 1, MPI_INT, i, TERMINATION_T, MPI_COMM_WORLD, &recvTerminationReqReturn[i]);
            MPI_Recv_init(&recvTerminationReqReturn_buffer[i], 1, MPI_INT, i, TERMINATION_T, MPI_COMM_WORLD, &recvTerminationReqReturn[i]);
            MPI_Start(&recvTerminationReqReturn[i]);
            
        }
    }

    void destroyMPICommunication(int number_of_processes) {
        for (int i = 0; i < number_of_processes; i++) { 
            MPI_Request_free(&sendNodeReturn[i]);
            MPI_Request_free(&sendNodeRequest[i]);
            MPI_Request_free(&recvNodeRequestReturn[i]);
            MPI_Request_free(&sendBestTourCost[i]);
            MPI_Request_free(&recvBestTourCost[i]);
            MPI_Request_free(&sendTerminationReq[i]);
            MPI_Request_free(&sendTerminationReturn[i]);
            MPI_Request_free(&recvTerminationReqReturn[i]);
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
        //std::cout << "Cities min1 and min2 info:" << std::endl;
        //for (int i = 0; i < _numberOfCities; i++) {
        //    std::cout << " -> City " << i << " | min1: " << _minPairs[i].first << "; min2: " << _minPairs[i].second << std::endl;
        //}
        return sum / 2;
    }

    double newLowerBound(int origin, int destiny, double LB, double cost) {
        double cf = cost >= _minPairs[origin].second ? _minPairs[origin].second : _minPairs[origin].first;
        double ct = cost >= _minPairs[destiny].second ? _minPairs[destiny].second : _minPairs[destiny].first;
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
        std::cout << _pid << " marking " << id << " as terminated" << std::endl;
        processes_state[id] = TERMINATED;
    }

    void mark_as_ready_to_terminate() {
        processes_state[_pid] = READY_TO_TERMINATE;
    }

    void revoke_terminated(int id) {
        std::cout << _pid << " revoking " << id <<  std::endl;
        processes_state[id] = NOT_TERMINATED;
    }

    void update_process_state(int id, int state) {
        std::cout << _pid << " marking " << id << " as " << state << std::endl;
        processes_state[id] = state;
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
            // std::cout << _pid << " running over the tour " << node_request_object->tour[node] << std::endl; 
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
                double newBestLowerBound = recv_bestTourCostBuffer[senderPid];
                if (newBestLowerBound < _bestLowerBound) {
                    _bestLowerBound = newBestLowerBound;
                }
                std::cout <<_pid << " Received LB " << newBestLowerBound << " From " << senderPid  << " " << std::endl;
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
        if (destPid == _pid) return; // don't send to itself
        node_request* node_request_object = getNodeReqFromByteArray(sendNodeBuffer_m, destPid);
        node_request_object->terminated = get_process_state();
        node_request_object->timestamp = ++sendTimestamps[destPid];
        node_request_object->type = NODERET_WITHOUTNODE_T;
        if (get_process_state() == NOT_TERMINATED && queue.size() > NODE_LIMIT) {
            node_request_object->type = NODERET_WITHNODE_T;
            std::shared_ptr<VisitedCity> city = queue.pop();
            node_request_object->cost = city->getCost();
            node_request_object->lowerBound = city->getLB();
            node_request_object->length = city->getLength();
            copyTourIntoNodeRequest(node_request_object, city->getTour(), city->getLength());
            revoke_terminated(destPid);
        }
        std::cout << _pid << " sending node return to " << destPid << std::endl;
        // std::cout << _pid << " sending flag as " << node_request_object->terminated << " to " << destPid  << std::endl;
        MPI_Wait(&sendNodeReturn[destPid], &status); // must wait before issuing next send
        MPI_Start(&sendNodeReturn[destPid]);
    }

    void mpi_testRecvNodeRequestReturn() {
        MPI_Status status;
        int finished;
        for (int i = 0; i < _numberOfProcesses; i++) {
            MPI_Test(&recvNodeRequestReturn[i], &finished, &status);
            if (finished) {
                int senderPid = status.MPI_SOURCE;
                if (senderPid != _pid) {
                    node_request* request = getNodeReqFromByteArray(recvNodeReqRetBuffer, senderPid);
                    if (request->type == NODEREQ_T) {
                        std::cout << _pid << " received a node request from " << status.MPI_SOURCE << " with state=" << request->terminated << std::endl;
                        update_state_if_new_message(senderPid, request->terminated, request->timestamp);
                        handleRecvNodeRequest(senderPid);
                    }
                    else if (request->type == NODERET_WITHOUTNODE_T || request->type == NODERET_WITHNODE_T) {
                        responseCount[senderPid]++; // keep track of all responses 
                        std::cout << _pid << " received node from " << status.MPI_SOURCE << " count is " << responseCount[senderPid] << std::endl;
                        update_state_if_new_message(senderPid, request->terminated, request->timestamp);
                        if (request->type == NODERET_WITHNODE_T && !request->terminated) {
                            addNodeToQueue(request);
                            revoke_terminated(_pid);
                        }
                    }
                }
                MPI_Start(&recvNodeRequestReturn[senderPid]);
            }
        }
    }

    void update_state_if_new_message(int senderPid, int state, int timestamp) {
        if (timestamp > recvTimestamps[senderPid]) {
            update_process_state(senderPid, state);
            recvTimestamps[senderPid] = timestamp;
        }
    }

    void broadcastNodeReq() {
        //std::cout << _pid << " broadcast start " << std::endl;
        std::vector<MPI_Status> statuses(_numberOfProcesses);
        sendNodeRequestBuffer->terminated = get_process_state();
        sendNodeRequestBuffer->type = NODEREQ_T;
        for (int i = 0; i < _numberOfProcesses; i++) {
            // attention lack of & before requestObj on start because of & in parameters
            if (i == _pid) continue;
            requestCount[i]++; // keep track of all sent messages
            sendNodeRequestBuffer->timestamp = ++sendTimestamps[i];
            MPI_Start(&sendNodeRequest[i]);
            std::cout << _pid << " requesting to " << i << " - count is " << requestCount[i] << std::endl;
        }
        int err_code = MPI_Waitall(_numberOfProcesses, sendNodeRequest.data(), statuses.data());
        if (err_code != MPI_SUCCESS) std::cout << _pid << " MPI_Waitall error " << std::endl;
        
        //std::cout << _pid << " broadcast end" << std::endl;
    }

        void broadcastBestTourCost() {
        //std::cout << _pid << " broadcast start " << std::endl;
        std::vector<MPI_Status> statuses(_numberOfProcesses);
        sendNodeRequestBuffer->terminated = get_process_state();
        for (int i = 0; i < _numberOfProcesses; i++) {
            // attention lack of & before requestObj on start because of & in parameters
            if (i == _pid) continue;
            sendNodeRequestBuffer->timestamp = ++sendTimestamps[i];
            MPI_Start(&sendBestTourCost[i]);
            std::cout << _pid << " sending to " << i << " with state= " << sendNodeRequestBuffer->terminated << std::endl;
        }
        int err_code = MPI_Waitall(_numberOfProcesses, sendBestTourCost.data(), statuses.data());
        if (err_code != MPI_SUCCESS) std::cout << _pid << " MPI_Waitall error " << std::endl;
        
        //std::cout << _pid << " broadcast end" << std::endl;
    }

    void mpi_shareBestToursAndDecideBest() {
        std::cout << _pid << " START gather results" << std::endl;
        char* recvBuffer = new char[NODE_REQUEST_SIZE * _numberOfProcesses];
        // create send buffer
        finalBestTourInfo->cost = _bestTourCost;
        std::cout << _pid << " cost = " << finalBestTourInfo->cost << std::endl;

        if (_bestTourCost < std::numeric_limits<double>::infinity()) { // might not have encoutered a full route
            std::cout << _pid << " copying tour" << std::endl;
            copyTourIntoNodeRequest(finalBestTourInfo, _bestTour, _numberOfCities + 1); // best tour starts and ends at 0, thus n_cities + 1
        }
        // gather
        std::cout << _pid << " gather results" << std::endl;
        MPI_Gather((char*) finalBestTourInfo, NODE_REQUEST_SIZE, MPI_BYTE, recvBuffer,
            NODE_REQUEST_SIZE , MPI_BYTE, 0, MPI_COMM_WORLD);
        std::cout << _pid << " after gather results" << std::endl;
        // choose the best from all received
        if (_pid == 0) {
            for (int i = 0; i < _numberOfProcesses; i++) {
                node_request* request = getNodeReqFromByteArray(recvBuffer, i);
                std::cout << _pid << " received cost of " << request->cost << " from " << i << std::endl;
                if ( request->cost < finalBestTourInfo->cost) {
                    std::cout << _pid << " choose from " << i << std::endl;
                    memcpy(finalBestTourInfo, request, NODE_REQUEST_SIZE);
                    std::cout << _pid << " has chosen a cost of " << finalBestTourInfo->cost << std::endl;
                }
            }
        }
        std::cout << _pid << " gather results end" << std::endl;
    }

    bool receivedAllSentNodes() {
        for (int i = 0 ; i < _numberOfProcesses ; i++) {
            std::cout << _pid << " respCount for " << i << " is " << responseCount[i] << std::endl;
            std::cout << _pid << " requestCount for " << i << " is " << requestCount[i] << std::endl;
            if (responseCount[i] != requestCount[i]) return false;
        }
        return true;
    }


    bool try_to_agree_on_terminate() {
        std::cout << _pid << " trying to terminate " << std::endl;
        mpi_testRecvNodeRequestReturn();
        mark_as_ready_to_terminate();
        MPI_Status status;
        bool canTerminate = true;
        // send to all
        for (int i = 0; i < _numberOfProcesses; i++) {
            if (i == _pid) continue;
            // std::cout << " A1 " << &sendTerminationReq[i] << std::endl;
            MPI_Wait(&sendTerminationReq[i], &status); // must wait before issuing next send
            MPI_Start(&sendTerminationReq[i]);
        }

        // wait to receive all

        for (int i = 0; i < _numberOfProcesses; i++) {
            if (i == _pid) continue;
            int finished;
            while(true) { // while waiting, still respond to other termination requests
                // mpi_testReceivedTerminationReq();
                // std::cout << " C1 " << &recvTerminationReturn[i] << std::endl;
                MPI_Test(&recvTerminationReqReturn[i], &finished, &status);
                if (finished) {
                    if (recvTerminationReqReturn_buffer[i] != READY_TO_TERMINATE) {
                        std::cout << _pid << " NOT READY TO TERMINATE " << std::endl; 
                        revoke_terminated(i);
                        canTerminate = false;
                    }
                    break;
                }
            }
            MPI_Start(&recvTerminationReqReturn[i]);
        }

        std::cout << _pid << " TERMINATED " << std::endl; 
        return canTerminate;
    }

    void mpi_testReceivedTerminationReq() {
        MPI_Status status;
        int finished;
        bool terminated;
        for (int i = 0; i < _numberOfProcesses; i++) {
            if (i == _pid) continue;
            MPI_Test(&recvTerminationReqReturn[i], &finished, &status);
            if (finished) {
                int senderPid = status.MPI_SOURCE;
                std::cout <<_pid << " Received terminationReq from " << senderPid  << " " << std::endl;
                sendTerminationReturn_buffer = get_process_state();
                MPI_Wait(&sendTerminationReturn[senderPid], &status); // must wait before issuing next send
                // std::cout << " D1 " << &sendTerminationReturn[senderPid] << std::endl;
                MPI_Start(&sendTerminationReturn[senderPid]);
                // std::cout << " E1 " << &recvTerminationReq[senderPid] << std::endl;
                MPI_Start(&recvTerminationReqReturn[senderPid]);
            }
        }
    }

    void findSolution() {
        std::shared_ptr<tour_t> tour;
        if (_pid == 0) {
            queue.push(std::make_shared<VisitedCity>(extendTour(nullptr, 0), 0, computeInitialLowerBound(), 1));
        }else computeInitialLowerBound();

        while(true) {
            // Branch and Bound main loop
            while (std::any_of(processes_state.begin(), processes_state.end(), [](short state) { return state == NOT_TERMINATED; })) {
                // for (int i = 0 ; i < _numberOfProcesses ; i++) {
                //     std::cout << _pid << " start while - state of " << i << " is " << processes_state[i] << std::endl;
                // }
                mpi_testReceivedTerminationReq();
                if (queue.empty()) {
                    // std::cout << _pid << " has empty queue" << std::endl;
                    mark_as_terminated(_pid);
                    // std::cout << _pid << " requesting nodes" << std::endl;
                    mpi_testRecvNodeRequestReturn();
                    broadcastNodeReq();
                    mpi_testRecvBestTourCost();

                    // Continue (while receiving messages and updating process_state)
                } else {

                    mpi_testRecvNodeRequestReturn();
                    mpi_testRecvBestTourCost();
                    std::cout << _pid << " queue size " << queue.size() << std::endl;
                    std::shared_ptr<VisitedCity> city = queue.pop();
                    tour = city->getTour();
                    double tourCost = city->getCost();
                    double bound = city->getLB();
                    int length = city->getLength();
                    int node = tour->node;
                    // std::cout << _pid << " has " << queue.size() << " nodes in queue" << std::endl;


                    //std::cout << ".. Visiting new city in the context of tour: (";
                    //showTour(city->getTour());
                    //std::cout << " ) | Node: " << node << "; Cost: " << tourCost <<
                    //    "; Lower-bound: " << bound << "; Length: " << length << "; Visited cities: ";
                    //printBits(tour->visitedCities);
                    //std::cout << std::endl;
                    std::cout << _pid << " looking at tour ending at " << node << " with length " << length << " and LB=" << bound << std::endl;
                    // std::cout << _pid << " bestLB=" << _bestLowerBound << " bestToutCost=" << _bestTourCost << std::endl;
                    if (bound >= _bestLowerBound) {
                        queue = PriorityQueue<std::shared_ptr<VisitedCity>, VisitedCity::CompareCityByLowerBound>();
                        // no need to set buffer (nullptr)
                        // std::cout << _pid << " bound >= bestTourCost" << std::endl;
                        mark_as_terminated(_pid);
                        mpi_testRecvNodeRequestReturn();
                        broadcastNodeReq();
                        continue;
                    }

                    // for (int i = 0 ; i < _numberOfProcesses ; i++) {
                    //     std::cout << _pid << " mid - state of " << i << " is " << processes_state[i] << std::endl;
                    // }
                    
                    if (length == _numberOfCities) {
                        // for (int i = 0 ; i < _numberOfProcesses ; i++) {
                        //     std::cout << _pid << " length==n_cities - state of " << i << " is " << processes_state[i] << std::endl;
                        // }
                        double costUntilEnd = tourCost + _roadsCost[node][0];
                        // std::cout << _pid << " length==n_cities. costUntilEnd= " << costUntilEnd << " with bestTourCost " << _bestTourCost << std::endl;
                        if (costUntilEnd < _bestTourCost) {
                            if (!_bestTour) {
                                _bestTour = std::make_shared<tour_t>();
                                _bestTour->node = 0;
                            }
                            _bestTour->prev = tour;
                            _bestTourCost = costUntilEnd;
                            if (_bestTourCost < _bestLowerBound) {
                                _bestLowerBound = _bestTourCost;
                                bestTourCostBuffer = _bestTourCost;
                                std::cout << _pid << " sending bestTourCost = " << bestTourCostBuffer << std::endl;
                                broadcastBestTourCost();
                            }
                            //showTour(_bestTour);
                            //std::cout << " -- Tour cost: " << _bestTourCost << std::endl;
                        }
                    }
                    else {
                        // for (int i = 0 ; i < _numberOfProcesses ; i++) {
                        //     std::cout << _pid << " begin else - state of " << i << " is " << processes_state[i] << std::endl;
                        // }      
                        for (const int& destiny : _citiesNeighbors[node]) {
                            if (hasNotVisitedCity(tour->visitedCities, destiny)) {

                                double cost = _roadsCost[node][destiny];
                                double newBound = newLowerBound(node, destiny, bound, cost);
                                if (newBound > _bestLowerBound) {
                                    // std::cout << _pid << " skip node " << destiny << " with newBound=" << newBound << std::endl;
                                    continue;
                                }
                                queue.push(std::make_shared<VisitedCity>(extendTour(tour, destiny), tourCost + cost, newBound, length + 1));
                            }
                        }
                        // for (int i = 0 ; i < _numberOfProcesses ; i++) {
                        //     std::cout << _pid << " end else - state of " << i << " is " << processes_state[i] << std::endl;
                        // }    
                    }
                    // std::cout << _pid << " has " << queue.size() << " nodes" << std::endl;
                    // std::cout << _pid << " condition is " << std::any_of(processes_state.begin(), processes_state.end(), [](short state) { return state == NOT_TERMINATED; }) << std::endl;
                    for (int i = 0 ; i < _numberOfProcesses ; i++) {
                        std::cout << _pid << " end while - state of " << i << " is " << processes_state[i] << std::endl;
                    }
                }
            }

            std::cout << _pid << " finished findsol" << std::endl;
            // wait for all delayed expected nodes
            while (!receivedAllSentNodes()) {
                mpi_testRecvNodeRequestReturn();
                mpi_testReceivedTerminationReq();
            }
            std::cout << _pid << " waited for all responses " << std::endl;

            // TODO initiate termination messages
            if (try_to_agree_on_terminate()) break;
        }

        std::cout << _pid << " stop " << std::endl;
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
    std::cout << pid << " exited findSol" << std::endl;
    elapsed_time += MPI_Wtime();

    tsp.destroyMPICommunication(number_of_processes);
    MPI_Finalize();

    if (pid == 0) {
        std::cout << "----- RESULTS ------" << std::endl;
        fprintf(stderr, "Time of execution: %.6fs\n", elapsed_time);
        tsp.print_result();
    }

    // Attention! please 

    return 0;
}
