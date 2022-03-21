/*
 * Copyright (c) 2008, Maxim Likhachev
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Carnegie Mellon University nor the names of its
 *       contributors may be used to endorse or promote products derived from
 *       this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * 
 * THE EDITS ARE MADE BY MOHAMED KHALIL KAHLAOUI ON MARCH 2022
 * THIS CODE SUPPORTS THE "Rapidly-exploring random tree" ALGORITHM WITH A DOUBLE TREE SEARCH AND GENERATES A SMOOTH PATH USING "Chaikin" ALGORITHM  
 *  
 */

#include <sbpl/planners/rrtplanner.h>

#include <random>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <tuple>
#include <math.h>
#include <limits>
#include <list>

#include <sbpl/discrete_space_information/environment.h>
#include <sbpl/discrete_space_information/environment_nav2D.h>
#include <sbpl/utils/heap.h>
#include <sbpl/utils/key.h>
#include <sbpl/utils/list.h>
#include <sbpl/utils/tree.hh>
#include <sbpl/utils/tree_util.hh>

using namespace std;
//GLBALE VARIABLES FOR THE DECLARATION
double GOAL;
double xX,yY;
bool noOBS=true;
int env_sizeX1;
int env_sizeY1;
std::vector<double> tempX;
std::vector<double> tempY;

std::vector<double> nodex;
std::vector<double> nodey;

std::vector<double> obsXmin{};
std::vector<double> obsYmin{};

std::vector<double> pathX{};
std::vector<double> pathY{};

std::vector<double> FinalpathX{};
std::vector<double> FinalpathY{};

std::vector<double> nx;
std::vector<double> ny;
std::vector<double> x11;
std::vector<double> y11;

void RRTPlanner::get_envinformation()
{
	int env_sizeX, env_sizeY, env_startX, env_startY, env_goalX, env_goalY; unsigned char env_obsthresh;
	((EnvironmentNAV2D *)environment_)->GetEnvParms(&env_sizeX, &env_sizeY, &env_startX, &env_startY, &env_goalX, &env_goalY, &env_obsthresh);
	std::cout << "Environment Header Information" << std::endl << std::endl;
	std::cout << "Size:  " << env_sizeX  << ", " << env_sizeY  << std::endl;
	std::cout << "Start: " << env_startX << ", " << env_startY << std::endl;
	std::cout << "Goal:  " << env_goalX  << ", " << env_goalY  << std::endl;
	std::cout << "Obstacle Thresh: " << (int)env_obsthresh << std::endl << std::endl;
	
	std::cout << "Environment Map Information" << std::endl << std::endl;
	int x, y;
	unsigned char cost;
	for (y=0; y<env_sizeY; y++) {
		for (x=0; x<env_sizeX; x++) {
			cost = ((EnvironmentNAV2D *)environment_)->GetMapCost(x,y);
			std::cout << (int)cost << " "; 
		}
		std::cout << std::endl;
	}

	std::cout << std::endl << "Start (S) + Goal (G) + Obstacle Information (*)" << std::endl << std::endl;
	bool isobs;
	for (y=0; y<env_sizeY; y++) {
		for (x=0; x<env_sizeX; x++) {
			if (x == env_startX && y == env_startY) {
				// Starting Point
				std::cout << "S ";
			}
			else if (x == env_goalX && y == env_goalY) {
				// Goal Point
				std::cout << "G ";
			}
			else {
				isobs = ((EnvironmentNAV2D *)environment_)->IsObstacle(x,y);
				if (isobs) {
					// Obstacle
					std::cout << "* ";
					obsXmin.push_back(x);
					obsYmin.push_back(y);				
				}
				else {
					// Free Space
					std::cout << "0 ";
				}
			}
		}
        std::cout << std::endl;
	}
    std::cout << std::endl;
}

//CHAIKIN PROCEDURE IS USED TO CREATE A SMOOTHER PATH FOR THE ROBOT 
void chaikin(double attr_chaikiniter, std::tuple<std::vector<double>, std::vector<double>> Finalpath){

    int lastX,lastY;
    float percent1;
    x11.push_back(FinalpathX.at(0));
    y11.push_back(FinalpathY.at(0));
    //PREPARING THE SERIE OF POINTS FOR THE CHAIKIN ALGORITHM 
    for(int i=0; i < FinalpathX.size(); i++){

        nx.push_back(FinalpathX.at(i));
        ny.push_back(FinalpathY.at(i));
        x11.push_back(FinalpathX.at(i));
        y11.push_back(FinalpathY.at(i));

        int firstX = FinalpathX.at(0);
        int firstY = FinalpathY.at(0);

         lastX = FinalpathX.at(FinalpathX.size() - 1);
         lastY = FinalpathY.at(FinalpathY.size() - 1);
         percent1 =0.5;
    }
    nx.clear(); ny.clear();
    //STARTING THE CHAIKIN ITERATIONS
    for (int i = 0; i < attr_chaikiniter; i++) {

        for (int i = 0; i < x11.size() - 1; i++) {

       int dx=x11.at(i+1)-x11.at(i);
       int dy=y11.at(i+1)-y11.at(i);

       int xvalue=x11.at(i)+dx*percent1;
       int xvalue1=x11.at(i)+dx*(1-percent1);

       int yvalue=(y11.at(i)+dy*percent1);
       int yvalue1=y11.at(i)+dy*(1-percent1);

       //CHECKING IF THERE AN OBSTACLE IN THE NEW POINTS OR NOT
         for (int j = 0; j < obsXmin.size(); j++)
        {
            if (obsXmin.at(j)==xvalue && obsYmin.at(j)==yvalue || obsXmin.at(j)==xvalue1 && obsYmin.at(j)==yvalue1)
            {
             noOBS=false;
            }
            else
            nx.push_back(x11.at(i)+dx*percent1);
            nx.push_back(x11.at(i)+dx*(1-percent1));

            ny.push_back(y11.at(i)+dy*percent1);
            ny.push_back(y11.at(i)+dy*(1-percent1)); 
        }  
     }
     nx.push_back(lastX);
     ny.push_back(lastY);
     
     }
}

//BRESENHAM PROCEDURE IS USED TO GENERATE A GROUP OF POINTS BETWEEN TWO SEPERATE POINTS IN A CARTISIAN SYSTEM (X,Y) 
void bresenham(int x1, int y1, int x2, int y2, int dx, int dy, int decide)
{
    //pk is initial decesion making parameter
    //Note:x1&y1,x2&y2, dx&dy values are intercganged
    //and passed in plotPixel function so
    //it can handle both cases when m>1 & m<1
    int pk = 2 * dy - dx;
    for (int i = 0; i <= dx; i++)
    {
        if(decide==0){
            pathX.push_back(x1);
            pathY.push_back(y1);
    }
        if(decide==1){
            pathX.push_back(y1);
            pathY.push_back(x1);
    }
        //checking either to decrement or increment the value
        //if we have to plot from (0,100) to (100,0)
        x1 < x2 ? x1++ : x1--;
        if (pk < 0)
        {
            //decesion value will decide to plot
            //either  x1 or y1 in x's position
            if (decide == 0)
            {
               // putpixel(x1, y1, RED);
                pk = pk + 2 * dy;
            }
            else
            {
                //(y1,x1) is passed in xt
               // putpixel(y1, x1, YELLOW);
                pk = pk + 2 * dy;
            }
        }
        else
        {
            y1 < y2 ? y1++ : y1--;
            if (decide == 0)
            {
 
                //putpixel(x1, y1, RED);
            }
            else
            {
              //  putpixel(y1, x1, YELLOW);
            }
            pk = pk + 2 * dy - 2 * dx;
        }
    }
}

struct Node
{

    double posX;
    double posY;
    Node *prev;
    Node *next;
};


class RRT
{

private:
    Node *start;
    Node *goal;
    std::vector<Node *> rrtNodes;

public:
    RRT(double startX, double startY, double goalX, double goalY, int env_sizeX, int env_sizeY)
    {

        Node *node = new Node;
        this->start = node;
        start->posX = startX;
        start->posY = startY;
        start->prev = nullptr;
        rrtNodes.push_back(node);

        node = new Node;
        this->goal = node;
        goal->posX = goalX;
        goal->posY = goalY;

        env_sizeX1=env_sizeX;
        env_sizeY1=env_sizeY;
    }
//THIS FUNCTION CHECK IF THERE IS AN OBSTACLE BETWEEN TWO NODES OR NOT
    bool checkObstacles(double x1, double y1, double x2, double y2, std::tuple<std::vector<double>, std::vector<double>> obs)
    {
         bool decision;
         int dx, dy, pk;
         int counter;
    //cin cout
    dx = abs(x2 - x1);
    dy = abs(y2 - y1);
    //If slope is less than one
    if (dx > dy)
    {
        //passing argument as 0 to plot(x,y)
        bresenham(x1, y1, x2, y2, dx, dy, 0);
    }
    //if slope is greater than or equal to 1
    else
    {
        //passing argument as 1 to plot (y,x)
        bresenham(y1, x1, y2, x2, dy, dx, 1);
    }
        counter=0;

        for (int j = 0; j < pathX.size(); j++)
        {

             for (int i = 0; i < obsXmin.size(); i++)
                {
            if (obsXmin.at(i)==pathX.at(j) && obsYmin.at(i)==pathY.at(j))
            {
                counter ++ ;
               
            }
            else
                decision = false;
                }
        }

        if(x1>=env_sizeX1 ||y1>=env_sizeY1 ||x2>=env_sizeX1 ||y2>=env_sizeY1){
            counter++;
        }
        if(counter!=0){
            decision= true;
        }
    pathX.clear();
    pathY.clear();
    return decision;
    }
//THIS FUNCTION IS TO CHECK THE NEAREST NODE TO THE NEW RANDOM ONE 
      Node *checkNearestNode(Node *new_node,double attr_stepsize)
    {

        Node *near_node = new Node;
        double minDistance = std::numeric_limits<double>::max();

        int corrX = 0.0;
        int corrY = 0.0;
        bool check_obstacle;

        for (auto &ii : rrtNodes)
        {
            double distance = std::sqrt(std::pow((new_node->posX - ii->posX), 2) + std::pow((new_node->posY - ii->posY), 2));
            if (distance < minDistance)
            {
                minDistance = distance;
                near_node = ii;   
            }
        }
        double dx = new_node->posX - near_node->posX;
        double dy = new_node->posY - near_node->posY;
        double angle = std::atan2(dy, dx) * 180 / M_PI;

        if (minDistance > attr_stepsize)
        {
            //CALCULATING THE NEW POSITION OF THE RANDOMLY GENERATE NEW NODE BECAUSE OF THE DISTANCE
            corrX = std::abs(near_node->posX + std::cos(angle) * attr_stepsize);
            corrY = std::abs(near_node->posY + std::sin(angle) * attr_stepsize);
        }

        if (minDistance <= attr_stepsize)
        {
            //MOVE THE RANDOMLY GENERATE NEW NODE BECAUSE OF THE DISTANCE
            corrX = new_node->posX;
            corrY = new_node->posY;
        }
        if (rrtNodes.size() > 0)
        {
            //CHECK IF THERE IS AN OBSTACLE BETWEEN TWO NODES
            check_obstacle = checkObstacles(near_node->posX, near_node->posY, corrX, corrY, std::make_tuple(obsXmin, obsYmin));
        }

        new_node->posX = corrX;
        new_node->posY = corrY;

        near_node->next = new_node;
        new_node->prev = near_node;

        if (rrtNodes.size() == 0)
        {
            //SAVING THE FIRST TWO NODE POSITIONS FOR THE 2 TREE SEARCH 
            new_node->prev = start;
            new_node->prev = goal;
        }

        if (check_obstacle == 0)
        {
            //SAVE THE NEW NODE IF THE THERE IS NO OBSTACLES
            rrtNodes.push_back(new_node);
        }

        if (((double)new_node->posX == (double)this->goal->posX) && ((double)new_node->posY == (double)this->goal->posY) && check_obstacle == 0)
        {
            std::cout << "The GOAL achive NODE POSITIONS:" << std::endl;
            GOAL = 1;
            while (new_node->prev != nullptr)
            {
                //WHILE THE GOAL IS REACHED THE FINAL NODES WOULD BE SAVED 
                std::cout << new_node->posX << " :?: " << new_node->posY << std::endl;
                new_node = new_node->prev;
                tempX.push_back(new_node->posX);
                tempY.push_back(new_node->posY);
            }

        }
        return new_node;
    }

    double lookForPath(double attr_maxplanningtime,double attr_stepsize)
    {
        std::random_device dev;
        std::mt19937 rng(dev());
        std::mt19937 rng1(dev());
        std::uniform_int_distribution<std::mt19937::result_type> dist100(0, env_sizeX1); // distribution in range [1, 6]
        std::uniform_int_distribution<std::mt19937::result_type> dist10(0, env_sizeY1);

           clock_t TimeStarted = clock(); 

    while ((clock() - TimeStarted) < (attr_maxplanningtime * (double)CLOCKS_PER_SEC)) {

            Node *random_node = new Node;
            Node *last_node = new Node;
            //GENERATE THE RANDOM NODES
            int randX = dist100(rng);
            int randY = dist100(rng1);
            random_node->posX = randX;
            random_node->posY = randY;

            //checkNearestNode a function that take a decision to obtain the recommended node or reject it 
            last_node = checkNearestNode(random_node,attr_stepsize);

            if (GOAL == 1)
            {
                goal->prev = last_node;
               std::cout << last_node->posX << " :: " << last_node->posY << std::endl;

                return -1;
            }
        }
        return 1;
    }
};



std::string RRTPlanner::TreeTypeToStr(errt_TreeType treetype)
{
    switch (treetype) {
		case SINGLE_TREE_TYPE:
			return std::string("Single Tree");
		case DOUBLE_TREE_TYPE:
			return std::string("Double Tree");
		default:
			return std::string("Invalid");
    }
}
 
std::string RRTPlanner::PathTypeToStr(errt_PathType pathtype)
{
    switch (pathtype) {
		case NOTSMOOTHED_PATH_TYPE:
			return std::string("Not Smoothed");
		case SMOOTHED_PATH_TYPE:
			return std::string("Smoothed");
		default:
			return std::string("Invalid");
    }
}



RRTPlanner::RRTPlanner(DiscreteSpaceInformation* environment, errt_TreeType treetype, errt_PathType pathtype,
                       double stepsize, double qgoalprob, int chaikiniter, double maxplanningtime)
{
	prevConstructor(environment, false);

	// Saving Parameters
	std::cout << std::endl << "**********" << std::endl << "Parameters" << std::endl << "**********" << std::endl;
	attr_treetype        = treetype;
	attr_pathtype        = pathtype;
	attr_stepsize        = stepsize;
	attr_qgoalprob       = qgoalprob;
	attr_chaikiniter     = chaikiniter;
	attr_maxplanningtime = maxplanningtime;
	std::cout << "Tree Type -> " << TreeTypeToStr(attr_treetype) << std::endl;
	std::cout << "Path Type -> " << PathTypeToStr(attr_pathtype) << std::endl;
	std::cout << "Step Size -> " << attr_stepsize    << std::endl;
	std::cout << "qgoal Probability   -> " << attr_qgoalprob   << std::endl;
	std::cout << "Chaikin Iterations  -> " << attr_chaikiniter << std::endl;
	std::cout << "Max. Planning Time  -> " << attr_maxplanningtime << " seconds" << std::endl;


	// Accessing the Environment Information
	std::cout << std::endl << "**************************************" << std::endl << "Accessing the Environment Information!" << std::endl << "**************************************" << std::endl;
	//get_envinformation();
}

void RRTPlanner::prevConstructor(DiscreteSpaceInformation* environment, bool bSearchForward)
{
	environment_   = environment;
    bforwardsearch = bSearchForward;

    bsearchuntilfirstsolution = false;
    finitial_eps              = ARA_DEFAULT_INITIAL_EPS;
    final_epsilon             = ARA_FINAL_EPS;
    dec_eps                   = ARA_DECREASE_EPS;
    use_repair_time           = false;
    repair_time               = INFINITECOST;
    searchexpands             = 0;
    MaxMemoryCounter          = 0;

#ifndef ROS
    const char* debug = "debug.txt";
#endif

    fDeb = SBPL_FOPEN(debug, "w");
    if (fDeb == NULL) {
        throw SBPL_Exception("ERROR: could not open planner debug file");
    }

    pSearchStateSpace_ = new RRTSearchStateSpace_t;

    //create the RRT planner
    if (CreateSearchStateSpace(pSearchStateSpace_) != 1) {
        SBPL_ERROR("ERROR: failed to create statespace\n");
        return;
    }

    //set the start and goal states
    if (InitializeSearchStateSpace(pSearchStateSpace_) != 1) {
        SBPL_ERROR("ERROR: failed to create statespace\n");
        return;
    }
    finitial_eps_planning_time      = -1.0;
    final_eps_planning_time         = -1.0;
    num_of_expands_initial_solution = 0;
    final_eps                       = -1.0;
}

void RRTPlanner::ex_randnumbergenerator(int setofnumbers)
{
	std::cout << std::endl << std::endl << "********************************" << std::endl << "Example: Random Number Generator" << std::endl << "********************************" << std::endl;
	srand(time(NULL));
	for (int i=0; i<setofnumbers; i++) {
		// rand() returns a random number in the range between 0 and RAND_MAX
		std::cout << std::fixed << std::setprecision(5) << rand()/(float)RAND_MAX << std::endl;
	}
}

void RRTPlanner::savepath(vector<RRTMapCoord *>* path_mapcoord_V, vector<int>* solution_stateIDs_V)
{
	// coordinates are converted into states
	RRTMapCoord *mapcoord; int stateID;
	solution_stateIDs_V->clear();
	for (int i=0; i<path_mapcoord_V->size(); i++) {
		mapcoord = path_mapcoord_V->at(i);
		stateID  = ((EnvironmentNAV2D *)environment_)->GetStateFromCoord(mapcoord->attr_x, mapcoord->attr_y);
		solution_stateIDs_V->push_back(stateID);
	}
}

int RRTPlanner::compute_pathcost(vector<RRTMapCoord *>* path_mapcoord_V)
{
	int pathcost=0, cellcost, cellcost1, cellcost2, cellcost3;
	RRTMapCoord *mapcoord1, *mapcoord2;
	
	mapcoord1 = path_mapcoord_V->at(0);
	for (int i=1; i<path_mapcoord_V->size(); i++) {
		mapcoord2 = path_mapcoord_V->at(i);

		if ((abs(mapcoord1->attr_x - mapcoord2->attr_x) + abs(mapcoord1->attr_y - mapcoord2->attr_y)) > 1)
		{
			// diagonal movement
			cellcost1 = ((EnvironmentNAV2D *)environment_)->GetMapCost(mapcoord2->attr_x, mapcoord2->attr_y);
			cellcost2 = ((EnvironmentNAV2D *)environment_)->GetMapCost(mapcoord2->attr_x, mapcoord1->attr_y);
			cellcost3 = ((EnvironmentNAV2D *)environment_)->GetMapCost(mapcoord1->attr_x, mapcoord2->attr_y);

			cellcost  = cellcost1;
			if (cellcost2 > cellcost) cellcost = cellcost2;
			if (cellcost3 > cellcost) cellcost = cellcost3;

			pathcost += ((cellcost + 1) * ENVNAV2D_COSTMULT * sqrt(2));
		}
		else {
			// horizontal/vertical movement
			cellcost  = ((EnvironmentNAV2D *)environment_)->GetMapCost(mapcoord2->attr_x, mapcoord2->attr_y);
			pathcost += ((cellcost + 1) * ENVNAV2D_COSTMULT);
		}

		mapcoord1 = mapcoord2;
	}
	return(pathcost);
}


int RRTPlanner::replan(vector<int>* solution_stateIDs_V, ReplanParams params, int* solcost)
{
    vector<RRTMapCoord *> path_mapcoord_V;
	int env_sizeX, env_sizeY, env_startX, env_startY, env_goalX, env_goalY; unsigned char env_obsthresh;
	((EnvironmentNAV2D *)environment_)->GetEnvParms(&env_sizeX, &env_sizeY, &env_startX, &env_startY, &env_goalX, &env_goalY, &env_obsthresh);
		get_envinformation();

        //Launching of the RRT algortihm with inicializing the starting and the goal points
        RRT rrt(env_startX, env_startY,env_goalX, env_goalY,env_sizeX, env_sizeY);
        //Start looking for the path 
        double checkpath = rrt.lookForPath(attr_maxplanningtime,attr_stepsize);


         int dx, dy;
        nodex.push_back(env_goalX);
        nodey.push_back(env_goalY);
        
         for(int i=0; i < tempX.size(); i++){
          
            nodex.push_back(tempX.at(i));
            nodey.push_back(tempY.at(i));
        }
            tempX=nodex;
            tempY=nodey;

for(int i=0; i < tempX.size()-1; i++){
// CREATING THE PATH THROUGH THE NODES
double x1=tempX.at(i);
double x2=tempX.at(i+1);

double y1=tempY.at(i);
double y2=tempY.at(i+1);

    dx = abs(x2 - x1);
    dy = abs(y2 - y1);
    //If slope is less than one
    if (dx > dy)
    {
        //passing argument as 0 to plot(x,y)
        bresenham(x1, y1, x2, y2, dx, dy, 0);
    }
    //if slope is greater than or equal to 1
    else
    {
        //passing argument as 1 to plot (y,x)
        bresenham(y1, x1, y2, x2, dy, dx, 1);
    }
for(int j=0; j < pathX.size()-1; j++){
    FinalpathX.push_back(pathX.at(j));
    FinalpathY.push_back(pathY.at(j));
}
    pathX.clear();
    pathY.clear();
}

if(checkpath==-1 && PathTypeToStr(attr_pathtype)=="Smoothed"){
//SOOMTHING THE PATH WE HAVE IF THE SMOOTHING WILL TOUCH AN OBSTACLE IT IS GOING TO BE SKIPPED
    chaikin(attr_chaikiniter, std::make_tuple(FinalpathX, FinalpathY));
  
    if(noOBS==true)
      for(int i=0; i < nx.size(); i++){
    path_mapcoord_V.push_back(new RRTMapCoord(nx.at(i),ny.at(i)));
    }
    }
string value=PathTypeToStr(attr_pathtype);
    if(noOBS==false){

    value="Not Smoothed" ;
    }
if(checkpath==-1 && value=="Not Smoothed" ){

    std::cout << "NO OB11 " << std::endl;
        for(int i=0; i < FinalpathX.size(); i++){
        path_mapcoord_V.push_back(new RRTMapCoord(FinalpathX.at(i),FinalpathY.at(i)));   
        }
    }    
	// if a path exists:
	// 		- save the path into the solution_stateIDs_V vector
	// 		- compute the path cost
	// 		- return 1

		// Example: it is assumed that the path goes through the following cells:
    if(checkpath==-1){
		savepath(&path_mapcoord_V, solution_stateIDs_V);
		*solcost = compute_pathcost(&path_mapcoord_V);  
	
		// free memory
        
		for (int i=0; i<path_mapcoord_V.size(); i++)
		{ delete(path_mapcoord_V.at(i)); }

		return 1;
	}
	// if a path does not exist:
	// 		- return 0
	else
		return 0;
    
}

// ************************************************
// IT IS NOT NECESSARY TO SEE THE REST OF THE CODE!
// ************************************************

RRTPlanner::~RRTPlanner()
{
    if (pSearchStateSpace_ != NULL) {
        //delete the statespace
        DeleteSearchStateSpace( pSearchStateSpace_);
        delete pSearchStateSpace_;
    }
    SBPL_FCLOSE( fDeb);
}

void RRTPlanner::Initialize_searchinfo(CMDPSTATE* state, RRTSearchStateSpace_t* pSearchStateSpace)
{
    RRTState* searchstateinfo = (RRTState*)state->PlannerSpecificData;

    searchstateinfo->MDPstate = state;
    InitializeSearchStateInfo(searchstateinfo, pSearchStateSpace);
}

CMDPSTATE* RRTPlanner::CreateState(int stateID, RRTSearchStateSpace_t* pSearchStateSpace)
{
    CMDPSTATE* state = NULL;

#if DEBUG
    if (environment_->StateID2IndexMapping[stateID][ARAMDP_STATEID2IND] != -1) {
        throw SBPL_Exception("ERROR in CreateState: state already created");
    }
#endif

    //adds to the tail a state
    state = pSearchStateSpace->searchMDP.AddState(stateID);

    //remember the index of the state
    environment_->StateID2IndexMapping[stateID][ARAMDP_STATEID2IND] =
            pSearchStateSpace->searchMDP.StateArray.size() - 1;

#if DEBUG
    if(state !=
       pSearchStateSpace->searchMDP.StateArray[environment_->StateID2IndexMapping[stateID][ARAMDP_STATEID2IND]])
    {
        throw SBPL_Exception("ERROR in CreateState: invalid state index");
    }
#endif

    //create search specific info
    state->PlannerSpecificData = (RRTState*)malloc(sizeof(RRTState));
    Initialize_searchinfo(state, pSearchStateSpace);
    MaxMemoryCounter += sizeof(RRTState);

    return state;
}

CMDPSTATE* RRTPlanner::GetState(int stateID, RRTSearchStateSpace_t* pSearchStateSpace)
{
    if (stateID >= (int)environment_->StateID2IndexMapping.size()) {
        std::stringstream ss("ERROR int GetState: stateID ");
        ss << stateID << " is invalid";
        throw SBPL_Exception(ss.str());
    }

    if (environment_->StateID2IndexMapping[stateID][ARAMDP_STATEID2IND] == -1)
        return CreateState(stateID, pSearchStateSpace);
    else
        return pSearchStateSpace->searchMDP.StateArray[environment_->StateID2IndexMapping[stateID][ARAMDP_STATEID2IND]];
}

//-----------------------------------------------------------------------------------------------------

int RRTPlanner::ComputeHeuristic(CMDPSTATE* MDPstate, RRTSearchStateSpace_t* pSearchStateSpace)
{
    //compute heuristic for search

    if (bforwardsearch) {

#if MEM_CHECK == 1
        //int WasEn = DisableMemCheck();
#endif

        //forward search: heur = distance from state to searchgoal which is Goal RRTState
        int retv = environment_->GetGoalHeuristic(MDPstate->StateID);

#if MEM_CHECK == 1
        //if (WasEn)
        //	EnableMemCheck();
#endif

        return retv;

    }
    else {
        //backward search: heur = distance from searchgoal to state
        return environment_->GetStartHeuristic(MDPstate->StateID);
    }
}

// initialization of a state
void RRTPlanner::InitializeSearchStateInfo(RRTState* state, RRTSearchStateSpace_t* pSearchStateSpace)
{
    state->g = INFINITECOST;
    state->v = INFINITECOST;
    state->iterationclosed = 0;
    state->callnumberaccessed = pSearchStateSpace->callnumber;
    state->bestnextstate = NULL;
    state->costtobestnextstate = INFINITECOST;
    state->heapindex = 0;
    state->listelem[ARA_INCONS_LIST_ID] = 0;
#if DEBUG
    state->numofexpands = 0;
#endif

    state->bestpredstate = NULL;

    //compute heuristics
#if USE_HEUR
    if(pSearchStateSpace->searchgoalstate != NULL)
        state->h = ComputeHeuristic(state->MDPstate, pSearchStateSpace);
    else
        state->h = 0;
#else
        state->h = 0;
#endif
}

//re-initialization of a state
void RRTPlanner::ReInitializeSearchStateInfo(RRTState* state, RRTSearchStateSpace_t* pSearchStateSpace)
{
    state->g = INFINITECOST;
    state->v = INFINITECOST;
    state->iterationclosed = 0;
    state->callnumberaccessed = pSearchStateSpace->callnumber;
    state->bestnextstate = NULL;
    state->costtobestnextstate = INFINITECOST;
    state->heapindex = 0;
    state->listelem[ARA_INCONS_LIST_ID] = 0;
#if DEBUG
    state->numofexpands = 0;
#endif

    state->bestpredstate = NULL;

    //compute heuristics
#if USE_HEUR

    if(pSearchStateSpace->searchgoalstate != NULL)
    {
        state->h = ComputeHeuristic(state->MDPstate, pSearchStateSpace);
    }
    else
    state->h = 0;

#else

    state->h = 0;

#endif
}

void RRTPlanner::DeleteSearchStateData(RRTState* state)
{
    //no memory was allocated
    MaxMemoryCounter = 0;
    return;
}

//used for backward search
void RRTPlanner::UpdatePreds(RRTState* state, RRTSearchStateSpace_t* pSearchStateSpace)
{
    vector<int> PredIDV;
    vector<int> CostV;
    CKey key;
    RRTState *predstate;

    environment_->GetPreds(state->MDPstate->StateID, &PredIDV, &CostV);

    //iterate through predecessors of s
    for (int pind = 0; pind < (int)PredIDV.size(); pind++) {
        CMDPSTATE* PredMDPState = GetState(PredIDV[pind], pSearchStateSpace);
        predstate = (RRTState*)(PredMDPState->PlannerSpecificData);
        if (predstate->callnumberaccessed != pSearchStateSpace->callnumber) {
            ReInitializeSearchStateInfo(predstate, pSearchStateSpace);
        }

        //see if we can improve the value of predstate
        if (predstate->g > state->v + CostV[pind]) {
            predstate->g = state->v + CostV[pind];
            predstate->bestnextstate = state->MDPstate;
            predstate->costtobestnextstate = CostV[pind];

            //re-insert into heap if not closed yet
            if (predstate->iterationclosed != pSearchStateSpace->searchiteration) {
                key.key[0] = predstate->g + (int)(pSearchStateSpace->eps * predstate->h);
                //key.key[1] = predstate->h;
                if (predstate->heapindex != 0)
                    pSearchStateSpace->heap->updateheap(predstate, key);
                else
                    pSearchStateSpace->heap->insertheap(predstate, key);
            }
            //take care of incons list
            else if (predstate->listelem[ARA_INCONS_LIST_ID] == NULL) {
                pSearchStateSpace->inconslist->insert(predstate, ARA_INCONS_LIST_ID);
            }
        }
    } //for predecessors
}

//used for forward search
void RRTPlanner::UpdateSuccs(RRTState* state, RRTSearchStateSpace_t* pSearchStateSpace)
{
    vector<int> SuccIDV;
    vector<int> CostV;
    CKey key;
    RRTState *succstate;

    environment_->GetSuccs(state->MDPstate->StateID, &SuccIDV, &CostV);

    //iterate through predecessors of s
    for (int sind = 0; sind < (int)SuccIDV.size(); sind++) {
        CMDPSTATE* SuccMDPState = GetState(SuccIDV[sind], pSearchStateSpace);
        int cost = CostV[sind];

        succstate = (RRTState*)(SuccMDPState->PlannerSpecificData);
        if (succstate->callnumberaccessed != pSearchStateSpace->callnumber) {
            ReInitializeSearchStateInfo(succstate, pSearchStateSpace);
        }

        //see if we can improve the value of succstate
        //taking into account the cost of action
        if (succstate->g > state->v + cost) {
            succstate->g = state->v + cost;
            succstate->bestpredstate = state->MDPstate;

            //re-insert into heap if not closed yet
            if (succstate->iterationclosed != pSearchStateSpace->searchiteration) {

                key.key[0] = succstate->g + (int)(pSearchStateSpace->eps * succstate->h);

                //key.key[1] = succstate->h;

                if (succstate->heapindex != 0)
                    pSearchStateSpace->heap->updateheap(succstate, key);
                else
                    pSearchStateSpace->heap->insertheap(succstate, key);
            }
            //take care of incons list
            else if (succstate->listelem[ARA_INCONS_LIST_ID] == NULL) {
                pSearchStateSpace->inconslist->insert(succstate, ARA_INCONS_LIST_ID);
            }
        } //check for cost improvement
    } //for actions
}

//TODO-debugmax - add obsthresh and other thresholds to other environments in 3dkin
int RRTPlanner::GetGVal(int StateID, RRTSearchStateSpace_t* pSearchStateSpace)
{
    CMDPSTATE* cmdp_state = GetState(StateID, pSearchStateSpace);
    RRTState* state = (RRTState*)cmdp_state->PlannerSpecificData;
    return state->g;
}

//returns 1 if the solution is found, 0 if the solution does not exist and 2 if it ran out of time
int RRTPlanner::ImprovePath(RRTSearchStateSpace_t* pSearchStateSpace, double MaxNumofSecs)
{
    int expands;
    RRTState *state, *searchgoalstate;
    CKey key, minkey;
    CKey goalkey;

    expands = 0;

    if (pSearchStateSpace->searchgoalstate == NULL) {
        throw SBPL_Exception("ERROR searching: no goal state is set");
    }

    //goal state
    searchgoalstate = (RRTState*)(pSearchStateSpace->searchgoalstate->PlannerSpecificData);
    if (searchgoalstate->callnumberaccessed != pSearchStateSpace->callnumber) {
        ReInitializeSearchStateInfo(searchgoalstate, pSearchStateSpace);
    }

    //set goal key
    goalkey.key[0] = searchgoalstate->g;
    //goalkey.key[1] = searchgoalstate->h;

    //expand states until done
    minkey = pSearchStateSpace->heap->getminkeyheap();
    CKey oldkey = minkey;
    while (!pSearchStateSpace->heap->emptyheap() && minkey.key[0] < INFINITECOST && goalkey > minkey &&
           (clock() - TimeStarted) < MaxNumofSecs * (double)CLOCKS_PER_SEC &&
               (pSearchStateSpace->eps_satisfied == INFINITECOST ||
               (clock() - TimeStarted) < repair_time * (double)CLOCKS_PER_SEC))
    {
        //get the state
        state = (RRTState*)pSearchStateSpace->heap->deleteminheap();

#if DEBUG
        SBPL_FPRINTF(fDeb, "expanding state(%d): h=%d g=%u key=%u v=%u iterc=%d callnuma=%d expands=%d (g(goal)=%u)\n",
                     state->MDPstate->StateID, state->h, state->g, state->g+(int)(pSearchStateSpace->eps*state->h),
                     state->v, state->iterationclosed, state->callnumberaccessed, state->numofexpands,
                     searchgoalstate->g);
        SBPL_FPRINTF(fDeb, "expanding: ");
        PrintSearchState(state, fDeb);
        if (state->listelem[ARA_INCONS_LIST_ID] != NULL) {
            SBPL_FPRINTF(fDeb, "ERROR: expanding a state from inconslist\n");
            throw SBPL_Exception("ERROR: expanding a state from inconslist");
        }
        //SBPL_FFLUSH(fDeb);
#endif

#if DEBUG
        if (minkey.key[0] < oldkey.key[0] && fabs(this->finitial_eps - 1.0) < ERR_EPS) {
            //throw SBPL_Exception("WARN in search: the sequence of keys decreases");
        }
        oldkey = minkey;
#endif

        if (state->v == state->g) {
            SBPL_ERROR("ERROR: consistent state is being expanded\n");
#if DEBUG
            SBPL_FPRINTF(fDeb, "ERROR: consistent state is being expanded\n");
            throw SBPL_Exception("ERROR: consistent state is being expanded");
#endif
        }

        //recompute state value
        state->v = state->g;
        state->iterationclosed = pSearchStateSpace->searchiteration;

        //new expand
        expands++;
#if DEBUG
        state->numofexpands++;
#endif

        if (bforwardsearch == false)
            UpdatePreds(state, pSearchStateSpace);
        else
            UpdateSuccs(state, pSearchStateSpace);

        //recompute minkey
        minkey = pSearchStateSpace->heap->getminkeyheap();

        //recompute goalkey if necessary
        if (goalkey.key[0] != (int)searchgoalstate->g) {
            //SBPL_PRINTF("re-computing goal key\n");
            //recompute the goal key (heuristics should be zero)
            goalkey.key[0] = searchgoalstate->g;
            //goalkey.key[1] = searchgoalstate->h;
        }

        if (expands % 100000 == 0 && expands > 0) {
            SBPL_PRINTF("expands so far=%u\n", expands);
        }
    }

    int retv = 1;
    if (searchgoalstate->g == INFINITECOST && pSearchStateSpace->heap->emptyheap()) {
        SBPL_PRINTF("solution does not exist: search exited because heap is empty\n");
        retv = 0;
    }
    else if (!pSearchStateSpace->heap->emptyheap() && goalkey > minkey) {
        SBPL_PRINTF("search exited because it ran out of time\n");
        retv = 2;
    }
    else if (searchgoalstate->g == INFINITECOST && !pSearchStateSpace->heap->emptyheap()) {
        SBPL_PRINTF("solution does not exist: search exited because all candidates for expansion have "
                    "infinite heuristics\n");
        retv = 0;
    }
    else {
        SBPL_PRINTF("search exited with a solution for eps=%.3f\n", pSearchStateSpace->eps);
        retv = 1;
    }

    //SBPL_FPRINTF(fDeb, "expanded=%d\n", expands);

    searchexpands += expands;

    return retv;
}

void RRTPlanner::BuildNewOPENList(RRTSearchStateSpace_t* pSearchStateSpace)
{
    RRTState *state;
    CKey key;
    CHeap* pheap = pSearchStateSpace->heap;
    CList* pinconslist = pSearchStateSpace->inconslist;

    //move incons into open
    while (pinconslist->firstelement != NULL) {
        state = (RRTState*)pinconslist->firstelement->liststate;

        //compute f-value
        key.key[0] = state->g + (int)(pSearchStateSpace->eps * state->h);
        //key.key[1] = state->h;

        //insert into OPEN
        pheap->insertheap(state, key);
        //remove from INCONS
        pinconslist->remove(state, ARA_INCONS_LIST_ID);
    }
}

void RRTPlanner::Reevaluatefvals(RRTSearchStateSpace_t* pSearchStateSpace)
{
    CKey key;
    int i;
    CHeap* pheap = pSearchStateSpace->heap;

    //recompute priorities for states in OPEN and reorder it
    for (i = 1; i <= pheap->currentsize; ++i) {
        RRTState* state = (RRTState*)pheap->heap[i].heapstate;
        pheap->heap[i].key.key[0] = state->g + (int)(pSearchStateSpace->eps * state->h);
        //pheap->heap[i].key.key[1] = state->h;
    }
    pheap->makeheap();

    pSearchStateSpace->bReevaluatefvals = false;
}

void RRTPlanner::Reevaluatehvals(RRTSearchStateSpace_t* pSearchStateSpace)
{
    for(int i = 0; i < (int)pSearchStateSpace->searchMDP.StateArray.size(); i++)
    {
        CMDPSTATE* MDPstate = pSearchStateSpace->searchMDP.StateArray[i];
        RRTState* state = (RRTState*)MDPstate->PlannerSpecificData;
        state->h = ComputeHeuristic(MDPstate, pSearchStateSpace);
    }
}

//creates (allocates memory) search state space
//does not initialize search statespace
int RRTPlanner::CreateSearchStateSpace(RRTSearchStateSpace_t* pSearchStateSpace)
{
    //create a heap
    pSearchStateSpace->heap = new CHeap;
    pSearchStateSpace->inconslist = new CList;
    MaxMemoryCounter += sizeof(CHeap);
    MaxMemoryCounter += sizeof(CList);

    pSearchStateSpace->searchgoalstate = NULL;
    pSearchStateSpace->searchstartstate = NULL;

    searchexpands = 0;
    num_of_expands_initial_solution = -1;

    pSearchStateSpace->bReinitializeSearchStateSpace = false;

    return 1;
}

//deallocates memory used by SearchStateSpace
void RRTPlanner::DeleteSearchStateSpace(RRTSearchStateSpace_t* pSearchStateSpace)
{
    if (pSearchStateSpace->heap != NULL) {
        pSearchStateSpace->heap->makeemptyheap();
        delete pSearchStateSpace->heap;
        pSearchStateSpace->heap = NULL;
    }

    if (pSearchStateSpace->inconslist != NULL) {
        pSearchStateSpace->inconslist->makeemptylist(ARA_INCONS_LIST_ID);
        delete pSearchStateSpace->inconslist;
        pSearchStateSpace->inconslist = NULL;
    }

    //delete the states themselves
    int iend = (int)pSearchStateSpace->searchMDP.StateArray.size();
    for (int i = 0; i < iend; i++) {
        CMDPSTATE* state = pSearchStateSpace->searchMDP.StateArray[i];
        if (state != NULL && state->PlannerSpecificData != NULL) {
            DeleteSearchStateData((RRTState*)state->PlannerSpecificData);
            free((RRTState*)state->PlannerSpecificData);
            state->PlannerSpecificData = NULL;
        }
    }
    pSearchStateSpace->searchMDP.Delete();
}

//reset properly search state space
//needs to be done before deleting states
int RRTPlanner::ResetSearchStateSpace(RRTSearchStateSpace_t* pSearchStateSpace)
{
    pSearchStateSpace->heap->makeemptyheap();
    pSearchStateSpace->inconslist->makeemptylist(ARA_INCONS_LIST_ID);

    return 1;
}

//initialization before each search
void RRTPlanner::ReInitializeSearchStateSpace(RRTSearchStateSpace_t* pSearchStateSpace)
{
    CKey key;

    //increase callnumber
    pSearchStateSpace->callnumber++;

    //reset iteration
    pSearchStateSpace->searchiteration = 0;
    pSearchStateSpace->bNewSearchIteration = true;

#if DEBUG
    SBPL_FPRINTF(fDeb, "reinitializing search state-space (new call number=%d search iter=%d)\n",
        pSearchStateSpace->callnumber,pSearchStateSpace->searchiteration );
#endif

    pSearchStateSpace->heap->makeemptyheap();
    pSearchStateSpace->inconslist->makeemptylist(ARA_INCONS_LIST_ID);

    //reset
    pSearchStateSpace->eps = this->finitial_eps;
    pSearchStateSpace->eps_satisfied = INFINITECOST;

    //initialize start state
    RRTState* startstateinfo = (RRTState*)(pSearchStateSpace->searchstartstate->PlannerSpecificData);
    if (startstateinfo->callnumberaccessed != pSearchStateSpace->callnumber) {
        ReInitializeSearchStateInfo(startstateinfo, pSearchStateSpace);
    }
    startstateinfo->g = 0;

    //initialize goal state
    RRTState* searchgoalstate = (RRTState*)(pSearchStateSpace->searchgoalstate->PlannerSpecificData);
    if (searchgoalstate->callnumberaccessed != pSearchStateSpace->callnumber) {
        ReInitializeSearchStateInfo(searchgoalstate, pSearchStateSpace);
    }

    //insert start state into the heap
    key.key[0] = (long int)(pSearchStateSpace->eps * startstateinfo->h);
    //key.key[1] = startstateinfo->h;
    pSearchStateSpace->heap->insertheap(startstateinfo, key);

    pSearchStateSpace->bReinitializeSearchStateSpace = false;
    pSearchStateSpace->bReevaluatefvals = false;
}

//very first initialization
int RRTPlanner::InitializeSearchStateSpace(RRTSearchStateSpace_t* pSearchStateSpace)
{
    if (pSearchStateSpace->heap->currentsize != 0 || pSearchStateSpace->inconslist->currentsize != 0) {
        throw SBPL_Exception("ERROR in InitializeSearchStateSpace: heap or list is not empty");
    }

    pSearchStateSpace->eps = this->finitial_eps;
    pSearchStateSpace->eps_satisfied = INFINITECOST;
    pSearchStateSpace->searchiteration = 0;
    pSearchStateSpace->bNewSearchIteration = true;
    pSearchStateSpace->callnumber = 0;
    pSearchStateSpace->bReevaluatefvals = false;

    //create and set the search start state
    pSearchStateSpace->searchgoalstate = NULL;
    //pSearchStateSpace->searchstartstate = GetState(SearchStartStateID, pSearchStateSpace);
    pSearchStateSpace->searchstartstate = NULL;

    pSearchStateSpace->bReinitializeSearchStateSpace = true;

    return 1;
}

int RRTPlanner::SetSearchGoalState(int SearchGoalStateID, RRTSearchStateSpace_t* pSearchStateSpace)
{
    if (pSearchStateSpace->searchgoalstate == NULL ||
        pSearchStateSpace->searchgoalstate->StateID != SearchGoalStateID)
    {
        pSearchStateSpace->searchgoalstate = GetState(SearchGoalStateID, pSearchStateSpace);

        //should be new search iteration
        pSearchStateSpace->eps_satisfied = INFINITECOST;
        pSearchStateSpace->bNewSearchIteration = true;
        pSearchStateSpace_->eps = this->finitial_eps;

#if USE_HEUR
        //recompute heuristic for the heap if heuristics is used
        pSearchStateSpace->bReevaluatefvals = true;
#endif
    }

    return 1;
}

int RRTPlanner::SetSearchStartState(int SearchStartStateID, RRTSearchStateSpace_t* pSearchStateSpace)
{
    CMDPSTATE* MDPstate = GetState(SearchStartStateID, pSearchStateSpace);

    if (MDPstate != pSearchStateSpace->searchstartstate) {
        pSearchStateSpace->searchstartstate = MDPstate;
        pSearchStateSpace->bReinitializeSearchStateSpace = true;
    }

    return 1;
}

int RRTPlanner::ReconstructPath(RRTSearchStateSpace_t* pSearchStateSpace)
{
    if (bforwardsearch) //nothing to do, if search is backward
    {
        CMDPSTATE* MDPstate = pSearchStateSpace->searchgoalstate;
        CMDPSTATE* PredMDPstate;
        RRTState *predstateinfo, *stateinfo;

#if DEBUG
        SBPL_FPRINTF(fDeb, "reconstructing a path:\n");
#endif

        while (MDPstate != pSearchStateSpace->searchstartstate) {
            stateinfo = (RRTState*)MDPstate->PlannerSpecificData;

#if DEBUG
            PrintSearchState(stateinfo, fDeb);
#endif
            if (stateinfo->g == INFINITECOST) {
                //throw SBPL_Exception("ERROR in ReconstructPath: g of the state on the path is INFINITE");
                return -1;
            }

            if (stateinfo->bestpredstate == NULL) {
                SBPL_ERROR("ERROR in ReconstructPath: bestpred is NULL\n");
                throw SBPL_Exception("ERROR in ReconstructPath: bestpred is NULL");
            }

            //get the parent state
            PredMDPstate = stateinfo->bestpredstate;
            predstateinfo = (RRTState*)PredMDPstate->PlannerSpecificData;

            //set its best next info
            predstateinfo->bestnextstate = MDPstate;

            //check the decrease of g-values along the path
            if (predstateinfo->v >= stateinfo->g) {
                SBPL_ERROR("ERROR in ReconstructPath: g-values are non-decreasing\n");
                PrintSearchState(predstateinfo, fDeb);
                throw SBPL_Exception("ERROR in ReconstructPath: g-values are non-decreasing");
            }

            //transition back
            MDPstate = PredMDPstate;
        }
    }

    return 1;
}

void RRTPlanner::PrintSearchPath(RRTSearchStateSpace_t* pSearchStateSpace, FILE* fOut)
{
    RRTState* searchstateinfo;
    CMDPSTATE* state, *gstate;
    int goalID;
    int PathCost;

    if (bforwardsearch) {
        state  = pSearchStateSpace->searchstartstate;
        gstate = pSearchStateSpace->searchgoalstate;
    }
    else {
        state  = pSearchStateSpace->searchgoalstate;
        gstate = pSearchStateSpace->searchstartstate;
    }
    goalID = gstate->StateID;

    if (fOut == NULL) fOut = stdout;

    PathCost = ((RRTState*)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g;

    SBPL_FPRINTF(fOut, "\n\nPrinting a path from the state %d ", gstate->StateID);
    environment_->PrintState(gstate->StateID, true, fOut);
    SBPL_FPRINTF(fOut, " to the state %d ", state->StateID);
    environment_->PrintState(state->StateID, true, fOut);
	SBPL_FPRINTF(fOut, "\nEpsilon   = %.2f\n", pSearchStateSpace->eps);
    SBPL_FPRINTF(fOut, "Path cost = %d\n\n", PathCost);

    // int costFromStart = 0;
    while (state->StateID != goalID) {
        SBPL_FPRINTF(fOut, "state %4d ", state->StateID);
        environment_->PrintState(state->StateID, true, fOut);

        if (state->PlannerSpecificData == NULL) {
            SBPL_FPRINTF(fOut, "path does not exist since search data does not exist\n");
            break;
        }

        searchstateinfo = (RRTState*)state->PlannerSpecificData;

        if (searchstateinfo->bestnextstate == NULL) {
            SBPL_FPRINTF(fOut, "path does not exist since bestnextstate == NULL\n");
            break;
        }
        if (searchstateinfo->g == INFINITECOST) {
            SBPL_FPRINTF(fOut, "path does not exist since bestnextstate == NULL\n");
            break;
        }

        // int costToGoal = PathCost - costFromStart;
        // int transcost = searchstateinfo->g - ((RRTState*)(searchstateinfo->bestnextstate->PlannerSpecificData))->v;
        // if (bforwardsearch) transcost = -transcost;
        // costFromStart += transcost;

        SBPL_FPRINTF(fOut, " <-- g=%8d e*h=%11.2f <-- state %4d ", searchstateinfo->g, searchstateinfo->h * pSearchStateSpace->eps,
                     searchstateinfo->bestnextstate->StateID);

        state = searchstateinfo->bestnextstate;

        environment_->PrintState(state->StateID, true, fOut); SBPL_FPRINTF(fOut, "\n");
    }
    SBPL_FPRINTF(fOut, "\n\n");
}

void RRTPlanner::PrintSearchState(RRTState* state, FILE* fOut)
{
#if DEBUG
    SBPL_FPRINTF(fOut, "state %d: h=%d g=%u v=%u iterc=%d callnuma=%d expands=%d heapind=%d inconslist=%d\n",
                 state->MDPstate->StateID, state->h, state->g, state->v,
                 state->iterationclosed, state->callnumberaccessed, state->numofexpands, state->heapindex,
                 state->listelem[ARA_INCONS_LIST_ID] ? 1 : 0);
#else
    SBPL_FPRINTF(fOut, "state %d: h=%d g=%u v=%u iterc=%d callnuma=%d heapind=%d inconslist=%d\n",
                 state->MDPstate->StateID, state->h, state->g, state->v, state->iterationclosed,
                 state->callnumberaccessed, state->heapindex, state->listelem[ARA_INCONS_LIST_ID] ? 1 : 0);
#endif
    environment_->PrintState(state->MDPstate->StateID, true, fOut);
}

int RRTPlanner::getHeurValue(RRTSearchStateSpace_t* pSearchStateSpace, int StateID)
{
    CMDPSTATE* MDPstate = GetState(StateID, pSearchStateSpace);
    RRTState* searchstateinfo = (RRTState*)MDPstate->PlannerSpecificData;
    return searchstateinfo->h;
}

vector<int> RRTPlanner::GetSearchPath(RRTSearchStateSpace_t* pSearchStateSpace, int& solcost)
{
    vector<int> SuccIDV;
    vector<int> CostV;
    vector<int> wholePathIds;
    RRTState* searchstateinfo;
    CMDPSTATE* state = NULL;
    CMDPSTATE* goalstate = NULL;
    CMDPSTATE* startstate = NULL;

    if (bforwardsearch) {
        startstate = pSearchStateSpace->searchstartstate;
        goalstate = pSearchStateSpace->searchgoalstate;

        //reconstruct the path by setting bestnextstate pointers appropriately
        ReconstructPath(pSearchStateSpace);
    }
    else {
        startstate = pSearchStateSpace->searchgoalstate;
        goalstate = pSearchStateSpace->searchstartstate;
    }

    state = startstate;

    wholePathIds.push_back(state->StateID);
    solcost = 0;

    FILE* fOut = stdout;
    if (fOut == NULL) {
        throw SBPL_Exception("ERROR: could not open file");
    }
    while (state->StateID != goalstate->StateID) {
        if (state->PlannerSpecificData == NULL) {
            SBPL_FPRINTF(fOut, "path does not exist since search data does not exist\n");
            break;
        }

        searchstateinfo = (RRTState*)state->PlannerSpecificData;

        if (searchstateinfo->bestnextstate == NULL) {
            SBPL_FPRINTF(fOut, "path does not exist since bestnextstate == NULL\n");
            break;
        }
        if (searchstateinfo->g == INFINITECOST) {
            SBPL_FPRINTF(fOut, "path does not exist since bestnextstate == NULL\n");
            break;
        }

        environment_->GetSuccs(state->StateID, &SuccIDV, &CostV);
        int actioncost = INFINITECOST;
        for (int i = 0; i < (int)SuccIDV.size(); i++) {
            if (SuccIDV.at(i) == searchstateinfo->bestnextstate->StateID && CostV.at(i) < actioncost) {
                actioncost = CostV.at(i);
            }
        }
        if (actioncost == INFINITECOST) SBPL_PRINTF("WARNING: actioncost = %d\n", actioncost);

        solcost += actioncost;

        //SBPL_FPRINTF(fDeb, "actioncost=%d between states %d and %d\n",
        //        actioncost, state->StateID, searchstateinfo->bestnextstate->StateID);
        //environment_->PrintState(state->StateID, false, fDeb);
        //environment_->PrintState(searchstateinfo->bestnextstate->StateID, false, fDeb);

#if DEBUG
        RRTState* nextstateinfo = (RRTState*)(searchstateinfo->bestnextstate->PlannerSpecificData);
        if (actioncost != abs((int)(searchstateinfo->g - nextstateinfo->g)) &&
            pSearchStateSpace->eps_satisfied <= 1.001)
        {
            SBPL_FPRINTF(fDeb, "ERROR: actioncost=%d is not matching the difference in g-values of %d\n",
                         actioncost, abs((int)(searchstateinfo->g - nextstateinfo->g)));
            SBPL_ERROR("ERROR: actioncost=%d is not matching the difference in g-values of %d\n",
                       actioncost,abs((int)(searchstateinfo->g - nextstateinfo->g)));
            PrintSearchState(searchstateinfo, fDeb);
            PrintSearchState(nextstateinfo, fDeb);
        }
#endif

        state = searchstateinfo->bestnextstate;

        wholePathIds.push_back(state->StateID);
    }

    return wholePathIds;
}

bool RRTPlanner::Search(RRTSearchStateSpace_t* pSearchStateSpace, vector<int>& pathIds, int & PathCost,
                        bool bFirstSolution, bool bOptimalSolution, double MaxNumofSecs)
{
    CKey key;
    TimeStarted = clock();
    searchexpands = 0;
    num_of_expands_initial_solution = -1;
    double old_repair_time = repair_time;
    if (!use_repair_time)
        repair_time = MaxNumofSecs;

#if DEBUG
    SBPL_FPRINTF(fDeb, "new search call (call number=%d)\n", pSearchStateSpace->callnumber);
#endif

    if (pSearchStateSpace->bReevaluatefvals) {
        // costs have changed or a new goal has been set
        environment_->EnsureHeuristicsUpdated(bforwardsearch);
        Reevaluatehvals(pSearchStateSpace);
    }

    if (pSearchStateSpace->bReinitializeSearchStateSpace) {
        //re-initialize state space
        ReInitializeSearchStateSpace(pSearchStateSpace);
    }

    if (bOptimalSolution) {
        pSearchStateSpace->eps = 1;
        MaxNumofSecs = INFINITECOST;
        repair_time = INFINITECOST;
    }
    else if (bFirstSolution) {
        MaxNumofSecs = INFINITECOST;
        repair_time = INFINITECOST;
    }

    //the main loop of ARA*
    stats.clear();
    int prevexpands = 0;
    clock_t loop_time;
    while (pSearchStateSpace->eps_satisfied > final_epsilon &&
           (clock() - TimeStarted) < MaxNumofSecs * (double)CLOCKS_PER_SEC &&
               (pSearchStateSpace->eps_satisfied == INFINITECOST ||
               (clock() - TimeStarted) < repair_time * (double)CLOCKS_PER_SEC))
    {
        loop_time = clock();
        //decrease eps for all subsequent iterations
        if (fabs(pSearchStateSpace->eps_satisfied - pSearchStateSpace->eps) < ERR_EPS && !bFirstSolution) {
            pSearchStateSpace->eps = pSearchStateSpace->eps - dec_eps;
            if (pSearchStateSpace->eps < final_epsilon)
                pSearchStateSpace->eps = final_epsilon;

            //the priorities need to be updated
            pSearchStateSpace->bReevaluatefvals = true;

            //it will be a new search
            pSearchStateSpace->bNewSearchIteration = true;
        }

        if (pSearchStateSpace->bNewSearchIteration) {
            pSearchStateSpace->searchiteration++;
            pSearchStateSpace->bNewSearchIteration = false;
            BuildNewOPENList(pSearchStateSpace);
        }

        //re-compute f-values if necessary and reorder the heap
        if (pSearchStateSpace->bReevaluatefvals)
            Reevaluatefvals(pSearchStateSpace);

        //improve or compute path
        if (ImprovePath(pSearchStateSpace, MaxNumofSecs) == 1) {
            pSearchStateSpace->eps_satisfied = pSearchStateSpace->eps;
        }

        //print the solution cost and eps bound
        SBPL_PRINTF("eps=%f expands=%d g(searchgoal)=%d time=%.3f\n\n", pSearchStateSpace->eps_satisfied,
                    searchexpands - prevexpands,
                    ((RRTState*)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g,
                    double(clock() - loop_time) / CLOCKS_PER_SEC);

        if (pSearchStateSpace->eps_satisfied == finitial_eps && pSearchStateSpace->eps == finitial_eps) {
            finitial_eps_planning_time = double(clock() - loop_time) / CLOCKS_PER_SEC;
            num_of_expands_initial_solution = searchexpands - prevexpands;
        }

        if (stats.empty() || pSearchStateSpace->eps_satisfied != stats.back().eps) {
            PlannerStats tempStat;
            tempStat.eps = pSearchStateSpace->eps_satisfied;
            tempStat.expands = searchexpands - prevexpands;
            tempStat.time = double(clock() - loop_time) / CLOCKS_PER_SEC;
            tempStat.cost = ((RRTState*)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g;
            stats.push_back(tempStat);
        }

#if DEBUG
        SBPL_FPRINTF(fDeb, "eps=%f expands=%d g(searchgoal)=%d time=%.3f\n", pSearchStateSpace->eps_satisfied,
                     searchexpands - prevexpands,
                     ((RRTState*)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g,
                     double(clock()-loop_time)/CLOCKS_PER_SEC);
        PrintSearchState((RRTState*)pSearchStateSpace->searchgoalstate->PlannerSpecificData, fDeb);
        print_searchpath(fDeb);
#endif
        prevexpands = searchexpands;

        //if just the first solution then we are done
        if (bFirstSolution)
            break;

        //no solution exists
        if (((RRTState*)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g == INFINITECOST)
            break;
    }
    repair_time = old_repair_time;

#if DEBUG
    SBPL_FFLUSH(fDeb);
#endif

    PathCost = ((RRTState*)pSearchStateSpace->searchgoalstate->PlannerSpecificData)->g;
    MaxMemoryCounter += environment_->StateID2IndexMapping.size() * sizeof(int);

    // SBPL_PRINTF("MaxMemoryCounter = %d\n", MaxMemoryCounter);

    int solcost = INFINITECOST;
    bool ret = false;
    if (PathCost == INFINITECOST) {
        // SBPL_PRINTF("could not find a solution\n");
        ret = false;
    }
    else {
        // SBPL_PRINTF("solution is found\n");
        pathIds = GetSearchPath(pSearchStateSpace, solcost);
        ret = true;
    }

    // SBPL_PRINTF("total expands this call = %d, planning time = %.3f secs, solution cost=%d\n",
    //             searchexpands, (clock() - TimeStarted) / ((double)CLOCKS_PER_SEC), solcost);
    final_eps_planning_time = (clock() - TimeStarted) / ((double)CLOCKS_PER_SEC);
    final_eps = pSearchStateSpace->eps_satisfied;
    //SBPL_FPRINTF(fStat, "%d %d\n", searchexpands, solcost);

    return ret;
}

//-----------------------------Interface function-----------------------------------------------------

int RRTPlanner::replan(vector<int>* solution_stateIDs_V, ReplanParams params)
{
    int solcost;
    return replan(solution_stateIDs_V, params, &solcost);
}

// int RRTPlanner::replan(vector<int>* solution_stateIDs_V, ReplanParams params, int* solcost)
// {
//    finitial_eps              = params.initial_eps;
//    final_epsilon             = params.final_eps;
//    dec_eps                   = params.dec_eps;
//    bsearchuntilfirstsolution = params.return_first_solution;
//    use_repair_time           = params.repair_time > 0;
//    repair_time               = params.repair_time;
//    return replan(params.max_time, solution_stateIDs_V, solcost);
// }

//returns 1 if found a solution, and 0 otherwise
int RRTPlanner::replan(double allocated_time_secs, vector<int>* solution_stateIDs_V)
{
    int solcost;

    return replan(allocated_time_secs, solution_stateIDs_V, &solcost);
}

//returns 1 if found a solution, and 0 otherwise
int RRTPlanner::replan(double allocated_time_secs, vector<int>* solution_stateIDs_V, int* psolcost)
{
    vector<int> pathIds;
    bool bFound = false;
    int PathCost;
    bool bFirstSolution = this->bsearchuntilfirstsolution;
    bool bOptimalSolution = false;
    *psolcost = 0;

    SBPL_PRINTF("planner: replan called (bFirstSol=%d)\n\n", bFirstSolution);

    //plan
    bFound = Search(pSearchStateSpace_, pathIds, PathCost, bFirstSolution, bOptimalSolution, allocated_time_secs);
    if (!bFound)
    {
        SBPL_PRINTF("failed to find a solution\n");
    }

    //copy the solution
    *solution_stateIDs_V = pathIds;
    *psolcost = PathCost;

    return (int)bFound;
}

int RRTPlanner::set_goal(int goal_stateID)
{
    SBPL_PRINTF("planner: setting goal to %d\n", goal_stateID);
    environment_->PrintState(goal_stateID, true, stdout);

    if (bforwardsearch) {
        if (SetSearchGoalState(goal_stateID, pSearchStateSpace_) != 1) {
            SBPL_ERROR("ERROR: failed to set search goal state\n");
            return 0;
        }
    }
    else {
        if (SetSearchStartState(goal_stateID, pSearchStateSpace_) != 1) {
            SBPL_ERROR("ERROR: failed to set search start state\n");
            return 0;
        }
    }

    return 1;
}

int RRTPlanner::set_start(int start_stateID)
{
    SBPL_PRINTF("planner: setting start to %d\n", start_stateID);
    environment_->PrintState(start_stateID, true, stdout);

    if (bforwardsearch) {
        if (SetSearchStartState(start_stateID, pSearchStateSpace_) != 1) {
            SBPL_ERROR("ERROR: failed to set search start state\n");
            return 0;
        }
    }
    else {
        if (SetSearchGoalState(start_stateID, pSearchStateSpace_) != 1) {
            SBPL_ERROR("ERROR: failed to set search goal state\n");
            return 0;
        }
    }

    return 1;
}

void RRTPlanner::costs_changed(StateChangeQuery const & stateChange)
{
    pSearchStateSpace_->bReevaluatefvals = true;
    pSearchStateSpace_->bReinitializeSearchStateSpace = true;
}

void RRTPlanner::costs_changed()
{
    pSearchStateSpace_->bReevaluatefvals = true;
    pSearchStateSpace_->bReinitializeSearchStateSpace = true;
}

int RRTPlanner::force_planning_from_scratch()
{
    SBPL_PRINTF("planner: forceplanfromscratch set\n");

    pSearchStateSpace_->bReinitializeSearchStateSpace = true;

    return 1;
}

int RRTPlanner::force_planning_from_scratch_and_free_memory()
{
    SBPL_PRINTF("planner: forceplanfromscratch set\n");
    int start_id = -1;
    int goal_id = -1;
    if (pSearchStateSpace_->searchstartstate)
        start_id = pSearchStateSpace_->searchstartstate->StateID;
    if (pSearchStateSpace_->searchgoalstate)
        goal_id = pSearchStateSpace_->searchgoalstate->StateID;

    if (!bforwardsearch) {
        int temp = start_id;
        start_id = goal_id;
        goal_id = temp;
    }

    DeleteSearchStateSpace(pSearchStateSpace_);
    CreateSearchStateSpace(pSearchStateSpace_);
    InitializeSearchStateSpace(pSearchStateSpace_);
    for (unsigned int i = 0; i < environment_->StateID2IndexMapping.size(); i++)
        for (int j = 0; j < NUMOFINDICES_STATEID2IND; j++)
            environment_->StateID2IndexMapping[i][j] = -1;

    if (start_id >= 0)
        set_start(start_id);
    if (goal_id >= 0)
        set_goal(goal_id);
    return 1;
}

int RRTPlanner::set_search_mode(bool bSearchUntilFirstSolution)
{
    SBPL_PRINTF("planner: search mode set to %d\n", bSearchUntilFirstSolution);

    bsearchuntilfirstsolution = bSearchUntilFirstSolution;

    return 1;
}

void RRTPlanner::print_searchpath(FILE* fOut)
{
    PrintSearchPath(pSearchStateSpace_, fOut);
}

double RRTPlanner::compute_suboptimality()
{
    if (!pSearchStateSpace_) {
        SBPL_ERROR("Search state space is NULL in compute_suboptimality()");
        return -1.0;
    }

    // find the min cost from all paths that go through states in the incons
    // list
    int inconsListMin = std::numeric_limits<int>::max();
    if (pSearchStateSpace_->inconslist) {
        AbstractSearchState* currElem =
                pSearchStateSpace_->inconslist->getfirst();
        while (currElem) {
            // get the stateID of this state
            RRTState* state = (RRTState*)currElem;
            assert(state);

            // int stateID = state->MDPstate->StateID;

            // update the min to the f-value of this state if necessary
            int h = state->h;
            int g = state->g;
            if (g + h < inconsListMin) {
                inconsListMin = g + h;
            }

            currElem = pSearchStateSpace_->inconslist->getnext(
                    currElem, ARAMDP_STATEID2IND);
        }
    }

    // find the min cost from all paths that go through states still in the open
    // list
    int openListMin = std::numeric_limits<int>::max();
    if (pSearchStateSpace_->heap) {
        for (int i = 1; i < pSearchStateSpace_->heap->currentsize; i++) {
            AbstractSearchState* abstractState =
                    pSearchStateSpace_->heap->heap[i].heapstate;
            if (!abstractState) {
                SBPL_ERROR("heap element with keys %d and %d has NULL AbstractSearchState\n", (int)pSearchStateSpace_->heap->heap[i].key[0], (int)pSearchStateSpace_->heap->heap[i].key[1]);
                continue; // return -1.0 ?
            }

            RRTState* state2 = (RRTState*)abstractState;

            int h = state2->h;
            int g = state2->g;
            if (g + h < openListMin) {
                openListMin = g + h;
            }
        }
    }

    SBPL_DEBUG("Done looking through INCONS and OPEN lists for state with minimum f-value\n");

    int overallMin = openListMin < inconsListMin ? openListMin : inconsListMin;
    SBPL_DEBUG("f_min = min(f_open_min = %d, f_incons_min = %d) = %d\n", openListMin, inconsListMin, overallMin);

    if (overallMin == std::numeric_limits<int>::max()) {
        SBPL_ERROR("Couldn't find a min f-value. Empty incons list or open list.\n");
        return -1.0;
    }

    if (!pSearchStateSpace_->searchgoalstate) {
        return -1.0;
    }

    int goalGValue = GetGVal(
            pSearchStateSpace_->searchgoalstate->StateID, pSearchStateSpace_);
    double lowerBound = std::numeric_limits<double>::max();
    if (overallMin != 0) {
        lowerBound = double(goalGValue) / double(overallMin);
    }

    SBPL_DEBUG("Lower Bound = %d / %d = %0.3f", goalGValue, overallMin, lowerBound);
    SBPL_DEBUG("Eps Satisfied = %0.3f", pSearchStateSpace_->eps_satisfied);

    return std::max(1.0, std::min(pSearchStateSpace_->eps_satisfied, lowerBound));
}

void RRTPlanner::get_search_stats(vector<PlannerStats>* s)
{
    s->clear();
    s->reserve(stats.size());
    for (unsigned int i = 0; i < stats.size(); i++) {
        s->push_back(stats[i]);
    }
}
