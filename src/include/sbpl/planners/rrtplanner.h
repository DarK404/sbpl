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
 */

#ifndef __RRTPLANNER_H_
#define __RRTPLANNER_H_

#include <cstdio>
#include <ctime>
#include <vector>
#include <ostream>
#include <sbpl/planners/planner.h>
#include <sbpl/utils/mdp.h>

//---configuration----
//control of EPS
//initial suboptimality bound (cost solution <= cost(eps*cost optimal solution)
#define ARA_DEFAULT_INITIAL_EPS	    5.0
//as planning time exist, ARA* decreases epsilon bound
#define ARA_DECREASE_EPS    0.2
//final epsilon bound
#define ARA_FINAL_EPS	    1.0
//---------------------
#define ARA_INCONS_LIST_ID 0

class CHeap;
class CList;
class DiscreteSpaceInformation;
class MDPConfig;
class StateChangeQuery;

enum errt_TreeType
{
    INVALID_TREE_TYPE = -1, SINGLE_TREE_TYPE, DOUBLE_TREE_TYPE,

    NUM_TREE_TYPES
};

enum errt_PathType
{
    INVALID_PATH_TYPE = -1, NOTSMOOTHED_PATH_TYPE, SMOOTHED_PATH_TYPE,

    NUM_PATH_TYPES
};

//-------------------------------------------------------------

/**
 * \brief Node of a Tree
 */

class RRTNode
{
	// Atributos
public:
	std::string aID;
	int aX;
	int aY;
	
	// Métodos
public:
	RRTNode() { }

	RRTNode(std::string ID, int X, int Y) {
		aID = ID;
		aX  = X;
		aY  = Y;
	}

	friend std::ostream& operator << (std::ostream& os, const RRTNode& node)
	{
		os << node.aID << "<<" << node.aX << "," << node.aY << ">>";
		return os;
	}
};

/**
 * \brief Coordinates of a Map Cell
 */
class RRTMapCoord
{
public:
    int attr_x;
	int attr_y;

	RRTMapCoord(int x, int y) {
		attr_x = x;
		attr_y = y;
	}
};

/**
 * \brief state structure used in ARA* search tree
 */
typedef class RRTSEARCHSTATEDATA : public AbstractSearchState
{
public:
    /**
     * \brief the MDP state itself
     */
    CMDPSTATE* MDPstate;
    /**
     * \brief ARA* relevant data
     */
    unsigned int v;
    /**
     * \brief ARA* relevant data
     */
    unsigned int g;
    /**
     * \brief ARA* relevant data
     */
    short unsigned int iterationclosed;
    /**
     * \brief ARA* relevant data
     */
    short unsigned int callnumberaccessed;

#if DEBUG
    /**
     * \brief ARA* relevant data
     */
    short unsigned int numofexpands;
#endif

    /**
     * \brief best predecessor and the action from it, used only in forward searches
     */
    CMDPSTATE *bestpredstate;
    /**
     * \brief the next state if executing best action
     */
    CMDPSTATE *bestnextstate;
    unsigned int costtobestnextstate;
    int h;

public:
    RRTSEARCHSTATEDATA() { }
    ~RRTSEARCHSTATEDATA() { }
} RRTState;

/**
 * \brief the statespace of ARA*
 */
typedef struct RRTSEARCHSTATESPACE
{
    double eps;
    double eps_satisfied;
    CHeap* heap;
    CList* inconslist;
    short unsigned int searchiteration;
    short unsigned int callnumber;
    CMDPSTATE* searchgoalstate;
    CMDPSTATE* searchstartstate;

    CMDP searchMDP;

    bool bReevaluatefvals;
    bool bReinitializeSearchStateSpace;
    bool bNewSearchIteration;
} RRTSearchStateSpace_t;

/**
 * \brief RRT planner
 */
class RRTPlanner : public SBPLPlanner
{
public:
    /**
     * \brief replan a path within the allocated time, return the solution in the vector
     */
    virtual int replan(double allocated_time_secs, std::vector<int>* solution_stateIDs_V);

    /**
     * \brief replan a path within the allocated time, return the solution in the vector, also returns solution cost
     */
    virtual int replan(double allocated_time_sec, std::vector<int>* solution_stateIDs_V, int* solcost);

    /**
     * \brief works same as replan function with time and solution states, but
     *        it let's you fill out all the parameters for the search
     */
    virtual int replan(std::vector<int>* solution_stateIDs_V, ReplanParams params);

    /**
     * \brief works same as replan function with time, solution states, and
     *        cost, but it let's you fill out all the parameters for the search
     */
    virtual int replan(std::vector<int>* solution_stateIDs_V, ReplanParams params, int* solcost);

    /**
     * \brief set the goal state
     */
    virtual int set_goal(int goal_stateID);

    /**
     * \brief set the start state
     */
    virtual int set_start(int start_stateID);

    /**
     * \brief inform the search about the new edge costs
     */
    virtual void costs_changed(StateChangeQuery const & stateChange);

    /**
     * \brief inform the search about the new edge costs -
     * \note since ARA* is non-incremental, it is sufficient (and more
     *       efficient) to just inform ARA* of the fact that some costs changed
     */
    virtual void costs_changed();

    /**
     * \brief set a flag to get rid of the previous search efforts, and
     *        re-initialize the search, when the next replan is called
     */
    virtual int force_planning_from_scratch();

    /**
     * \brief Gets rid of the previous search efforts, release the memory and re-initialize the search.
     */
    virtual int force_planning_from_scratch_and_free_memory();

    /**
     * \brief you can either search forwards or backwards
     */
    virtual int set_search_mode(bool bSearchUntilFirstSolution);

    /**
     * \brief returns the suboptimality bound on the currently found solution
     */
    virtual double get_solution_eps() const { return pSearchStateSpace_->eps_satisfied; }

    /**
     * \brief returns the number of states expanded so far
     */
    virtual int get_n_expands() const { return searchexpands; }

    /**
     * \brief sets the value of the initial epsilon (suboptimality bound) used
     */
    virtual void set_initialsolution_eps(double initialsolution_eps) { finitial_eps = initialsolution_eps; }

    /**
     * \brief sets the value to decrease from eps at each iteration
     */
    virtual void set_eps_step(double eps) { dec_eps = eps; }

    /**
     * \brief prints out the search path into a file
     */
    virtual void print_searchpath(FILE* fOut);

    /**
     * \brief Compute the suboptimality bound for the most recent solution.
     *
     * The suboptimality bound of the solution may be provably less than the
     * value of epsilon satisfied during the most recent planning iteration.
     * This suboptimality bound is computed as the ratio between the current
     * g-value for the goal and the minimum un-weighted f-value of a locally
     * inconsistent state.
     * 
     * \return The suboptimality bound of the most recently computed solution
     */
    double compute_suboptimality();

    /**
     * \brief constructor
     */
	RRTPlanner(DiscreteSpaceInformation* environment, errt_TreeType treetype, errt_PathType pathtype,
	           double stepsize, double qgoalprob, int chaikiniter, double maxplanningtime);
    void prevConstructor(DiscreteSpaceInformation* environment, bool bforwardsearch);

    /**
     * \brief destructor
     */
    ~RRTPlanner();

    /**
     * \brief returns the initial epsilon
     */
    virtual double get_initial_eps() { return finitial_eps; }

    /**
     * \brief returns the time taken to find the first solution
     */
    virtual double get_initial_eps_planning_time() { return finitial_eps_planning_time; }

    /**
     * \brief returns the time taken to get the final solution
     */
    virtual double get_final_eps_planning_time() { return final_eps_planning_time; }

    /**
     * \brief Return the number of expands to find the first solution or -1 if no solution has been found.
     */
    virtual int get_n_expands_init_solution() { return num_of_expands_initial_solution; }

    /**
     * \brief returns the final epsilon achieved during the search
     */
    virtual double get_final_epsilon() { return final_eps; }

    /**
     * \brief fills out a vector of stats from the search
     */
    virtual void get_search_stats(std::vector<PlannerStats>* s);

protected:
    //member variables
    double finitial_eps, finitial_eps_planning_time, final_eps_planning_time, final_eps, dec_eps, final_epsilon;
    double repair_time;
    bool use_repair_time;

    std::vector<PlannerStats> stats;

    int num_of_expands_initial_solution;

    MDPConfig* MDPCfg_;

    bool bforwardsearch; //if true, then search proceeds forward, otherwise backward

    bool bsearchuntilfirstsolution; //if true, then search until first solution only (see planner.h for search modes)

    RRTSearchStateSpace_t* pSearchStateSpace_;

    unsigned int searchexpands;
    int MaxMemoryCounter;
    clock_t TimeStarted;
    FILE *fDeb;

    //member functions
    virtual void Initialize_searchinfo(CMDPSTATE* state, RRTSearchStateSpace_t* pSearchStateSpace);

    virtual CMDPSTATE* CreateState(int stateID, RRTSearchStateSpace_t* pSearchStateSpace);

    virtual CMDPSTATE* GetState(int stateID, RRTSearchStateSpace_t* pSearchStateSpace);

    virtual int ComputeHeuristic(CMDPSTATE* MDPstate, RRTSearchStateSpace_t* pSearchStateSpace);

    //initialization of a state
    virtual void InitializeSearchStateInfo(RRTState* state, RRTSearchStateSpace_t* pSearchStateSpace);

    //re-initialization of a state
    virtual void ReInitializeSearchStateInfo(RRTState* state, RRTSearchStateSpace_t* pSearchStateSpace);

    virtual void DeleteSearchStateData(RRTState* state);

    //used for backward search
    virtual void UpdatePreds(RRTState* state, RRTSearchStateSpace_t* pSearchStateSpace);

    //used for forward search
    virtual void UpdateSuccs(RRTState* state, RRTSearchStateSpace_t* pSearchStateSpace);

    virtual int GetGVal(int StateID, RRTSearchStateSpace_t* pSearchStateSpace);

    //returns 1 if the solution is found, 0 if the solution does not exist and 2 if it ran out of time
    virtual int ImprovePath(RRTSearchStateSpace_t* pSearchStateSpace, double MaxNumofSecs);

    virtual void BuildNewOPENList(RRTSearchStateSpace_t* pSearchStateSpace);

    virtual void Reevaluatefvals(RRTSearchStateSpace_t* pSearchStateSpace);
    virtual void Reevaluatehvals(RRTSearchStateSpace_t* pSearchStateSpace);

    //creates (allocates memory) search state space
    //does not initialize search statespace
    virtual int CreateSearchStateSpace(RRTSearchStateSpace_t* pSearchStateSpace);

    //deallocates memory used by SearchStateSpace
    virtual void DeleteSearchStateSpace(RRTSearchStateSpace_t* pSearchStateSpace);

    //debugging
    virtual void PrintSearchState(RRTState* state, FILE* fOut);

    //reset properly search state space
    //needs to be done before deleting states
    virtual int ResetSearchStateSpace(RRTSearchStateSpace_t* pSearchStateSpace);

    //initialization before each search
    virtual void ReInitializeSearchStateSpace(RRTSearchStateSpace_t* pSearchStateSpace);

    //very first initialization
    virtual int InitializeSearchStateSpace(RRTSearchStateSpace_t* pSearchStateSpace);

    virtual int SetSearchGoalState(int SearchGoalStateID, RRTSearchStateSpace_t* pSearchStateSpace);

    virtual int SetSearchStartState(int SearchStartStateID, RRTSearchStateSpace_t* pSearchStateSpace);

    //reconstruct path functions are only relevant for forward search
    virtual int ReconstructPath(RRTSearchStateSpace_t* pSearchStateSpace);

    virtual void PrintSearchPath(RRTSearchStateSpace_t* pSearchStateSpace, FILE* fOut);

    virtual int getHeurValue(RRTSearchStateSpace_t* pSearchStateSpace, int StateID);

    //get path
    virtual std::vector<int> GetSearchPath(RRTSearchStateSpace_t* pSearchStateSpace, int& solcost);

    virtual bool Search(RRTSearchStateSpace_t* pSearchStateSpace, std::vector<int>& pathIds, int & PathCost,
                        bool bFirstSolution, bool bOptimalSolution, double MaxNumofSecs);
						
private:
    // member variables
	errt_TreeType attr_treetype;
	errt_PathType attr_pathtype;
	double attr_stepsize;
	double attr_qgoalprob;
	int    attr_chaikiniter;
	double attr_maxplanningtime;

	// member functions
	std::string TreeTypeToStr(errt_TreeType treetype);
	std::string PathTypeToStr(errt_PathType pathtype);	

	void test_treelibrary();
	void get_envinformation();
	void ex_randnumbergenerator(int setofnumbers);
	void savepath(std::vector<RRTMapCoord *>* path_mapcoord_V, std::vector<int>* solution_stateIDs_V);
	int  compute_pathcost(std::vector<RRTMapCoord *>* path_mapcoord_V);
};

#endif
