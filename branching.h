
/* c++11 implementation of branching process simulator
 *
 * General scheme for simulating cell division processes
 * with arbitrary waiting time and state transitions 
 */

#ifndef BRANCHING_H
#define BRANCHING_H

#include <iostream>
#include <memory>
#include <vector>
#include <limits>
#include <queue>
#include <cmath>

/****************
 * CELL CLASSES *
 ****************/
/*
 * Simple Cell Class Interface
 * We only require 3 things to be implemented:
 * -- Constructor that initializes next_event_time
 * -- perform_next_event()
 * -- next_event_time
 */
template <class WorkingCell>
class Cell {
public:
  /* Cell constructor should somehow record current simulation time
   * AND should also properly intialize next_event_time!
   */

  /* Cell logic contained in perform_next_event
   * Perform pre-determined event for individual cell
   *  -- currently pure virtual that needs to be implemented
   */
  virtual std::vector< std::shared_ptr<WorkingCell> > perform_next_event() = 0;
  // Store time of next event to be stored -- should ALWAYS have a value
  double next_event_time;
};

/*
 * Test Cell for testing Cell interface
 */
class TestCell: public Cell<TestCell> {
public:
  TestCell(double t=0.0) : time(t) {get_next_event();}
  std::vector< std::shared_ptr<TestCell> > perform_next_event();

  // States are useful but not necessary
  // enums actually seem like more trouble than they're worth
  std::string get_state() {return state;}

private:  
  double time;
  // probabilistic decisions made in get_next_event
  void get_next_event() {next_event_time = time + 1;}
  // default starting state
  std::string state = "A";
};

/* 
 * TestCell Implementation - simply makes two cells upon division
 */
std::vector< std::shared_ptr<TestCell> > TestCell::perform_next_event() 
{
  time = next_event_time;
  std::vector< std::shared_ptr<TestCell> > new_cells;
  new_cells.emplace_back(std::make_shared<TestCell>(time));
  new_cells.emplace_back(std::make_shared<TestCell>(time));
  
  // get next_event for current cell
  get_next_event();
  return new_cells;
}


/********************
 * LISTENER CLASSES *
 ********************/
/*
 * Listener interface for iteratively updating simulation record data
 * We require 3 things to be implemented:
 * -- init(time,cells): initialize listener with time and starting cells
 * -- pop_event(time,cell): provide time and cell of a cell event
 * -- push_event(time,cells): provide time and result of the last pop_event
 */
template <class WorkingCell>
class Listener {
public:
  virtual void init(double time,std::vector< std::shared_ptr<WorkingCell> > &cells) = 0;
  virtual void pop_event(double time,std::shared_ptr<WorkingCell> c) = 0;
  virtual void push_event(double time,std::vector< std::shared_ptr<WorkingCell> > &new_cells) = 0;
};

/***************************
 * BRANCHING PROCESS CLASS *
 ***************************/
/* 
 * This Branching Process class holds the main loop logic 
 * for our simulation -- identifies times of cell events and calls them
 *
 * instantiation of BProcess requires a Cell class template for use
 * -- requires: next_event_time and perform_next_event()
 *
 * DefaultCell is simply the default Cell for easy process construction
 */
template <class WorkingCell,class DefaultCell = WorkingCell>
class BProcess {
public:
  // Simple initialization with provided number of default cells
  // Variadic template forwarding of constructor arguments for default cell class
  template <typename ...Args>
  BProcess(int num_init = 1,Args&&... params) 
  {
    // populate initial event heap with num_init cells
    for (int i = 0; i < num_init; ++i){
      // creates new default cells and puts their ptrs in the heap
      EHeap.emplace(std::make_shared<DefaultCell>(std::forward<Args>(params)...));
    }
  }

  // Fancier initialization with provided Cell vector
  BProcess(std::vector< std::shared_ptr<WorkingCell> >& initial_cells) 
  {
    for (auto ci : initial_cells){
      EHeap.push(ci);
    }
  }

  // main loop 
  void run(double TMAX = std::numeric_limits<double>::max(),
	   unsigned int NMAX = std::numeric_limits<unsigned int>::max());

  unsigned int num_cells(){return EHeap.size();}

  // only add listener, will initialize on run
  void add_listener(std::shared_ptr< Listener<WorkingCell> > lst)
  {
    LArray.push_back(lst);
  }

private:
  // cell event time compare for event ordering
  // minimum next_event_time based ordering
  struct CellComp
  {
    CellComp(){};
    bool operator() (std::shared_ptr<WorkingCell> a,std::shared_ptr<WorkingCell> b) const
    {
      return (a->next_event_time > b->next_event_time);
    }
  };

  // Central data structure for ordering events
  std::priority_queue< std::shared_ptr<WorkingCell>,
		       std::vector< std::shared_ptr<WorkingCell> >,
		       CellComp > EHeap;
  
  // vector for holding simulation listeners
  std::vector< std::shared_ptr< Listener<WorkingCell> > > LArray;
  void init_listeners(double time); 
};

// note these template member function initializations need
// to be in the header otherwise linking will fail
/*
 * Branching Process Implementation - initialize listeners
 */
template <class WorkingCell, class DefaultCell>
void BProcess<WorkingCell,DefaultCell>::init_listeners(double time)
{
  // kind of ugly solution but shouldn't be too much overhead for initialization
  auto hcopy = EHeap; // copy of EHeap
  std::vector< std::shared_ptr<WorkingCell> > init_cells;
  while (!hcopy.empty()){
    init_cells.push_back(hcopy.top());
    hcopy.pop();
  }
  
  // initialize each listener with this new vector
  for (auto l : LArray) {
    l->init(time,init_cells);
  }
}

/*
 * Branching Process Implementation - main simulation loop
 */
template <class WorkingCell, class DefaultCell>
void BProcess<WorkingCell,DefaultCell>::run(double TMAX,unsigned int NMAX)
{
  // current time in simulation 
  double current_time = 0.0;
  init_listeners(current_time); // initialize listeners
  while (current_time < TMAX && num_cells() < NMAX) {
    std::shared_ptr<WorkingCell> next_cell = EHeap.top(); // grab next cell event
    // update current time to next event time
    current_time = next_cell->next_event_time;
    for (auto l : LArray) {l->pop_event(current_time,next_cell);}
    EHeap.pop(); // remove that element

    std::vector< std::shared_ptr<WorkingCell> > new_cells = next_cell->perform_next_event();
    for (auto l: LArray) {l->push_event(current_time,new_cells);}
    for (auto new_cell : new_cells) {EHeap.push(new_cell);}
  }
}

#endif
