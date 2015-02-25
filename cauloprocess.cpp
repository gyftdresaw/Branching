
/*
 * Some Caulobacter Project specific Cell/Listener Classes
 * for branching process simulation
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <functional>
#include <random>

// for use with branching header library
#include "branching.h"

/*
 * Basic Cell 
 * -- arbitrary waiting time distribution
 * -- n progeny per division
 * -- no cell states 
 *
 * enable_shared_from_this inheritance is for retrieving a smart ptr to the object
 */

class BasicCell : public Cell, public std::enable_shared_from_this<BasicCell> {
public:
  BasicCell(std::function<double()> &w,std::function<int()> &p,double t=0.0) : 
    waiting(w),progeny(p),birth_time(t) {get_next_event();}
  std::vector< std::shared_ptr<Cell> > perform_next_event();
private:
  double birth_time;
  int nprogeny;
  std::function<double()> &waiting; // function REFERENCES -- crucial to maintain state
  std::function<int()> &progeny;
  void get_next_event();
};

// simply update time of next event and number of progeny to be produced
void BasicCell::get_next_event() 
{
  next_event_time = birth_time + waiting();
  nprogeny = progeny();
}

/*
 * BasicCell Implementation - produce n new cells given by progeny function
 */
std::vector< std::shared_ptr<Cell> > BasicCell::perform_next_event()
{
  std::vector< std::shared_ptr<Cell> > new_cells;
  // keep current BasicCell and produce more
  if (nprogeny >= 1){
    // update current cell and place in vector
    birth_time = next_event_time;
    get_next_event();
    new_cells.push_back(shared_from_this());
    for (int i = 1; i < nprogeny; ++i){
      new_cells.emplace_back(std::make_shared<BasicCell>(waiting,progeny,birth_time));
    }
  }
  
  return new_cells;
}

/*
 * Record number of cells at each event time 
 *
 * [TODO I think the time/event based inputs could be refactored 
 * or even inherited more coherently]
 */
class NCellListener : public Listener {
public:
  // Constructors
  NCellListener(double PRC=1e-15) : prc(PRC) {times_set = false;}
  NCellListener(std::vector<double> &ts,double PRC=1e-15) : times(ts), prc(PRC) {times_set = true;}
  // intialize according to whether times are given or to be determined
  void init(double time,std::vector< std::shared_ptr<Cell> > &cells)
  {
    // record times aren't given a priori 
    if (!times_set) {
      times.push_back(time); // simply set starting time as the only record time
    }
    N = std::vector<unsigned int>(times.size(),0);
    while (times[tindex] < time && !AlmostEqual(times[tindex],time,prc)) {++tindex;}
    // current index is start
    N[tindex] += cells.size();
  }

  // will actually simply ignore pop_event and account for cell removed in push_event
  void pop_event(double time,std::shared_ptr<Cell> c) {};
  void push_event(double time,std::vector< std::shared_ptr<Cell> > &new_cells)
  {
    if (!times_set) {
      // update by event
      if (!AlmostEqual(time,times.back(),prc)){
	times.push_back(time); // add new time entry
	N.push_back(N.back()); // add new N entry
      }
    }
    // update according to time
    while (times[tindex] < time && !AlmostEqual(times[tindex],time,prc) &&
	   tindex < int(times.size()) - 1) {
      ++tindex;
      N[tindex] = N[tindex - 1];
    }
    if (times[tindex] > time || AlmostEqual(times[tindex],time,prc)) {
      N[tindex] += new_cells.size() - 1;
    }
  }

  // output helpers
  void print()
  {
    for (auto t : times){std::cout << t << ' ';}
    std::cout << std::endl;
    for (auto n : N){std::cout << n << ' ';}
    std::cout << std::endl;
  }
  void write(std::string filename,bool include_times=true,bool append=false)
  {
    std::ofstream to_file;
    if (append) {
      to_file.open(filename,std::ios::out | std::ios::app);
    }
    else {
      to_file.open(filename,std::ios::out);
    }
    if (to_file.is_open()) {
      if (include_times){
	for (auto t : times) {to_file << t << '\t';}
	to_file << std::endl;
      }
      for (auto n : N) {to_file << n << '\t';}
      to_file << std::endl;
    }
  }
private:
  std::vector<double> times;
  std::vector<unsigned int> N;
  // for event based recording
  double prc; // time difference precision
  bool AlmostEqual(double a,double b,double EPS=1.0e-15){return std::abs(a - b) < EPS;}
  // for set time based recording
  bool times_set;
  int tindex = 0; // tindex is always pointing to the next record time
};

/*
 * Record Age Distribution of cells at set times
 */
class AgeListener : public Listener {
public:

private:
  std::vector<double> times;
  // need to store histograms at each time
  
  // for set time based recording
  bool times_set;
  int tindex = 0;
};

// basic testing
int main(int argc, char const ** argv){
  // for random numbers
  std::random_device rd;
  std::mt19937_64 gen(rd());

  // construct waiting time function 
  std::exponential_distribution<double> exp(1.0);
  std::function<double()> exp_wt = [&](){return exp(gen);};
  std::gamma_distribution<double> gam(6.0,0.2);
  std::function<double()> gam_wt = [&](){return gam(gen);};
  std::function<double()> default_waiting = [](){return 1.0;};

  // construct simple number of progeny produced per division
  std::function<int()> default_progeny = [](){return 3;};

  // TESTING
  auto Nlst = std::make_shared<NCellListener>(1e-6);
  BProcess<BasicCell> bp(1,default_waiting,default_progeny);
  bp.add_listener(Nlst);
  bp.run(4.0,1000);
  Nlst->print();
  
  /*
  // construct vector for times to be recorded
  std::vector<double> times;
  double dmax = 14.0;
  for (double d = 0.0; d < dmax || std::abs(d-dmax) < 1e-6; d += 0.01) {
    times.push_back(d);
  }

  int ntrials = 200;
  std::string filename = "results/basic_ncell_gam6_02_p3.txt";
  for (int i = 0; i < ntrials; ++i) {
    auto Nlst = std::make_shared<NCellListener>(times,1e-6);
    BProcess<BasicCell> bp(1,gam_wt,default_progeny);
    bp.add_listener(Nlst);
  
    bp.run(dmax,1e8);
    // Nlst->print();
    if (i == 0) {
      Nlst->write(filename); // first run start new file
    }
    else {
      Nlst->write(filename,false,true); // append runs after first
    }
    std::cout << i << std::endl;
  }
  */
  return 0;
}
