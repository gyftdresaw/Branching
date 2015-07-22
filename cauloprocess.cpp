
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
#include <algorithm>

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
class BasicCell : public Cell<BasicCell>, public std::enable_shared_from_this<BasicCell> {
public:
  BasicCell(std::function<double()> &w,std::function<int()> &p,double t=0.0) : 
    waiting(w),progeny(p),birth_time(t) {get_next_event();}
  std::vector< std::shared_ptr<BasicCell> > perform_next_event();
  double get_age(double t) {return t - birth_time;}
  std::string get_state() {return "";} // need this to work with state age listener
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
std::vector< std::shared_ptr<BasicCell> > BasicCell::perform_next_event()
{
  std::vector< std::shared_ptr<BasicCell> > new_cells;
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
 * Asymmetric Cell 
 * -- arbitrary waiting time distributions for division and swarmer->stalk transition
 * -- n progeny per division
 * -- two states - swarmer and stalk
 *
 * enable_shared_from_this inheritance is for retrieving a smart ptr to the object
 */
class AsymmetricCell : public Cell<AsymmetricCell>, public std::enable_shared_from_this<AsymmetricCell> {
public:
  // instantiated with three functions - two waiting, one progeny
  AsymmetricCell(std::function<double()> &w,std::function<double()> &tw,std::function<int()> &p,double t=0.0,std::string s="stalk") : 
    waiting(w),transition(tw),progeny(p),last_time(t),state(s) {get_next_event();}
  std::vector< std::shared_ptr<AsymmetricCell> > perform_next_event();
  double get_age(double t) {return t - last_time;}
  std::string get_state() {return state;} // need this to work with state age listener
private:
  double last_time;
  int nprogeny;
  std::function<double()> &waiting; // function REFERENCES -- crucial to maintain state
  std::function<double()> &transition;
  std::function<int()> &progeny;
  std::string state;
  void get_next_event();
};

// simply update time of next event and number of progeny to be produced
void AsymmetricCell::get_next_event() 
{
  if (state == "stalk") {
    next_event_time = last_time + waiting();
    nprogeny = progeny();
  }
  else if (state == "swarmer")
    next_event_time = last_time + transition();
  else
    std::cout << "state not found!" << std::endl;
}

/*
 * AsymmetricCell Implementation - produce n new cells given by progeny function
 */
std::vector< std::shared_ptr<AsymmetricCell> > AsymmetricCell::perform_next_event()
{
  std::vector< std::shared_ptr<AsymmetricCell> > new_cells;
  // keep current AsymmetricCell and produce more if needed
  if (state == "stalk") {
    // division based on nprogeny
    if (nprogeny >= 1){
      // update current cell and place in vector
      last_time = next_event_time;
      get_next_event();
      new_cells.push_back(shared_from_this());
      for (int i = 1; i < nprogeny; ++i){
	// generate new swarmer cells 
	new_cells.emplace_back(std::make_shared<AsymmetricCell>(waiting,transition,progeny,last_time,"swarmer"));
      }
    }
  }
  else if (state == "swarmer") {
    // transition to stalk state
    state = "stalk";
    last_time = next_event_time;
    get_next_event();
    new_cells.push_back(shared_from_this()); // add current cell
  }
  else
    std::cout << "state not found!" << std::endl;
  return new_cells;
}

/*
 * TimeKeeper is a helper class for keeping track of time event logic 
 * in various listener classes
 */
class TimeKeeper {
public:
  /* two methods of construction:
   * -- keeping time by event
   * -- keeping times given 
   */
  TimeKeeper(std::vector<double> &t,double PRC=1e-15) : prc(PRC) {times_set = false;}
  TimeKeeper(std::vector<double> &t,std::vector<double> &ts,double PRC=1e-15) : prc(PRC)
  {
    t = ts;
    times_set = true;
  }
  
   // for use in init, time is starting time
  void init_times(std::vector<double> &t,double time)
  {
    // if event recording, need to add starting time
    if (!times_set) {
      t.push_back(time);
    }
    
    // increment tindex until we reach first recording time
    while (t[tindex] < time && !AlmostEqual(t[tindex],time,prc)) {++tindex;}
  }

  // for new event, bool returned indicates if we need a new record item
  // add new time if needed as well
  bool new_entry(std::vector<double> &t,double time)
  {
    if (!times_set) {
      // update by event
      if (!AlmostEqual(time,t.back(),prc)){
	t.push_back(time); // add new time entry
	return true;
      }
    }
    return false;
  }

  // step tindex if record time is less than or equal to next event time
  // if true is returned, tindex has been incremented and record items should be updated
  // allowing tindex to go past time vector
  bool step_time(std::vector<double> &t,double time)
  {
    if (tindex < int(t.size())) {
      if (t[tindex] < time && !AlmostEqual(t[tindex],time,prc)) {
	++tindex;
	return true;
      }
    }
    return false;
  }

  // condition for valid tindex
  bool in_range(std::vector<double> &t) {return tindex < int(t.size());}
  
  // condition for recording items from new event
  bool record(std::vector<double> &t,double time)
  {
    if (tindex < int(t.size()))
      return (t[tindex] > time || AlmostEqual(t[tindex],time,prc));
    else
      return false;
  }
  int tindex = 0;
private:
  double prc;
  bool times_set;
  bool AlmostEqual(double a,double b,double EPS=1.0e-15){return std::abs(a - b) < EPS;}
};

/*
 * Record number of cells at each event time 
 */
template <class WorkingCell>
class NCellListener : public Listener<WorkingCell> {
public:
  // Constructors
  NCellListener(double PRC=1e-15) : tkeeper(times,PRC) {}
  NCellListener(std::vector<double> &ts,double PRC=1e-15) : tkeeper(times,ts,PRC) {}
  // intialize according to whether times are given or to be determined
  void init(double time,std::vector< std::shared_ptr<WorkingCell> > &cells)
  {
    tkeeper.init_times(times,time); // initialize times and tindex 
    N = std::vector<unsigned int>(times.size(),0); // initialize N records
    // current tindex is start
    N[tkeeper.tindex] += cells.size();
  }

  // will actually simply ignore pop_event and account for cell removed in push_event
  void pop_event(double time,std::shared_ptr<WorkingCell> c) {};
  void push_event(double time,std::vector< std::shared_ptr<WorkingCell> > &new_cells)
  {
    // for event based recording, tells us whether we need a new record
    if (tkeeper.new_entry(times,time)) {
      N.push_back(N.back()); // add new N entry
    }

    // step time and record state until we pass event
    while (tkeeper.step_time(times,time)) {
      if (tkeeper.in_range(times))
	N[tkeeper.tindex] = N[tkeeper.tindex - 1];
    }

    // if event is passed, record result
    if (tkeeper.record(times,time)) {
      N[tkeeper.tindex] += new_cells.size() - 1;
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
  TimeKeeper tkeeper;
};


/*
 * Record Age Distribution of cells at set times -- UNDER CONSTRUCTION
 */
/*
class AgeListener : public Listener {
public:
  AgeListener(double PRC=1e-15) {}
  AgeListener(std::vector<double> &ts,double PRC=1e-15) {}
  void init(double time,std::vector< std::shared_ptr<Cell> > &cells) {}
  void pop_event(double time,std::shared_ptr<Cell> c) {}
  void push_event(double time,std::vector< std::shared_ptr<Cell> > &new_cells) {}
private:
  std::vector<double> times;
  // need to store histograms at each time
  
  TimeKeeper tkeeper;
};
*/

/*
 * Record ages of all current cells at each given time point
 *
 * template class to guarantee that the cell has a get_age method
 * if states given, assumes cell has get_state method as well
 */
template <class WorkingCell>
class FullAgeListener : public Listener<WorkingCell> {
public:
  FullAgeListener(std::vector<double> &ts,std::vector<std::string> s,double PRC=1e-15) : tkeeper(times,ts,PRC),states(s) {}
  // if states not given, initialized with empty vector
  FullAgeListener(std::vector<double> &ts,double PRC=1e-15) : FullAgeListener(ts,std::vector<std::string>(),PRC) {}
  
  void init(double time,std::vector< std::shared_ptr<WorkingCell> > &cells)
  {
    tkeeper.init_times(times,time); // initialize times and tindex

    // initialize ages records -- empty vectors for each time
    std::size_t NStates = ((states.size() == 0) ? 1 : states.size());
    ages = std::vector< std::vector< std::vector<double> > >
      (times.size(),std::vector< std::vector<double> >(NStates,std::vector<double>()));

    // initialize current cells
    current_cells = cells;
  }

  void pop_event(double time,std::shared_ptr<WorkingCell> c)
  {
    // step time and record state until we pass event
    while (tkeeper.step_time(times,time)) {
      record_ages(times[tkeeper.tindex-1],ages[tkeeper.tindex-1]);
    }

    // need to remove cell from current cell list
    // should be guaranteed to be present
    auto cpos = std::find(current_cells.begin(),current_cells.end(),c);
    if (cpos == current_cells.end()) {
      std::cout << "Warning: couldn't find cell in full age listener" << std::endl;
    }
    current_cells.erase(cpos);
  }
  void push_event(double time,std::vector< std::shared_ptr<WorkingCell> >&new_cells)
  {
    // add new cells to current cell list
    for (auto c : new_cells) {
      current_cells.push_back(c);
    }
  }
  
  // output helpers
  void print()
  {
    for (auto t : times){std::cout << t << ' ';}
    std::cout << std::endl;
    for (auto s : states){std::cout << s << ' ';}
    std::cout << std::endl;
    for (auto sv : ages)
    {
      for (auto avector : sv)
      {
	for (auto a : avector) {std::cout << a << ' ';}
	std::cout << std::endl;
      }
    }
  }
  void write(std::string filename,bool include_header=true,bool append=false)
  {
    std::ofstream to_file;
    if (append) {
      to_file.open(filename,std::ios::out | std::ios::app);
    }
    else {
      to_file.open(filename,std::ios::out);
    }
    if (to_file.is_open()) {
      if (include_header){
	for (auto t : times) {to_file << t << '\t';}
	to_file << std::endl;
	if (states.size() != 0) {
	  for (auto s : states) {to_file << s << '\t';}
	  to_file << std::endl;
	}
      }
      for (auto sv : ages)
      {
	for (auto avector : sv)
	{
	  for (auto a : avector) {to_file << a << '\t';}
	  to_file << std::endl;
	}
      }
    }
  }

private:
  void record_ages(double time,std::vector< std::vector<double> > &dest)
  {
    // empty destination vector to be safe
    std::size_t NStates = ((states.size() == 0) ? 1 : states.size());
    dest = std::vector< std::vector<double> >(NStates,std::vector<double>());
    // we just need to add ages of all current cells
    for (auto c : current_cells) {
      // need to check state
      int dindex = ((states.size() == 0) ? 0 : state_index(c->get_state()));
      dest[dindex].push_back(c->get_age(time));
    }
  }
  int state_index(std::string s)
  {
    auto sit = std::find(states.begin(),states.end(),s);
    if (sit == states.end()) {
      std::cout << "State not found!" << std::endl;
      return -1;
    }
    else
      return std::distance(states.begin(),sit);
  }
  // hold age info at designated times
  // at each time point, set of age vectors for each state
  std::vector< std::vector< std::vector<double> > > ages;
  // continuously update cell list
  std::vector< std::shared_ptr<WorkingCell> > current_cells;
  // for keeping track of states
  std::vector<std::string> states;
  // for time keeping 
  std::vector<double> times;
  TimeKeeper tkeeper;
};

/*
 * Simply storing routine to run simulations that record N cells in time
 *
 * Needs some functions/rngs to be defined to work
 */
/*
void run_ncells() {
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
}
*/

/*
 * Storing routine to run basic cell simulation (normal symmetric division)
 *
 * Needs some functions/rngs to be defined to work
 */
/*
void run_basic() {
  // TESTING
  std::vector<double> times;
  double dmax = 14.0;
  for (double d = 0.0; d < dmax || std::abs(d-dmax) < 1e-6; d += 0.1) {
    times.push_back(d);
  }

  int ntrials = 200;
  std::string filename = "results/basic_fullage_gam5_02_p2.txt";
  for (int i = 0; i < ntrials; ++i) {
    auto FAlst = std::make_shared< FullAgeListener<BasicCell> >(times,1e-6);
    BProcess<BasicCell> bp(1,gam_wt,default_progeny);
    bp.add_listener(FAlst);
    bp.run(dmax,1e8);
    if (i == 0) {
      FAlst->write(filename); // first run start new file
    }
    else {
      FAlst->write(filename,false,true); // append runs after first
    }
    std::cout << i << std::endl;
  }  
}
*/

#include <iomanip>
#include <map>

// basic testing
int main(int argc, char const ** argv){
  // for random numbers
  std::random_device rd;
  std::mt19937_64 gen(rd());

  // construct waiting time function 
  std::exponential_distribution<double> exp(1.0);
  std::function<double()> exp_wt = [&](){return exp(gen);};
  std::gamma_distribution<double> gam(3.0,1.0/3.0);
  std::function<double()> gam_wt = [&](){return gam(gen);};
  std::function<double()> default_waiting = [](){return 1.0;};
  // second gamma for swarmers
  std::gamma_distribution<double> gam2(5.0,0.2);
  std::function<double()> gam2_wt = [&](){return gam2(gen);};
  // second exp for swarmers
  std::exponential_distribution<double> exp2(0.5);
  std::function<double()> exp2_wt = [&](){return exp2(gen);};

  // bimodal waiting time
  std::bernoulli_distribution bern(0.5);
  std::function<double()> bimod_wt = [&](){return bern(gen) ? exp(gen) : gam(gen);};
  
  // construct simple number of progeny produced per division
  std::function<int()> default_progeny = [](){return 2;};

  // TESTING
  std::vector<double> times;
  double dmax = 10.0; // 14.0
  for (double d = 0.0; d < dmax || std::abs(d-dmax) < 1e-6; d += 0.1) {
    times.push_back(d);
  }

  // ASYMMETRIC CASE
  /*
  int ntrials = 2000;
  std::string filename = "results/asymmetric_fullage_exp1_exp2_p2_2000trajectories.txt";
  for (int i = 0; i < ntrials; ++i) {
    auto FAlst = std::make_shared< FullAgeListener<AsymmetricCell> >(times,std::vector<std::string>{"stalk","swarmer"},1e-6);
    BProcess<AsymmetricCell> bp(1,exp_wt,exp2_wt,default_progeny,0.0,"stalk");
    bp.add_listener(FAlst);
    bp.run(dmax,1e8);
    if (i == 0) {
      FAlst->write(filename); // first run start new file
    }
    else {
      FAlst->write(filename,false,true); // append runs after first
    }
    std::cout << i << std::endl;
  }
  */
  
  // BASIC CASE
  int ntrials = 2000;
  std::string filename = "results/basic_fullage_gam3_03_p2_2000trajectories_t10.txt";
  for (int i = 0; i < ntrials; ++i) {
    auto FAlst = std::make_shared< FullAgeListener<BasicCell> >(times,1e-6);
    BProcess<BasicCell> bp(1,gam_wt,default_progeny);
    bp.add_listener(FAlst);
    bp.run(dmax,1e8);
    if (i == 0) {
      FAlst->write(filename); // first run start new file
    }
    else {
      FAlst->write(filename,false,true); // append runs after first
    }
    std::cout << i << std::endl;
  }  
  
  return 0;
}
