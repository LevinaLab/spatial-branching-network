// --------------------------------------------------------- //
//                  IBNM                                     //
//   An implementation of a branching-like process in 2D     //
//   Original implemtantion by Roxana Zeraati, translation   //
//   to C++, Victor Buendia                                  //
//   Last revision May 2022                                  //
// --------------------------------------------------------- //

#include<iostream>
#include<cstdlib>
#include<vector>
#include<random>
#include<cmath>
#include<fstream>
#include<algorithm>

using namespace std;

// --- Preprocessor magic ---

//Different uses of the program
#define AVAL 0
#define DIAGRAM 1 
#define AVAL_SHAPE 2
#define TSERIES 3
#define SUBSAMPLING 4

//Default mode
#ifndef MODE
#define MODE SUBSAMPLING 
#endif

//4 or 8 neighbours
#define MOORE 0
#define MANHATTAN 1

#ifndef LATTICE
#define LATTICE MOORE
#endif

//Should we rewire the network?
#define FALSE 0
#define TRUE 1

#ifndef MF
#define MF FALSE 
#endif

// --- Global constants ---

int L = 64;
int N = L*L;    // System size N=L*L
double sub_ratio = 0.05;

int k;          //Neighbourhood size
int num_neighs;

double p_s;   //Self-inducing probability
double p_n;   //Interaction with neighbours

double p_rewire; //Probability to rewire

const double external_field = 1e-6; //External field  

//Standard simulation times (for non-avalanches)
const int tf = 1000;//100000;
const int t_relax = 0;
const int t_meas = 10;//100;
const int nruns = 100;//1000;
//Commented: values for paper. current, for mf check

string filename;    //Where to store results

//RNG
mt19937 gen(151234553);    
uniform_real_distribution<double> ran_u(0.0, 1.0);  

// --- Function declarations ---

int compute_numneighs();
int manh_distance(const int x1, const int y1, const int x2, const int y2);
void create_neighbour_table(vector<vector<int>> &neighs);
void get_indices_subsampling(const bool random_idxs, const int Lsub, vector<int> &subindex);
int update_lattice(vector<bool> &lattice, const vector<vector<int>> &neighs);
int update_lattice(vector<bool> &lattice, const vector<vector<int>> &neighs, const double external_field);
int update_lattice_sub(vector<bool> &lattice, const vector<vector<int>> &neighs, const vector<int> &subindex, int &sub_activity);
int update_lattice_mf(vector<bool> &lattice, const int n_active); 

void make_avalanches(vector<bool> &lattice, const vector<vector<int>> &neighs, const int naval, string avalpath);
void make_avalanches_subasample(const int Lsub, vector<bool> &lattice, const vector<vector<int>> &neighs, const vector<int> &subindices, const int naval, string avalpath);
void make_avalanches_mf(vector<bool> &lattice, const vector<vector<int>> &neighs, const int naval, string avalpath);

void time_series(vector<bool> &lattice, const vector<vector<int>> &neighs, const string outpath);
void make_diagram(vector<bool> &lattice, const vector<vector<int>> &neighs, const double q0, const double qf, const double nq, const string pdpath);
void make_diagram_mf(vector<bool> &lattice, const vector<vector<int>> &neighs, const double q0, const double qf, const double nq, const string pdpath);
void avalanche_profile(vector<bool> &lattice, const vector<vector<int>> &neighs, const int naval, const string avalpath);

// --- Main call ---

int main(int argc, char *argv[])
{
    //State of the system and neighbour list
    vector<bool> lattice;
    vector< vector<int>> neighs;

    int id;
    int naval;
    //Computation of avalanches
    #if MODE==AVAL || MODE==AVAL_SHAPE
        //Get all arguments from above
        if (argc == 9)
        {
            L   = stoi(argv[1]);
            p_s = stod(argv[2]);
            p_n = stod(argv[3]);
            k   = stoi(argv[4]);
            naval = stoi(argv[5]);
            p_rewire = stod(argv[6]);
            filename = string(argv[7]);
            id = stoi(argv[8]);
        }
        else 
        {
            cout << "Wrong number of parameters: call as" << endl << "./exe L p_s p_n k naval filename index" << endl;
            return EXIT_SUCCESS;
        }

        //Random seed depending on thread
        gen.seed(56123455 + id*id*15312);

        //Define lattice
        N = L*L;
        num_neighs = compute_numneighs();
        create_neighbour_table(neighs);

        //Core computation
        #if MODE==AVAL
            #if MF==TRUE
                make_avalanches_mf(lattice, neighs, naval, filename);
            #else
                make_avalanches(lattice, neighs, naval, filename);
            #endif
        #elif MODE==AVAL_SHAPE
            avalanche_profile(lattice, neighs, naval, filename);
        #endif

    #elif MODE == DIAGRAM
        double q0,qf;
        int nq;
        //Get all arguments from above
        if (argc == 9)
        {
            L   = stoi(argv[1]);
            p_s = stod(argv[2]);
            k   = stoi(argv[3]);
            q0 = stod(argv[4]);
            qf = stod(argv[5]);
            nq = stoi(argv[6]);
            p_rewire = stod(argv[7]);
            filename = string(argv[8]);
        }
        else 
        {
            cout << "Wrong number of parameters: call as" << endl << "./exe L p_s p_n k naval filename" << endl;
            return EXIT_SUCCESS;
        }

        //Define lattice
        N = L*L;
        num_neighs = compute_numneighs();
        create_neighbour_table(neighs);
        lattice = vector<bool>(N);

        //Core computation
        #if MF==TRUE 
            make_diagram_mf(lattice, neighs, q0, qf, nq, filename); 
        #else
            make_diagram(lattice, neighs, q0, qf, nq, filename); 
        #endif

    #elif MODE == TSERIES 
        //Get all arguments from above
        if (argc == 7)
        {
            L   = stoi(argv[1]);
            p_s = stod(argv[2]);
            p_n = stod(argv[3]);
            k   = stoi(argv[4]);
            p_rewire = stod(argv[5]);
            filename = string(argv[6]);
        }
        else 
        {
            cout << "Wrong number of parameters: call as" << endl << "./exe L p_s p_n k naval filename" << endl;
            return EXIT_SUCCESS;
        }

        //Define lattice
        N = L*L;
        num_neighs = compute_numneighs();
        create_neighbour_table(neighs);
        lattice = vector<bool>(N);

        time_series(lattice, neighs, filename);

        //Core computation

    #elif MODE==SUBSAMPLING

        vector<int> subindices;
        int Lsub;

        //Get all arguments from above
        if (argc == 10)
        {
            L   = stoi(argv[1]);
            p_s = stod(argv[2]);
            p_n = stod(argv[3]);
            k   = stoi(argv[4]);
            sub_ratio = stod(argv[5]);
            naval = stoi(argv[6]);
            p_rewire = stod(argv[7]);
            filename = string(argv[8]);
            id = stoi(argv[9]);
        }
        else 
        {
            cout << "Wrong number of parameters: call as" << endl << "./exe L p_s p_n k naval filename index" << endl;
            return EXIT_SUCCESS;
        }

        //Random seed depending on thread
        gen.seed(56123455 + id*id*15312);

        //Define lattice
        N = L*L;
        num_neighs = compute_numneighs();
        create_neighbour_table(neighs);

        //Initialize the subsampled network
        Lsub = int(sqrt(N*sub_ratio));
        get_indices_subsampling(true, Lsub, subindices);

        //Core computation
        make_avalanches_subasample(Lsub, lattice, neighs, subindices, naval, filename);

    #endif

    return EXIT_SUCCESS;
}

int compute_numneighs()
{
    #if LATTICE==MOORE
    return 4*k*(k+1);
    #else
    return 2*k*(k+1);
    #endif
}

int manh_distance(const int x1, const int y1, const int x2, const int y2)
{
    int dx = min(abs(x2 - x1) , L - abs(x2 - x1));
    int dy = min(abs(y2 - y1) , L - abs(y2 - y1));
    return dx + dy;
}

void create_neighbour_table(vector<vector<int>> &neighs)
{
    int i,j,x,y;    //Count indices
    int index;      //Selected node
    int nindex;     //Selected neighbour
    int h,v;        //To select horizontal and vertical with respect to node

    int random_neigh; //Acount for a random index
    bool already_connected;
    int curr_neighs;

    neighs = vector<vector<int>>(N, vector<int>(num_neighs)); 

    //For each node in the lattice...
    for (x=0; x < L; x++)
    {
        for (y=0; y < L; y++)
        {
            //...get its index..
            index = x + y*L;

            //...and explore its neighbourhood around...
            nindex = 0;
            for (i = -k; i <= k; i++)
            {

                //Periodic horizontal boundaries
                h = x+i;
                if (h < 0) h += L;
                else if (h >= L) h -= L;

                for (j= -k; j <= k; j++)
                {
                    //Skip (0,0), which is yourself, as neighbour
                    if (i == 0 && j ==0) continue;

                    //Periodic vertical boundaries
                    v = y+j;
                    if (v < 0) v += L;
                    else if (v >= L) v -= L;

                    if (ran_u(gen) > p_rewire) //Not rewiring
                    {
                        //Neighbour index
                        #if LATTICE==MOORE
                        neighs[index][nindex] = h + v*L;
                        nindex++;
                        #else
                        if (manh_distance(x,y, h,v) <= k)
                        {
                            neighs[index][nindex] = h + v*L;
                            nindex++;
                        }
                        #endif
                    }
                    else //Rewiring
                    {
                        do
                        {
                            random_neigh = int(ran_u(gen) * N);
                            curr_neighs = 0;
                            already_connected = false;
                            while (curr_neighs < nindex && !already_connected)
                            {
                                already_connected = random_neigh == neighs[index][curr_neighs] || random_neigh == index;
                                curr_neighs++;
                            }
                            
                        } while (already_connected);
                        
                        neighs[index][nindex] = random_neigh;
                        nindex++;
                    }
                }
            } //end neighbourhood 
        }
    } //end lattice sweep
    return;
}

//Indices of the nodes we are gona subsample
void get_indices_subsampling(const bool random_idxs, const int Lsub, vector<int> &subindex)
{
    int x, y;
    int n_sub = Lsub*Lsub;

    if (random_idxs)
    {
        //First generate all numbers from 0 to N-1, then shuffle then, finally cut at n_sub
        //This gives non-repeated random indices
        vector<int> allnumbers(N);
        iota(allnumbers.begin(), allnumbers.end(), 0); 
        shuffle(allnumbers.begin(), allnumbers.end(), gen);
        subindex = vector<int>(allnumbers.begin(), allnumbers.begin()+n_sub);
        //Sort the vector because then it is easier in the update function to check whether if a node is in the subsample indices
        sort(subindex.begin(), subindex.end());
    }
    else
    {
        //Grab subsampled thing in a small square 
        subindex = vector<int>(n_sub);
        int k = 0;

        for (y=0; y < Lsub; y++)
        {
            for (x=0; x < Lsub; x++)
            {
                subindex[k] = x + Lsub*y; 
                k++;
            }
        }
    }
    return;
}

int update_lattice(vector<bool> &lattice, const vector<vector<int>> &neighs)
{

    //Counters
    int i,j;

    //To check the active neighbours
    bool neigh_state;

    //My probability to avoid decay
    double p_not_active;

    //Lattice total activity
    int activity;

    //Copy state to avoid any overwriting
    //vector<bool> old_lattice(lattice);
    vector<bool> old_lattice(N);
    old_lattice.swap(lattice);

    //Update each site
    activity = 0;
    for (i=0; i < N; i++)
    {
        //If I am active, I have some p of decaying.
        //If not, in principle I will be inactive for sure 
        p_not_active = old_lattice[i] ? 1.0 - p_s : 1.0;

        //Get some probability to activate using active neighbours
        for (j=0; j < num_neighs; j++)          
        {
            neigh_state = old_lattice[neighs[i][j]];
            if (neigh_state) p_not_active *= 1-p_n; 
        }

        //Update node, keep track of current activity
        lattice[i] = ran_u(gen) < 1.0 - p_not_active;
        activity += lattice[i];
    }

    return activity;
}


int update_lattice(vector<bool> &lattice, const vector<vector<int>> &neighs, const double external_field)
{

    //Counters
    int i,j;

    //To check the active neighbours
    bool neigh_state;

    //My probability to avoid decay
    double p_not_active;

    //Lattice total activity
    int activity;

    //Copy state to avoid any overwriting
    vector<bool> old_lattice(lattice);

    //Update each site
    activity = 0;
    for (i=0; i < N; i++)
    {
        //If I am active, I have some p of decaying.
        //If not, in principle I will be inactive for sure 
        p_not_active = old_lattice[i] ? 1.0 - p_s : 1.0;
        p_not_active *= 1.0 - external_field;

        //Get some probability to activate using active neighbours
        for (j=0; j < num_neighs; j++)          
        {
            neigh_state = old_lattice[neighs[i][j]];
            if (neigh_state) p_not_active *= 1-p_n; 
        }

        //Update node, keep track of current activity
        lattice[i] = ran_u(gen) < 1.0 - p_not_active;
        activity += lattice[i];
    }

    return activity;
}

int update_lattice_sub(vector<bool> &lattice, const vector<vector<int>> &neighs, const vector<int> &subindex, int &sub_activity)
{

    //Counters
    int i,j,k;

    //To check the active neighbours
    bool neigh_state;

    //My probability to avoid decay
    double p_not_active;

    //Lattice total activity
    int activity;

    //Copy state to avoid any overwriting
    vector<bool> old_lattice(lattice);

    //Update each site
    activity = 0;
    sub_activity = 0;
    k = 0;
    for (i=0; i < N; i++)
    {
        //If I am active, I have some p of decaying.
        //If not, in principle I will be inactive for sure 
        p_not_active = old_lattice[i] ? 1.0 - p_s : 1.0;

        //Get some probability to activate using active neighbours
        for (j=0; j < num_neighs; j++)          
        {
            neigh_state = old_lattice[neighs[i][j]];
            if (neigh_state) p_not_active *= 1-p_n; 
        }

        //Update node, keep track of current activity
        lattice[i] = ran_u(gen) < 1.0 - p_not_active;
        activity += lattice[i];

        //Subsampled activity. Use the subindex vector is ordered
        if (k < subindex.size() && i == subindex[k])
        {
            sub_activity += lattice[i];
            k++;
        }
    }

    return activity;
}

int update_lattice_mf(vector<bool> &lattice, const int n_activos)
{

    //Counters
    int i;

    //My probability to avoid decay
    double p_not_active;
    const double act_density = n_activos / (1.0*N);
    double prob_density = pow(1.0 - p_n, act_density);

    //Lattice total activity
    int activity;

    //Copy state to avoid any overwriting
    //vector<bool> old_lattice(lattice);
    vector<bool> old_lattice(N);
    old_lattice.swap(lattice);

    //Update each site
    activity = 0;
    for (i=0; i < N; i++)
    {
        //If I am active, I have some p of decaying.
        //If not, in principle I will be inactive for sure 
        p_not_active = old_lattice[i] ? 1.0 - p_s : 1.0;
        p_not_active *= prob_density;

        //Update node, keep track of current activity
        lattice[i] = ran_u(gen) < 1.0 - p_not_active;
        activity += lattice[i];
    }

    return activity;
}


void make_avalanches(vector<bool> &lattice, const vector<vector<int>> &neighs, const int naval, const string avalpath)
{
    int i;                          //Counter
    int center = L/2 + L/2 * L;     //Center of the lattice

    //Store size and duration of the avalanche
    int rho;
    int size, t;

    //Where to store results
    ofstream output;

    output.open(avalpath);
    for (i=0; i < naval; i++)
    {
        //Initial conditions, single seed
        lattice = vector<bool>(N, 0);
        lattice[center] = 1;

        //Initialize counters
        t = 0;
        size = 1;
        rho = 1;

        //Run until absorbing state or make sure we are in supercritical state 
        while (rho > 0 && size < 1e10)
        {
            rho = update_lattice(lattice, neighs);
            t++;

            size += rho;
        }        

        //Store avalanche to file
        output << size << " " << t << endl;

    }
    output.close();

    return;
}

void make_avalanches_subasample(const int Lsub, vector<bool> &lattice, const vector<vector<int>> &neighs, const vector<int> &subindices, const int naval, string avalpath)
{
    int i;                          //Counter
    int xstart, ystart; 

    //Store size and duration of the avalanche
    int rho;
    int sub_activity = 0;
    int size, t;
    int total_size, total_t;

    //Where to store results
    ofstream output;

    output.open(avalpath);
    for (i=0; i < naval; i++)
    {
        //Initial conditions, single seed
        lattice = vector<bool>(N, 0);
        xstart = Lsub * ran_u(gen);
        ystart = Lsub * ran_u(gen);
        lattice[xstart + L*ystart] = 1;

        //Initialize counters
        t = 0;
        size = 1;
        rho = 1;

        total_size = 1;
        total_t = 0;

        //Run until absorbing state or make sure we are in supercritical state 
        while (rho > 0 && total_size < 1e10)
        {
            rho = update_lattice_sub(lattice, neighs, subindices, sub_activity);
            
            //Subsampled
            size += sub_activity;
            t += sub_activity > 0;

            //Total
            total_size += rho;
            total_t++;

            //Check end of subsampled avalanche    
            if (sub_activity == 0 && size > 0)
            {
                //Store avalanche to file
                output << size << " " << t << endl;
                size = 0;
                t = 0;
            }
        }        

        //Output real size and time. Make negative to differentiate 
        //with the subsampled ones! 
        output << -total_size << " " << -total_t << endl;
    }
    output.close();

    return;
}

void make_avalanches_mf(vector<bool> &lattice, const vector<vector<int>> &neighs, const int naval, const string avalpath)
{
    int i;                          //Counter
    int center = L/2 + L/2 * L;     //Center of the lattice

    //Store size and duration of the avalanche
    int rho;
    int size, t;

    //Where to store results
    ofstream output;

    output.open(avalpath);
    for (i=0; i < naval; i++)
    {
        //Initial conditions, single seed
        lattice = vector<bool>(N, 0);
        lattice[center] = 1;

        //Initialize counters
        t = 0;
        size = 1;
        rho = 1;

        //Run until absorbing state
        while (rho > 0 && size < 1e10)
        {
            rho = update_lattice_mf(lattice, rho);
            t++;

            size += rho;
        }        

        //Store avalanche to file
        output << size << " " << t << endl;

    }
    output.close();

    return;
}

void time_series(vector<bool> &lattice, const vector<vector<int>> &neighs, const string outpath)
{
    int i;  //Counter

    //Number of active individuals and time
    int nactive;
    int t;

    //Where to store results
    ofstream output;

    //Initial conditions
    lattice = vector<bool>(N, 0);

    #if MF==TRUE
        //Initial density 15%
        for (i=0; i < N; i++) 
        {
            lattice[i] = ran_u(gen) < 0.15;
            nactive += lattice[i];
        }
        for (t=0; t < t_relax; t++) nactive = update_lattice_mf(lattice, nactive);
    #else
        for (i=0; i < N; i++) lattice[i] = ran_u(gen) < 0.15;
        for (t=0; t < t_relax; t++) nactive = update_lattice(lattice, neighs, external_field);
    #endif

    nactive = 0.15 * N;
    //Run for some time
    output.open(outpath);
    #if MF==TRUE
        for (t=0; t < tf; t++)
        {
            output << t << " " << nactive / (N * 1.0) << endl;
            nactive = update_lattice_mf(lattice, nactive);
        }        
    #else
        for (t=0; t < tf; t++)
        {
            output << t << " " << nactive / (N * 1.0) << endl;
            nactive = update_lattice(lattice, neighs, external_field);
        }        
    #endif
    output.close();

    return;
}

void make_diagram(vector<bool> &lattice, const vector<vector<int>> &neighs, const double q0, const double qf, const double nq, const string pdpath)
{
    //Get discretization in branching parameter
    double dq = (qf - q0) / (1.0 * nq);
    double q;

    //To get order parameters
    double nactive;
    double rho, rho2;
    double avrho, avrho2;
    int survival_runs;
    int number_measures;

    //Output and counters
    int t;
    ofstream output;
    int i,j,run;

    //Start simulations
    output.open(pdpath);
    for (q=q0; q <= qf; q += dq)
    {
        //Compute the infection prob using the branching ratio
        p_n = (q - p_s) / num_neighs;

        //Then let's do some run, look for the survival ones
        survival_runs = 0;
        avrho = avrho2 = 0.0;
        for (run=0; run < nruns; run++)
        //while (survival_runs < nruns)
        {
            nactive = 0;
            number_measures = 0;

            //Initial density 15%
            for (i=0; i < N; i++) 
            {
                lattice[i] = ran_u(gen) < 0.15;
                nactive += lattice[i];
            }

            //Thermalise
            t = 0;
            while (nactive > 0 && t < t_relax) 
            {
                nactive = update_lattice(lattice, neighs);
                t++;
            }

            //If after thermalization there is still activity, measure 
            if (nactive > 0)
            {
                //Set up stuff
                t = 0;
                rho = rho2 = 0.0;
                
                //Update and measure
                while (nactive > 0 && t < tf) 
                {
                    nactive = update_lattice(lattice, neighs);

                    //Make some measurements from time to time
                    if (t % t_meas == 0)
                    {
                        rho += nactive;
                        rho2 += nactive * nactive;
                        number_measures++;
                    }

                    t++;
                }

                //Finish averages and write IF the run was sucessful
                if (nactive > 0)
                {
                    avrho += rho / (1.0 * number_measures * N);
                    avrho2 += rho2 / (1.0 * number_measures * N*N);
                    survival_runs++;
                }
            }
        } //End of run loop

        if (survival_runs > 0)
        {
            avrho /= survival_runs;
            avrho2 /= survival_runs;
            output << q << " " << avrho << " " << avrho2 - avrho*avrho << " " << survival_runs << endl; 
        }
        else output << q << " " << 0.0 << " " << 0.0 << " " << 0 << endl;
        
        
    }
    output.close();

}

void make_diagram_mf(vector<bool> &lattice, const vector<vector<int>> &neighs, const double q0, const double qf, const double nq, const string pdpath)
{
    //Get discretization in branching parameter
    double dq = (qf - q0) / (1.0 * nq);
    double q;

    //To get order parameters
    double nactive;
    double rho, rho2;
    double avrho, avrho2;
    int survival_runs;
    int number_measures;

    //Output and counters
    int t;
    ofstream output;
    int i,j,run;

    //Start simulations
    output.open(pdpath);
    for (q=q0; q <= qf; q += dq)
    {
        //Compute the infection prob using the branching ratio
        p_n = q;

        //Then let's do some run, look for the survival ones
        survival_runs = 0;
        avrho = avrho2 = 0.0;
        for (run=0; run < nruns; run++)
        //while (survival_runs < nruns)
        {
            nactive = 0;
            number_measures = 0;

            //Initial density 15%
            for (i=0; i < N; i++) 
            {
                lattice[i] = ran_u(gen) < 0.15;
                nactive += lattice[i];
            }

            //Thermalise
            t = 0;
            while (nactive > 0 && t < t_relax) 
            {
                nactive = update_lattice_mf(lattice, nactive);
                t++;
            }

            //If after thermalization there is still activity, measure 
            if (nactive > 0)
            {
                //Set up stuff
                t = 0;
                rho = rho2 = 0.0;
                
                //Update and measure
                while (nactive > 0 && t < tf) 
                {
                    nactive = update_lattice_mf(lattice, nactive);

                    //Make some measurements from time to time
                    if (t % t_meas == 0)
                    {
                        rho += nactive;
                        rho2 += nactive * nactive;
                        number_measures++;
                    }

                    t++;
                }

                //Finish averages and write IF the run was sucessful
                if (nactive > 0)
                {
                    avrho += rho / (1.0 * number_measures * N);
                    avrho2 += rho2 / (1.0 * number_measures * N*N);
                    survival_runs++;
                }
            }
        } //End of run loop

        if (survival_runs > 0)
        {
            avrho /= survival_runs;
            avrho2 /= survival_runs;
            output << q << " " << avrho << " " << avrho2 - avrho*avrho << " " << survival_runs << endl; 
        }
        else output << q << " " << 0.0 << " " << 0.0 << " " << 0 << endl;
        
    }
    output.close();
}

void avalanche_profile(vector<bool> &lattice, const vector<vector<int>> &neighs, const int naval, const string avalpath)
{
    int i;                          //Counter
    int center = L/2 + L/2 * L;     //Center of the lattice

    //Store size and duration of the avalanche
    int rho;
    int size, t;

    //Where to store results
    ofstream output;

    output.open(avalpath);
    for (i=0; i < naval; i++)
    {
        //Initial conditions, single seed
        lattice = vector<bool>(N, 0);
        lattice[center] = 1;

        //Initialize counters
        t = 0;
        size = 1;
        rho = 1;

        //Run until absorbing state.
        //Each row has the temporal evolution of a single avalanche
        while (rho > 0)
        {
            output << rho << " ";
            rho = update_lattice(lattice, neighs);
            t++;

        }        
        output << endl; //For next avalanche
    }
    output.close();

    return;
}

