
//                  IBNM                                     //
//   An implementation of a branching-like process in 2D     //
//   Original implemtantion by Roxana Zeraati, translation   //
//   to C++, Victor Buendia                                  //
//   Last revision March 2021                                //
// --------------------------------------------------------- //

#include<iostream>
#include<cstdlib>
#include<vector>
#include<random>
#include<fstream>

using namespace std;

// --- Preprocessor magic ---

//4 or 8 neighbours
#define MOORE 0
#define MANHATTAN 1

#ifndef LATTICE
#define LATTICE MOORE
#endif

// --- Global constants ---

int L = 64;
int N = L*L;    // System size N=L*L

int k;          //Neighbourhood size
int num_neighs;

double p_rewire; //Probability to rewire


string filename;    //Where to store results

//RNG
mt19937 gen(151234553);    
uniform_real_distribution<double> ran_u(0.0, 1.0);  

// --- Function declarations ---

int compute_numneighs();
int manh_distance(const int x1, const int y1, const int x2, const int y2);
void create_neighbour_table(const string filename);

// --- Main call ---

int main(int argc, char *argv[])
{
    //State of the system and neighbour list
    vector<bool> lattice;

    vector<double> L_list = {16,32,64,128};
    vector<double> p_list = {0.001, 0.01, 0.1, 0.3};
    vector<double> k_list = {1,2,3,4,5,6,7};

    int j;

    k=1;
    p_rewire=0.0;
    cout << "building L" << endl;
    for (j=0; j < L_list.size(); j++)
    {
        L = L_list[j];

        N = L*L;
        num_neighs = compute_numneighs();
        filename = "../networks/L" + to_string(L);
        create_neighbour_table(filename);
    }

    L=64;
    N = L*L;
    k=1;
    cout << "building p" << endl;
    for (j=0; j < p_list.size(); j++)
    {
        p_rewire = p_list[j];
        filename = "../networks/p" + to_string(j);
        num_neighs = compute_numneighs();
        create_neighbour_table(filename);
    }


    L=64;
    N = L*L;
    p_rewire=0.0;
    cout << "building k" << endl;
    for (j=0; j < k_list.size(); j++)
    {
        k = k_list[j];
        filename = "../networks/k" + to_string(k);
        num_neighs = compute_numneighs();
        create_neighbour_table(filename);
    }

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

void create_neighbour_table(const string filename)
{
    int i,j,x,y;    //Count indices
    int index;      //Selected node
    int nindex;     //Selected neighbour
    int h,v;        //To select horizontal and vertical with respect to node

    int random_neigh; //Acount for a random index
    bool already_connected;
    int curr_neighs;

    vector< vector<int>> neighs = vector<vector<int>>(N, vector<int>(num_neighs)); 

    ofstream output(filename);

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
                        output << index << " " << h + v*L << endl; 
                        nindex++;
                        #else
                        if (manh_distance(x,y, h,v) <= k)
                        {
                            neighs[index][nindex] = h + v*L;
                            output << index << " " << h + v*L << endl; 
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
                        output << index << " " << random_neigh << endl;
                        nindex++;
                    }
                }
            } //end neighbourhood 
        }
    } //end lattice sweep

    output.close();
    return;
}