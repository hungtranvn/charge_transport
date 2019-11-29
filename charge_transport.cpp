/* TRAN VAN HUNG
   Vietnam National University,
   University of Engineering and Technology,
   Nano Materials and Devices */

// Charge transport in polymer (disordered) materials

#include <iostream>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <cstdlib> 
#include <ctime>

using namespace std;

//file to output data into
ofstream DATA("test1D.txt",ios::out);

const int size = 100; // kich thuoc mang
const int n = size; // number of sites on lattice
float position[size+1]; // 1d lattice for sites
float energy[size+1]; // // 1d lattice for energy of each site
const float d = 1.; // lattice const
const float v0 = 1.; // the attempt-to-escape frequency
const float a = 0.2; // the localized length

//maximum number of particles
#define NMAX 100000
//maximum number of nearest neighbors in rate matrix
#define BMAX 100

//for ran0(&idum)
long int idum;
//rate matrix
float rate[NMAX][BMAX];
//number of neighbors
int nind[NMAX];
//identity of neighbors
int rind[NMAX][BMAX];
//rate total
float ktot[NMAX];
//temperature
float kT = 0.6;
float maxkT = 0.7;
float kTchange = 0.1;
//electric field
float eF = 0.001;
float maxeF = 1;
float eFchange = 0.05;

// ran0
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

float ran0(long *idum)
{
    long k;
    float ans;
    *idum ^= MASK;
    k = (*idum) / IQ;
    *idum = IA*(*idum - k*IQ) - IR*k;
    if (*idum < 0)*idum += IM;
    ans = AM*(*idum);
    *idum ^= MASK;
    return ans;
}

/*period parameters*/
#define N 624
#define M 397
#define MATRIX_A 0x9908b0df /*constant vector a*/
#define UPPER_MASK 0x80000000 /*most significant w-r bits*/
#define LOWER_MASK 0x7fffffff /*least significant r bit*/

/*tempering parameters*/
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y) (y >> 11)
#define TEMPERING_SHIFT_S(y) (y << 7)
#define TEMPERING_SHIFT_T(y) (y << 15)
#define TEMPERING_SHIFT_L(y) (y >> 18)

static unsigned long mt[N]; //the array for the state vector
static int mti = N + 1;

// initializing the array with a NONZERO seed
//insigned long seed;

void sgenrand(unsigned long seed)
{
    mt[0] = seed & 0xffffffff;
    for (mti = 1; mti < N; mti++)
        mt[mti] = (69069 * mt[mti - 1]) & 0xffffffff;
}

float genrand()
{
    unsigned long y;
    static unsigned long mag01[2] = { 0x0, MATRIX_A }; //mag01[x] = x * MATRIX_A for x = 0,1
    if (mti >= N)
    {
        if (mti == N + 1)
            sgenrand(4357);
        for (int kk = 0; kk < N - M; kk++)
        {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (int kk; kk < N - 1; kk++)
        {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1];
        mti = 0;
    }
    
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);
    return ((double)y / (unsigned long)0xffffffff); //reals

}

// generate gaussian distribution
float gaussrand(float m, float s)
{
    float x1, x2, w, y1, y2;

    do {
        x1 = 2.0 * genrand() - 1.0;
        x2 = 2.0 * genrand() - 1.0;
        w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );

    if(w == 0)
    {
        y1 = y2 = 0;
    }
    else
    {
        w = sqrt( (-2.0 * log( w ) ) / w );
        y1 = x1 * w;
        y2 = x2 * w;
    }

    return m + y2 * s;
}

// generate exponential distribution
float exprand(double m)
{
    return -log(1 - genrand())/m;
}

//function for random initialization of lattice
void initialize()
{
    for(int i = 1; i <= size; i++)
    {
        position[i] = genrand()*size;
        energy[i] = gaussrand(0, 2);
    }
}

// ouput of lattice configuration to the screen
void output()
{
    for(int i = 1; i <= size; i++)
    {
        cout << "[" << i << "]= " << position[i] << "\t" << energy[i] << endl;
    }
    cout << endl;

}

get_rate(int index)
{
    float dr;
    float rsq,k;
    int count = 0;
    ktot[index] = 0.0;

    rind[index][count] = index;
    count++;
    
    for(int i = 1; i <= size; i++)
    {
        if(index != i)
        {
            //cout << i << endl;
            
            dr = position[i] - position[index];
            //periodic boundaries in x
            if(dr >= size/2.) dr -= size;
            if(dr < -size/2.) dr += size;
            
            rsq = dr*dr;
            
            //within cutoff radius?
            if(rsq <= 10)
            {
                float vij, delE, BTZ_factor;
                delE = energy[i] - energy[index] - eF*0.5*dr;

                if (delE > 0)
                {
                    BTZ_factor = exp(-delE/kT*0.5);
                }
                else
                {
                    BTZ_factor = 1.0;
                }
                
                //compute rate
                vij = exp(-10*sqrt(rsq))*BTZ_factor;
                
                rate[index][count] = vij;
                
                rind[index][count] = i; //luu cac vi tri co the nhay toi
                ktot[index] += rate[index][count];
                
                //cout << endl << i << "\t" << position[i] << "\t" << energy[i] << "\t" << rate[index][count] << endl;
                count++;
            }
        }
    }
    nind[index] = count;
    
    float temp;
    int h;
    int i, j;
    for(int i = 1; i < count-1; i++)
    {
        for(j = i+1; j < count; j++)
        {
            if(rate[index][i] > rate[index][j])
            {
                temp = rate[index][i];
                rate[index][i] = rate[index][j];
                rate[index][j] = temp;
                
                k = rind[index][i];
                rind[index][i] = rind[index][j];
                rind[index][j] = k;
            }
        }
        
    }
    
    //for(int i = 1; i < count; i++)
    //{
    //  cout << endl << rind[index][i] << "\t" << rate[index][i] << endl;
    //}
}


//main function
int main()
{
    long tstop = 5e7;
    idum = -time(NULL);
    
    initialize();
    output();
    int istate = 10;
    get_rate(istate);
    for(; eF <= maxeF; eF = eF + eFchange)
    {
        bool alive = true;
        int count;
        float time = 0;
        float ftemp[3];
        int state, newstate;
        state = istate;
        float msq = 0;
        while(alive)
        {
            //cout << state << "\t" << time << "\t\t" << position[state] << "\t\t" << energy[state] << "\t\t" << msq << endl;
            //DATA << state << "\t" << time << "\t\t" << position[state] << "\t\t" << energy[state] << "\t\t" << msq << endl;
            get_rate(state);
            ftemp[0] = genrand(); // tao so ngau nhien trong khoang 0-1 (xac suat ngau nhien)
            ftemp[1] = ktot[state]*ftemp[0]; // tinh phan tram xac suat tong
            ftemp[2] = 0.0; // gan xac suat ban dau bang 0
            count = 0;
            
            //identify next step
            while(ftemp[2] < ftemp[1])
            {
                count++;
                ftemp[2] += rate[state][count];
                //cout << ftemp[2] << "\t" << ftemp[1] << endl;
            }
            
            time += exprand(1)/ktot[state];
            newstate = rind[state][count];
            state = newstate;
            
            if (time > tstop)
            {
                alive = false;
            }
            
            //indentify mean square displacement
            msq += (position[state] - position[istate])*(position[state] - position[istate]);
        }
        cout << msq/tstop << endl;
        DATA << eF << "\t" << msq/tstop << endl;
    }
        
    system("pause");
    return 0;
}
