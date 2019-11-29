/* TRAN VAN HUNG
   Vietnam National University,
   University of Engineering and Technology,
   Engineering Physics and Nanotechnology */

#include <iostream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <fstream>

using namespace std;

//file to output data into
ofstream DATA("DATA.GaussHistt.txt",ios::out);

// Mersenne Twister algorithm

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

double genrand()
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

// box muller transform gaussian distribution
double boxmullergauss(double m, double s)
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
    
    return (m + y2 * s);
}

int main()
{
    int n = 10000;
    int *p;
    double min, max;
    
    //cout << "How many random numbers do you generate? \n Enter number of values n = ";
    
    double a[n]; // array to store random numbers
    
    cout << "\n*****************************";
    cout << "\nGenerating " << n << " random numbers\n";
    
    for(int i = 0; i < n; i++)
    {
        a[i] = boxmullergauss(0, 1);
    }
    
    cout << "\n\n\n********************************************";
    cout << "\nThong ke lai cac gia tri so ngau nhien da tao\n";
    
    double count[n];// mang luu bien dem
    max = a[0];
    min = a[0];
        
    for (int i = 0; i < sizeof(a)/sizeof(double); i++)
    {
        if (max < a[i]) max = a[i];
        if (min > a[i]) min = a[i];
    }
    
    cout << "\n\nLength of array: " << sizeof(a)/sizeof(double) << endl;
    
    cout << "Maximum value: " << max << "\nMinimum value: " << min << endl;
    
    double dx = (max - min)/100; // chia khoang gia tri thanh 100 khoang nho de thong ke
    
    cout << "Interval length: dx = " << (max - min)/100 << endl;
    
    for (int j = 0; j < 100; j++)
    {
        int d = 0;
        
        for (int i = 0; i < n; i++)
        {
            if((a[i] >= min) & (a[i] < min + dx)) d++;
        }
        
        if (j == 99) d = d+1;
        
        count[j] = d;
        min += dx;
        
        cout << "Interval " << j + 1 << " have " << count[j] << " value" << endl << endl;
        DATA << min << "\t" << count[j] << endl;
    }
    
    system("pause");
}
