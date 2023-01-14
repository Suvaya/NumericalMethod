#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
using namespace std;
/*
// Function to calculate the divided difference table
vector<vector<double>> dividedDifferenceTable(vector<double> x, vector<double> y) {
    int n = x.size();
    vector<vector<double>> table(n, vector<double>(n));
    for (int i = 0; i < n; i++) {
        table[i][0] = y[i];
    }
    for (int j = 1; j < n; j++) {
        for (int i = j; i < n; i++) {
            table[i][j] = (table[i][j-1] - table[i-1][j-1]) / (x[i] - x[i-j]);
        }
    }
    return table;
}

// Function to calculate the derivative of a function using the divided difference table
double derivative(vector<double> x, vector<double> y, double x0) {
    int n = x.size();
    vector<vector<double>> table = dividedDifferenceTable(x, y);
    double result = 0, term = 1;
    for (int i = 1; i < n; i++) {
        term *= (x0 - x[i-1]);
        result += term * table[i][i];
    }
    return result;
}

int main() {
    // Input data (x, y)
    vector<double> x = {1, 2, 3, 4, 5};
    vector<double> y = {2, 4, 5, 4, 5};

    // Point at which to calculate the derivative
    double x0 = 3.5;

    // Calculate the derivative
    double deriv = derivative(x, y, x0);
    cout << "Derivative at x = " << x0 << ": " << deriv << endl;

    return 0;
}
*/

/*
Trapezoidal rule program
#include<iostream>
#include<math.h>


#define f(x) 1/(1+pow(x,2))

float trap(float a, float b){
    float temp = (f(a)+f(b))/2;
    float res=(b-a)*temp;
    return res;
}

using namespace std;
int main()
{
 float lower, upper, integration=0.0, stepSize, k;
 int i, subInterval;

 cout<<"Enter lower limit of integration: ";
 cin>>lower;
 cout<<"Enter upper limit of integration: ";
 cin>>upper;


 cout<< endl<<"Required value of integration is: "<<trap(lower, upper)<<endl;

 return 0;
}
*/

/*
#define f(x) 1/(1+pow(x,2))

//composite trapezoidal rule
float trap(float a, float b, float n){
    float sum=0;
    float h = (b-a)/n;
    float integ = f(a)+f(b);
    for(int i=1; i<=n-1; i++){
        float steps = a+i*h;
        sum += 2*f(steps);
    }
    integ = (integ+sum)*(h/2);
    return integ;
}


int main()
{
 float lower, upper, integration=0.0, stepSize, k;
 int i, subInterval;

 cout<<"Enter lower limit of integration: ";
 cin>>lower;
 cout<<"Enter upper limit of integration: ";
 cin>>upper;
 cout<<"Enter number of sub intervals: ";
 cin>>subInterval;

 cout<< endl<<"Required value of integration is: "<< integration<<endl<<trap(lower, upper, subInterval);

 return 0;
}
*/

/*

#define f(x) 1/(1+pow(x,2))
//for simpson 1/3 where n=2
float simpson(float a, float b){
    float h=(b-a)/2;
    float sum =f(a)+4*f((a+b)/2)+f(b);
    return (sum*h)/3;

}

float composite_simpson(float a, float b,int n)
{
    float h = (b - a) / n;

    // Internal sample points, there should be n - 1 of them
    float sum_odds = 0.0;
    for (int i = 1; i < n; i += 2)
    {
        sum_odds += f(a + i * h);
    }
    float sum_evens = 0.0;
    for (int i = 2; i < n; i += 2)
    {
        sum_evens += f(a + i * h);
    }

    return (f(a) + f(b) + 2 * sum_evens + 4 * sum_odds) * h / 3;
}

int main()
{
 float lower, upper, integration=0.0, stepSize, k;
 int i, subInterval;

 cout<<"Enter lower limit of integration: ";
 cin>>lower;
 cout<<"Enter upper limit of integration: ";
 cin>>upper;
 cout<<"Enter number of sub intervals: ";
 cin>>subInterval;
 

 cout<< endl <<"Required value of integration is: "<<simpson(lower, upper);
 cout<< endl <<"Required value of integration using composite simpson is: "<<composite_simpson(lower,upper,subInterval)<<endl;

 return 0;
}

*/

/*

#define f(x) 1/(1+pow(x,2))

//simpson 1/3 rule for n sub intervals
float simpson_rule(float a, float b,int n)
{
    float h = (b - a) / n;

    // Internal sample points, there should be n - 1 of them
    float sum_odds = 0.0;
    for (int i = 1; i < n; i += 2)
    {
        sum_odds += f(a + i * h);
    }
    float sum_evens = 0.0;
    for (int i = 2; i < n; i += 2)
    {
        sum_evens += f(a + i * h);
    }

    return (f(a) + f(b) + 2 * sum_evens + 4 * sum_odds) * h / 3;
}

using namespace std;
int main()
{
 float lower, upper, integration=0.0, stepSize, k;
 int i, subInterval;

 cout<<"Enter lower limit of integration: ";
 cin>>lower;
 cout<<"Enter upper limit of integration: ";
 cin>>upper;
 cout<<"Enter number of sub intervals: ";
 cin>>subInterval;


 cout<< endl <<"Required value of integration is: "<<simpson_rule(lower,upper,subInterval)<<endl;

 return 0;
}

*/


/*

#define f(x) 1/(1+pow(x,2))

//simpson 3/8 rule for n=3
float simpson(float a, float b){
    float h=(b-a)/3;
    float sum = f(a) + f(b);

    for(int i=1; i<=3-1; i++){
        float k = a+i*h;
        sum+=3*f(k);
    }
    return sum*(3*h/8);
}


using namespace std;
int main()
{
 float lower, upper;

 cout<<"Enter lower limit of integration: ";
 cin>>lower;
 cout<<"Enter upper limit of integration: ";
 cin>>upper;

 cout<< endl <<"Required value of integration is: "<<simpson(lower, upper);

 return 0;
}

*/

/*

#define f(x) 1/(1+pow(x,2))

//function for simpson 3/8 rule
float simpson(float a, float b){
    float h=(b-a)/3;
    float sum = f(a) + f(b);

    for(int i=1; i<=3-1; i++){
        float k = a+i*h;
        sum+=3*f(k);
    }
    return sum*(3*h/8);
}

//function for composite simpson 3/8 rule
float Csimpson(float a, float b, float n){
    float h=(b-a)/n;
    float sum = f(a) + f(b);
    float sumD3 = 0;
    float sumND3 = 0;

    for(int i=1; i<=n-1; i++){
        float k = a+i*h;

        if(i%3==0){
            sumD3 += f(k);
        }
        else{
            sumND3 += f(k);
        }

    }
    return (sum+3*sumND3+2*sumD3)*(3*h/8);
}


int main()
{
 float lower, upper, integration=0.0, stepSize, k;
 int i, subInterval;

 cout<<"Enter lower limit of integration: ";
 cin>>lower;
 cout<<"Enter upper limit of integration: ";
 cin>>upper;
 cout<<"Enter number of sub intervals: ";
 cin>>subInterval;

 cout<< endl <<"Required value of integration using simple simpson and composite simpson 3/8 rule respectively is: "<<simpson(lower, upper)<<" and "<<Csimpson(lower, upper, subInterval);

 return 0;
}

*/


/*

#define f(x) pow(x,2)+1

//the intervals should be -1 to 1 if not operate by using substitution method and find new function


//for n=2;
float gauss2(){
    float w1=1, w2=1, x1=-1/sqrt(3), x2=1/sqrt(3);
    float integration = w1*f(x1)+w2*f(x2);
    return integration;
}

//for n=3
float gauss3(){
    float w1=5/9, w2=8/9, w3=5/9, x1=-sqrt(3/5), x2=0, x3=sqrt(3/5);
    float integration = w1*f(x1)+w2*f(x2)+w3*f(x3);
    return integration;
}


int main(){
 cout<<gauss2()<<endl<<gauss3();
 return 0;
}

*/

/*
int N = 5;

double romberg(double (*func)(double), double a, double b) {
    double h[N+1], r[N+1][N+1];
    for (int i = 1; i < N + 1; ++i) {
        h[i] = (b - a) / pow(2, i - 1);
    }
    r[1][1] = h[1] / 2 * (func(a) + func(b));
    for (int i = 2; i < N + 1; ++i) {
        double coeff = 0;
        for (int k = 1; k <= pow(2, i - 2); ++k) {
            coeff += func(a + (2 * k - 1) * h[i]);
        }
        r[i][1] = 0.5 * (r[i - 1][1] + h[i - 1] * coeff);
    }
    
    for (int i = 2; i < N + 1; ++i) {
        for (int j = 2; j <= i; ++j) {
            r[i][j] = r[i][j - 1] + (r[i][j - 1] - r[i - 1][j - 1]) / (pow(4, j - 1) - 1);
        }
    }
    return r[N][N];
}

double f(double x) {
    return 1/x;
}

int main()
{
    cout << romberg(f, 1, 10) << endl;
}
*/

/*

const int MAX_INTERVALS = 100;

double func(double x, double y) {
    return x * y;
}

double Dintegration(double x0, double x1, double y0, double y1, int intervals) {
    double dx = (x1 - x0) / intervals;
    double dy = (y1 - y0) / intervals;

    double result = 0;
    for (int i = 0; i < intervals; i++) {
        for (int j = 0; j < intervals; j++) {
            double x = x0 + i * dx;
            double y = y0 + j * dy;
            result += func(x, y) + func(x + dx, y) + func(x, y + dy) + func(x + dx, y + dy);
        }
    }
    result = result * dx * dy / 4;
    return result;
}

int main() {
    double x0, x1, y0, y1;
    cout << "Enter the limits of integration in the form x0 x1 y0 y1: ";
    cin >> x0 >> x1 >> y0 >> y1;
    int intervals;
    cout << "Enter the number of intervals: ";
    cin >> intervals;

    double result = Dintegration(x0, x1, y0, y1, intervals);
    cout << "The result of the double integration is: " << result << endl;

    return 0;
}

*/

