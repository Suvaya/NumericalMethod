#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

/*
QN 1.
// Function to find the product term
float proterm(int i, float value, float x[])
{
	float pro = 1;
	for (int j = 0; j < i; j++) {
		pro = pro * (value - x[j]);
	}
	return pro;
}

// Function for calculating
// divided difference table
void dividedDiffTable(float x[], float y[][10], int n)
{
	for (int i = 1; i < n; i++) {
		for (int j = 0; j < n - i; j++) {
			y[j][i] = (y[j][i - 1] - y[j + 1][i - 1]) / (x[j] - x[i + j]);
		}
	}
}

// Function for applying Newton's
// divided difference formula
float applyFormula(float value, float x[],
				float y[][10], int n)
{
	float sum = y[0][0];

	for (int i = 1; i < n; i++) {
	sum = sum + (proterm(i, value, x) * y[0][i]);
	}
	return sum;
}

// Function for displaying
// divided difference table
void printDiffTable(float y[][10],int n)
{
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n - i; j++) {
			cout << setprecision(4) <<y[i][j] << "\t ";
		}
		cout << "\n";
	}
}

// Driver Function
int main()
{
	// number of inputs given

	float value, sum, y[10][10];
	// float x[] = { 5, 6, 9, 11 };
    float x[4];
    int n= sizeof(x)/sizeof(int);

    cout<<"Enter the value for x:"<<endl;
    for(int i=0; i<n; i++){
        int k=i;
        cout<<++k<<":";
        cin>>x[i];
    };

	// y[][] is used for divided difference
	// table where y[][0] is used for input
    cout<<"Enter the vales for y: "<<endl;
    for(int i=0; i<n; i++){
        int k=i;
        cout<<++k<<":";
        cin>>y[i][0];
    }
	// y[0][0] = 12;
	// y[1][0] = 13;
	// y[2][0] = 14;
	// y[3][0] = 16;

    cout<<endl<<"Divided Difference Table"<<endl;

	// calculating divided difference table
	dividedDiffTable(x, y, n);

	// displaying divided difference table
	printDiffTable(y,n);

	// value to be interpolated
	value = 7;

	// printing the value
	cout << "\nValue at " << value << " is "
			<< applyFormula(value, x, y, n) << endl;
	return 0;
}
*/

/*
Qn 2.
/*
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


// //Qn 3
// #define f(x) 1/(1+pow(x,2))

// //composite trapezoidal rule
// float trap(float a, float b, float n){
//     float sum=0;
//     float h = (b-a)/n;
//     float integ = f(a)+f(b);
//     for(int i=1; i<=n-1; i++){
//         float steps = a+i*h;
//         sum += 2*f(steps);
//     }
//     integ = (integ+sum)*(h/2);
//     return integ;
// }


// int main()
// {
//  float lower, upper, integration=0.0, stepSize, k;
//  int i, subInterval;

//  /* Input */
//  cout<<"Enter lower limit of integration: ";
//  cin>>lower;
//  cout<<"Enter upper limit of integration: ";
//  cin>>upper;
//  cout<<"Enter number of sub intervals: ";
//  cin>>subInterval;

//  cout<< endl<<"Required value of integration is: "<< integration<<endl<<trap(lower, upper, subInterval);

//  return 0;
// }

//Qn 4
// #define f(x) 1/(1+pow(x,2))
// //for simpson 1/3 where n=2
// float simpson(float a, float b){
//     float h=(b-a)/2;
//     float sum =f(a)+4*f((a+b)/2)+f(b);
//     return (sum*h)/3;

// }

// int main()
// {
//  float lower, upper, integration=0.0, stepSize, k;
//  int i, subInterval;

//  /* Input */
//  cout<<"Enter lower limit of integration: ";
//  cin>>lower;
//  cout<<"Enter upper limit of integration: ";
//  cin>>upper;


//  cout<< endl <<"Required value of integration is: "<<simpson(lower, upper);

//  return 0;
// }



// Qn 5
// /* Define function here */
// #define f(x) 1/(1+pow(x,2))

// //simpson 1/3 rule for n sub intervals
// float simpson_rule(float a, float b,int n)
// {
//     float h = (b - a) / n;

//     // Internal sample points, there should be n - 1 of them
//     float sum_odds = 0.0;
//     for (int i = 1; i < n; i += 2)
//     {
//         sum_odds += f(a + i * h);
//     }
//     float sum_evens = 0.0;
//     for (int i = 2; i < n; i += 2)
//     {
//         sum_evens += f(a + i * h);
//     }

//     return (f(a) + f(b) + 2 * sum_evens + 4 * sum_odds) * h / 3;
// }


// int main()
// {
//  float lower, upper, integration=0.0, stepSize, k;
//  int i, subInterval;

//  /* Input */
//  cout<<"Enter lower limit of integration: ";
//  cin>>lower;
//  cout<<"Enter upper limit of integration: ";
//  cin>>upper;
//  cout<<"Enter number of sub intervals: ";
//  cin>>subInterval;


//  cout<< endl <<"Required value of integration is: "<<simpson_rule(lower,upper,subInterval)<<endl;

//  return 0;
// }


//Qn6
// #define f(x) 1/(1+pow(x,2))

// //simpson 3/8 rule for n=3
// float simpson(float a, float b){
//     float h=(b-a)/3;
//     float sum = f(a) + f(b);

//     for(int i=1; i<=3-1; i++){
//         float k = a+i*h;
//         sum+=3*f(k);
//     }
//     return sum*(3*h/8);
// }

// int main()
// {
//  float lower, upper;

//  /* Input */
//  cout<<"Enter lower limit of integration: ";
//  cin>>lower;
//  cout<<"Enter upper limit of integration: ";
//  cin>>upper;

//  cout<< endl <<"Required value of integration is: "<<simpson(lower, upper);

//  return 0;
// }


//Qn7
// #define f(x) 1/(1+pow(x,2))

// float simpson(float a, float b, float n){
//     float h=(b-a)/n;
//     float sum = f(a) + f(b);
//     float sumD3 = 0;
//     float sumND3 = 0;

//     for(int i=1; i<=n-1; i++){
//         float k = a+i*h;

//         if(i%3==0){
//             sumD3 += f(k);
//         }
//         else{
//             sumND3 += f(k);
//         }

//     }
//     return (sum+3*sumND3+2*sumD3)*(3*h/8);
// }


// int main()
// {
//  float lower, upper, integration=0.0, stepSize, k;
//  int i, subInterval;

//  /* Input */
//  cout<<"Enter lower limit of integration: ";
//  cin>>lower;
//  cout<<"Enter upper limit of integration: ";
//  cin>>upper;
//  cout<<"Enter number of sub intervals: ";
//  cin>>subInterval;

//  cout<< endl <<"Required value of integration is: "<<simpson(lower, upper, subInterval);

//  return 0;
// }


//QN8
// #define f(x) pow(x,2)+1
// //the intervals should be -1 to 1 if not operate by using substitution method and find new function


// //for n=2;
// float gauss2(){
//     float w1=1, w2=1, x1=-1/sqrt(3), x2=1/sqrt(3);
//     float integration = w1*f(x1)+w2*f(x2);
//     return integration;
// }

// //for n=3
// float gauss3(){
//     float w1=5/9, w2=8/9, w3=5/9, x1=-sqrt(3/5), x2=0, x3=sqrt(3/5);
//     float integration = w1*f(x1)+w2*f(x2)+w3*f(x3);
//     return integration;
// }





// int main(){
//  cout<<gauss2()<<endl<<gauss3();
//  return 0;
// }


//Qn 9
// int N = 5;

// double romberg(double (*func)(double), double a, double b) {
//     double h[N+1], r[N+1][N+1];
//     for (int i = 1; i < N + 1; ++i) {
//         h[i] = (b - a) / pow(2, i - 1);
//     }
//     r[1][1] = h[1] / 2 * (func(a) + func(b));
//     for (int i = 2; i < N + 1; ++i) {
//         double coeff = 0;
//         for (int k = 1; k <= pow(2, i - 2); ++k) {
//             coeff += func(a + (2 * k - 1) * h[i]);
//         }
//         r[i][1] = 0.5 * (r[i - 1][1] + h[i - 1] * coeff);
//     }

//     for (int i = 2; i < N + 1; ++i) {
//         for (int j = 2; j <= i; ++j) {
//             r[i][j] = r[i][j - 1] + (r[i][j - 1] - r[i - 1][j - 1]) / (pow(4, j - 1) - 1);
//         }
//     }
//     return r[N][N];
// }

// double f(double x) {
//     return 1/x;
// }

// int main()
// {
//     cout << romberg(f, 1, 10) << endl;
// }




