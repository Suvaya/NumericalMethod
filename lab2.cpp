//1
//#include<iostream>
//#include<conio.h>
//
//using namespace std;
//
//int main()
//{
//	 float x[100], y[100], xp, yp=0, p;
//	 int i,j,n;
//
//	 //User Input
//	 cout<<"Enter number of data: ";
//	 cin>>n;
//	 cout<<"Enter data:"<< endl;
//	 for(i=1;i<=n;i++)
//	 {
//		  cout<<"x["<< i<<"] = ";
//		  cin>>x[i];
//		  cout<<"y["<< i<<"] = ";
//		  cin>>y[i];
//	 }
//	 cout<<"Enter interpolation point: ";
//	 cin>>xp;
//
//	 //Implementing Lagrange Interpolation
//	 for(i=1;i<=n;i++)
//	 {
//		  p=1;
//		  for(j=1;j<=n;j++)
//		  {
//			   if(i!=j)
//			   {
//			    	p = p* (xp - x[j])/(x[i] - x[j]);
//			   }
//		  }
//		  yp = yp + p * y[i];
//	 }
//	 cout<< endl<<"Interpolated value at "<< xp<< " is "<< yp;
//
//	 return 0;
//}

// 2. CPP Program to interpolate using
//// newton forward interpolation
////#include <bits/stdc++.h>
//#include<iostream>
//#include<stdlib.h>
//#include <iomanip>
//using namespace std;
//
//// calculating u mentioned in the formula
//float u_cal(float u, int n)
//{
//	float temp = u;
//	for (int i = 1; i < n; i++)
//		temp = temp * (u - i);
//	return temp;
//}
//
//// calculating factorial of given number n
//int fact(int n)
//{
//	int f = 1;
//	for (int i = 2; i <= n; i++)
//		f *= i;
//	return f;
//}
//
//int main()
//{
//	// Number of values given
//	int n = 4;
//	//float x[] = { 45, 50, 55, 60 };
//	float x[n];
//
//	// y[][] is used for difference table
//	// with y[][0] used for input
//	float y[n][n];
//
//	cout << "Enter given values:" << endl;
//	for (int i = 0; i < n; i++) {
//		cin >> x[i];
//	}
//
//
//	cout << "Enter value for the differnece table:" << endl;
//	for (int i = 0; i < n; i++) {
//		cin >> y[i][0];
//	}
//
//	// Calculating the forward difference
//	// table
//	for (int i = 1; i < n; i++) {
//		for (int j = 0; j < n - i; j++)
//			y[j][i] = y[j + 1][i - 1] - y[j][i - 1];
//	}
//
//	// Displaying the forward difference table
//	for (int i = 0; i < n; i++) {
//		cout << setw(4) << x[i]
//			<< "\t";
//		for (int j = 0; j < n - i; j++)
//			cout << setw(4) << y[i][j]
//			<< "\t";
//		cout << endl;
//	}
//
//	// Value to interpolate at
//	float value = 52;
//
//	// initializing u and sum
//	float sum = y[0][0];
//	float u = (value - x[0]) / (x[1] - x[0]);
//	for (int i = 1; i < n; i++) {
//		sum = sum + (u_cal(u, i) * y[0][i]) /
//			fact(i);
//	}
//
//	cout << "\n Value at " << value << " is "
//		<< sum << endl;
//	return 0;
//}


// 3. CPP Program to interpolate using
//// newton backward interpolation
//#include <bits/stdc++.h>
//using namespace std;
//
//// Calculation of u mentioned in formula
//float u_cal(float u, int n)
//{
//	float temp = u;
//	for (int i = 1; i < n; i++)
//		temp = temp * (u + i);
//	return temp;
//}
//
//// Calculating factorial of given n
//int fact(int n)
//{
//	int f = 1;
//	for (int i = 2; i <= n; i++)
//		f *= i;
//	return f;
//}
//
//int main()
//{
//	// number of values given
//	int n = 5;
//	//float x[] = { 1891, 1901, 1911,1921, 1931 };
//	float x[n];
//	float y[n][n];
//
//	cout<<"Enter given values:"<<endl;
//	for(int i=0; i<n; i++){
//		cin>>x[i];
//	}
//
//				
//	// y[][] is used for difference
//	// table and y[][0] used for input
//	
//	// y[0][0] = 46;
//	// y[1][0] = 66;
//	// y[2][0] = 81;
//	// y[3][0] = 93;
//	// y[4][0] = 101;
//	cout<<"Enter value for the differnece table:"<<endl;
//	for(int i=0; i<n; i++){
//		cin>>y[i][0];
//	}
//
//	// Calculating the backward difference table
//	for (int i = 1; i < n; i++) {
//		for (int j = n - 1; j >= i; j--)
//			y[j][i] = y[j][i - 1] - y[j - 1][i - 1];
//	}
//
//	// Displaying the backward difference table
//	for (int i = 0; i < n; i++) {
//		for (int j = 0; j <= i; j++)
//			cout << setw(4) << y[i][j]
//				<< "\t";
//		cout << endl;
//	}
//
//	// Value to interpolate at
//	float value = 1925;
//
//	// Initializing u and sum
//	float sum = y[n - 1][0];
//	float u = (value - x[n - 1]) / (x[1] - x[0]);
//	for (int i = 1; i < n; i++) {
//		sum = sum + (u_cal(u, i) * y[n - 1][i]) /
//									fact(i);
//	}
//
//	cout << "\n Value at " << value << " is "
//		<< sum << endl;
//	return 0;
//}




// 4. CPP program for implementing
//// Newton divided difference formula
//#include <bits/stdc++.h>
//using namespace std;
//
//// Function to find the product term
//float proterm(int i, float value, float x[])
//{
//	float pro = 1;
//	for (int j = 0; j < i; j++) {
//		pro = pro * (value - x[j]);
//	}
//	return pro;
//}
//
//// Function for calculating
//// divided difference table
//void dividedDiffTable(float x[], float y[][10], int n)
//{
//	for (int i = 1; i < n; i++) {
//		for (int j = 0; j < n - i; j++) {
//			y[j][i] = (y[j][i - 1] - y[j + 1]
//						[i - 1]) / (x[j] - x[i + j]);
//		}
//	}
//}
//
//// Function for applying Newton's
//// divided difference formula
//float applyFormula(float value, float x[],
//				float y[][10], int n)
//{
//	float sum = y[0][0];
//
//	for (int i = 1; i < n; i++) {
//	sum = sum + (proterm(i, value, x) * y[0][i]);
//	}
//	return sum;
//}
//
//// Function for displaying
//// divided difference table
//void printDiffTable(float y[][10],int n)
//{
//	for (int i = 0; i < n; i++) {
//		for (int j = 0; j < n - i; j++) {
//			cout << setprecision(4) <<
//								y[i][j] << "\t ";
//		}
//		cout << "\n";
//	}
//}
//
//// Driver Function
//int main()
//{
//	// number of inputs given
//	int n = 4;
//	float value, sum, y[10][10];
//	float x[] = { 5, 6, 9, 11 };
//
//	// y[][] is used for divided difference
//	// table where y[][0] is used for input
//	y[0][0] = 12;
//	y[1][0] = 13;
//	y[2][0] = 14;
//	y[3][0] = 16;
//
//	// calculating divided difference table
//	dividedDiffTable(x, y, n);
//
//	// displaying divided difference table
//	printDiffTable(y,n);
//
//	// value to be interpolated
//	value = 7;
//
//	// printing the value
//	cout << "\nValue at " << value << " is "
//			<< applyFormula(value, x, y, n) << endl;
//	return 0;
//}
//


//5
//#include <iostream>
//#include <math.h>
//using namespace std;
//
//class cubic
//{
//public:
//    void askAll();
//    void askI();
//    void askX();
//    void askFi();
//    void askXi();
//
//    void findA();
//    void findH();
//    void findS();
//    void findU();
//
//    void solve();
//
//
//private:
//    double xi[20],f[20],h[19],s[19],a[20],u[20];
//    int Inumber;
//    double XX;
//
//};
//void cubic::askX()
//{
//    cout << "Enter X: ";
//    cin >> XX;
//}
//void cubic::askI()
//{
//    cout << "Enter number of Items: ";
//    cin >> Inumber;
//}
//
//void cubic::askXi()
//{
//    for(int i = 0; i < Inumber; i++)
//    {
//        cout << "Enter X" << i << ": ";
//        cin >> xi[i];
//    }
//}
//
//
//void cubic::askFi()
//{
//    for(int i = 0; i < Inumber; i++)
//    {
//        cout << "Enter F" << i << ": ";
//        cin >> f[i];
//    }
//}
//
//void cubic::findU()
//{
//    for(int i = 0; i < Inumber; i++)
//    {
//        u[i] = XX - xi[i];
//    }
//}
//
//void cubic::findH()
//{
//    for(int i = 1; i < Inumber; i++)
//    {
//        h[i] = xi[i] - xi[i-1];
//    }
//}
//void cubic::findA()
//{
//    for(int i = 1; i < Inumber-1; i++)
//    {
//        a[i] = ( 6*( (f[i+1]-f[i])/h[i+1] -  (f[i]-f[i-1])/h[i]  )  -  h[i]*a[i-1]  -  
//               h[i+1]*a[i+1]  ) / (2*( h[i] + h[i+1] ));
//    }
//}
//
//void cubic::findS()
//{
//    int k;
//    for(int i = 1; i < Inumber; i++)
//    {
//        if(XX>=xi[i-1]&& XX<=xi[i]){
//        s[i] =  a[i-1]/(6*h[i]) * ( pow(h[i],2)*u[i] - pow(u[i],2) )  +  a[i]/(6*h[i]) * 
//                ( pow(u[i-1],3)  - pow(h[i],2)*u[i-1]  )   +  1/h[i] * (f[i]*u[i-1]  -  
//                f[i-1]*u[i]);
//        k = i;
//    }
//    }
//    cout << "s "<< s[k];
//
//}
//
//void cubic::askAll()
//{
//    askI();
//    askXi();
//    askFi();
//    askX();
//    a[0] = 0;
//    a[Inumber-1] = 0;
//}
//
//void cubic::solve()
//{
//    findH();
//    findU();
//    findA();
//    findS();
//}
//
//int main()
//{
//
//    cubic c1;
//    c1.askAll();
//    c1.solve();
//    return 0;
//}




//7
//#include<iostream>
//#include<iomanip>
//#include<cmath>
//using namespace std;
//int main()
//{
//    int i,j,k,n;
//    cout<<"\nEnter the no. of data pairs to be entered:\n";        //To find the size of arrays
//    cin>>n;
//    double x[n],y[n],a,b;
//    cout<<"\nEnter the x-axis values:\n";                //Input x-values
//    for (i=0;i<n;i++)
//        cin>>x[i];
//    cout<<"\nEnter the y-axis values:\n";                //Input y-values
//    for (i=0;i<n;i++)
//        cin>>y[i];
//    double xsum=0,x2sum=0,ysum=0,xysum=0;                //variables for sums/sigma of xi,yi,xi^2,xiyi etc
//    for (i=0;i<n;i++)
//    {
//        xsum=xsum+x[i];                        //calculate sigma(xi)
//        ysum=ysum+y[i];                        //calculate sigma(yi)
//        x2sum=x2sum+pow(x[i],2);                //calculate sigma(x^2i)
//        xysum=xysum+x[i]*y[i];                    //calculate sigma(xi*yi)
//    }
//    a=(n*xysum-xsum*ysum)/(n*x2sum-xsum*xsum);            //calculate slope
//    b=(x2sum*ysum-xsum*xysum)/(x2sum*n-xsum*xsum);            //calculate intercept
//    double y_fit[n];                        //an array to store the new fitted values of y    
//    for (i=0;i<n;i++)
//        y_fit[i]=a*x[i]+b;                    //to calculate y(fitted) at given x points
//    cout<<"S.no"<<setw(5)<<"x"<<setw(19)<<"y(observed)"<<setw(19)<<"y(fitted)"<<endl;
//    cout<<"-----------------------------------------------------------------\n";
//    for (i=0;i<n;i++)
//        cout<<i+1<<"."<<setw(8)<<x[i]<<setw(15)<<y[i]<<setw(18)<<y_fit[i]<<endl;//print a table of x,y(obs.) and y(fit.)    
//    cout<<"\nThe linear fit line is of the form:\n\n"<<a<<"x + "<<b<<endl;        //print the best fit line
//    return 0;
//}  



//9. Polynomial Fit
//#include<iostream>
//#include<iomanip>
//#include<cmath>
//using namespace std;
//int main()
//{
//    int i,j,k,n,N;
//    cout.precision(4);                        //set precision
//    cout.setf(ios::fixed);
//    cout<<"\nEnter the no. of data pairs to be entered:\n";        //To find the size of arrays that will store x,y, and z values
//    cin>>N;
//    double x[N],y[N];
//    cout<<"\nEnter the x-axis values:\n";                //Input x-values
//    for (i=0;i<N;i++)
//        cin>>x[i];
//    cout<<"\nEnter the y-axis values:\n";                //Input y-values
//    for (i=0;i<N;i++)
//        cin>>y[i];
//    cout<<"\nWhat degree of Polynomial do you want to use for the fit?\n";
//    cin>>n;                                // n is the degree of Polynomial 
//    double X[2*n+1];                        //Array that will store the values of sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
//    for (i=0;i<2*n+1;i++)
//    {
//        X[i]=0;
//        for (j=0;j<N;j++)
//            X[i]=X[i]+pow(x[j],i);        //consecutive positions of the array will store N,sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
//    }
//    double B[n+1][n+2],a[n+1];            //B is the Normal matrix(augmented) that will store the equations, 'a' is for value of the final coefficients
//    for (i=0;i<=n;i++)
//        for (j=0;j<=n;j++)
//            B[i][j]=X[i+j];            //Build the Normal matrix by storing the corresponding coefficients at the right positions except the last column of the matrix
//    double Y[n+1];                    //Array to store the values of sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
//    for (i=0;i<n+1;i++)
//    {    
//        Y[i]=0;
//        for (j=0;j<N;j++)
//        Y[i]=Y[i]+pow(x[j],i)*y[j];        //consecutive positions will store sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
//    }
//    for (i=0;i<=n;i++)
//        B[i][n+1]=Y[i];                //load the values of Y as the last column of B(Normal Matrix but augmented)
//    n=n+1;                //n is made n+1 because the Gaussian Elimination part below was for n equations, but here n is the degree of polynomial and for n degree we get n+1 equations
//    cout<<"\nThe Normal(Augmented Matrix) is as follows:\n";    
//    for (i=0;i<n;i++)            //print the Normal-augmented matrix
//    {
//        for (j=0;j<=n;j++)
//            cout<<B[i][j]<<setw(16);
//        cout<<"\n";
//    }    
//    for (i=0;i<n;i++)                    //From now Gaussian Elimination starts(can be ignored) to solve the set of linear equations (Pivotisation)
//        for (k=i+1;k<n;k++)
//            if (B[i][i]<B[k][i])
//                for (j=0;j<=n;j++)
//                {
//                    double temp=B[i][j];
//                    B[i][j]=B[k][j];
//                    B[k][j]=temp;
//                }
//     
//    for (i=0;i<n-1;i++)            //loop to perform the gauss elimination
//        for (k=i+1;k<n;k++)
//            {
//                double t=B[k][i]/B[i][i];
//                for (j=0;j<=n;j++)
//                    B[k][j]=B[k][j]-t*B[i][j];    //make the elements below the pivot elements equal to zero or elimnate the variables
//            }
//    for (i=n-1;i>=0;i--)                //back-substitution
//    {                        //x is an array whose values correspond to the values of x,y,z..
//        a[i]=B[i][n];                //make the variable to be calculated equal to the rhs of the last equation
//        for (j=0;j<n;j++)
//            if (j!=i)            //then subtract all the lhs values except the coefficient of the variable whose value                                   is being calculated
//                a[i]=a[i]-B[i][j]*a[j];
//        a[i]=a[i]/B[i][i];            //now finally divide the rhs by the coefficient of the variable to be calculated
//    }
//    cout<<"\nThe values of the coefficients are as follows:\n";
//    for (i=0;i<n;i++)
//        cout<<"x^"<<i<<"="<<a[i]<<endl;            // Print the values of x^0,x^1,x^2,x^3,....    
//    cout<<"\nHence the fitted Polynomial is given by:\ny=";
//    for (i=0;i<n;i++)
//        cout<<" + ("<<a[i]<<")"<<"x^"<<i;
//    cout<<"\n";
//    return 0;
//}