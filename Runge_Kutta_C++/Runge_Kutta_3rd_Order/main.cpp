#include <iostream>
#include <vector>
#include <string>
#include <fstream>

//This function creates a truple class to return from a function
auto truple(double y1, double y2, double y3) {
  struct retVals {        // Declare a local structure
    double i1, i2, i3;
  };
  return retVals {y1,y2,y3}; // Return the local structure
}

//defining the 3st order ODE to be solved by converting the equation to a linear one
// y1 := y'
double f(  double a, double b, double c){
    return b;
}

// y2 := y''= y1'
double h(  double a, double b, double c){
    return c;
}

// y3 := y''' = -yy''+ 1 - y'^2.
double g( double a, double b, double c){
    return -a*c + b*b - 1.;
}

//implements an update according to Runge-Kutta
auto rk_update( double t, double a,double b, double c, double dt ){
    
    double k11,k12,k13,k21,k22,k23,k31,k32,k33,k41,k42,k43;
    double ti=t;
    
    k11 = dt*f(a, b, c);
    k12 = dt*h(a, b, c);
    k13 = dt*g(a, b, c);
    
    k21 = dt*f(a+0.5*k11, b+0.5*k12, c+0.5*k13);
    k22 = dt*h(a+0.5*k11, b+0.5*k12, c+0.5*k13);
    k23 = dt*g(a+0.5*k11, b+0.5*k12, c+0.5*k13);
    
    k31 = dt*f(a+0.5*k21, b+0.5*k22, c+0.5*k23);
    k32 = dt*h(a+0.5*k21, b+0.5*k22, c+0.5*k23);
    k33 = dt*g(a+0.5*k21, b+0.5*k22, c+0.5*k23);
    
    k41 = dt*f(a+k31, b+k32, c+k33);
    k42 = dt*h(a+k31, b+k32, c+k33);
    k43 = dt*g(a+k31, b+k32, c+k33);
    
    a = a + 1/6.0*(k11+(2.0*k21)+(2.0*k31)+k41);
    b = b + 1/6.0*(k12+(2.0*k22)+(2.0*k32)+k42);
    c = c + 1/6.0*(k13+(2.0*k23)+(2.0*k33)+k43);
    
    return truple(a,b,c);

}


int main()
{
    //Initial conditions
    double y0dd =  0., y0d = 0., y0 = 0., t0 = 0., tf = 4., dt = 0.0001;
    double a = y0, b = y0d, c = y0dd, t = t0;
    
    std::ofstream outfile ("graph3.csv"); //This writes a CSV file to be plotted
    outfile << "t" <<","<< "y" << std::endl; //Writes headings

    
    //loops over 1 RK update
    while(t <= tf){
        
        outfile << t <<","<< y0 << std::endl; //This writes a CSV file to be plotted
        auto[a, b, c] = rk_update(t, y0, y0d, y0dd, dt );
        y0dd =  c; y0d = b; y0 = a;
        t += dt;

    }
    
    outfile.close();
    
    std::cout << "x(t="<< t <<")="<< y0 <<"\n";
}

