#include <iostream>
#include <vector>
#include <string>
#include <fstream>

//defines the 1st order ODE to be solved {dx/dt := f}
double f( double xx, double tt){
    return -pow(xx,2)*pow(tt,3);
}

//implements an update according to Runge-Kutta
double rk_update( double t, double xi, double dt ){
    
    double k1,k2,k3,k4;
    double ti = t;
    
    k1 = dt*f(xi, ti);
    k2 = dt*f(xi + k1/2., ti + dt/2.);
    k3 = dt*f(xi + k2/2., ti + dt/2.);
    k4 = dt*f(xi + k3, ti + dt);
    
    return xi+ 1./6.* (k1 + 2*k2 + 2*k3 + k4);
}


int main()
{
    
    double x0 =  1. ,t0 = 0., tf = 4., dt = 0.0001;
    double t = t0, x = x0;

    std::ofstream outfile ("graph1.csv"); //Writing a CSV file for plotting
    outfile << "t" <<","<< "y" << std::endl; //Writes headings

    //loops over 1 RK update
    while(t <= tf){
        outfile << t <<","<< x << std::endl;
        x = rk_update(t, x0, dt); //Updates x according to RK
        x0 = x;
        t += dt;
    }
    
    outfile.close();
    std::cout << "x(t="<< t <<")="<< x<<"\n"; //Printing result

}
