#include <iostream>
#include <string>
#include <fstream>
using std::cout;
#include <iomanip>
#include <vector>
#include "Field.h"

int nx          = 100;
int ny          = nx;
// These indices keep track of what parts of the grid we leave untouched
std::vector<int> boundary_indices_x(ny*nx), boundary_indices_y(ny*nx);

//This function returns the y-derivative according to central difference approximation
double y_derivative(Field p, int i, int j){
    return ny * (p(i,j+1) - p(i,j-1)) / 2;
}

//This function initialises the field according to the BS's
auto field_initialise(){
    Field phi( nx, ny );// Creating field instance
    int k = 0 ;         // Index which keeps track of boundary conditions

    for( int j = 1; j < ny; j++ ){
        for( int i = 1; i < nx; i++ ){
            
            //turns indices to spatial representation on the grid
            double i_n = (i+0.)/(nx+0.); double j_n = (j+.0)/(ny+0.);
            
            //initialising boundary indices
            boundary_indices_x[k] = -1; boundary_indices_y[k] = -1;
            
            //L shape BC's phi=-1
            if( (i_n >= 0.5) && (i_n <= 0.8) && (j_n == 0.7) ){
                phi(i,j) = -1.;
                boundary_indices_x[k] = i; boundary_indices_y[k] = j;
            } else if ( (j_n >= 0.3) && (j_n <= 0.7) && (i_n == 0.8) ){
                phi(i,j) = -1.;
                boundary_indices_x[k] = i; boundary_indices_y[k] = j;
            }
            
            //Square BC's
            if( (i_n >= 0.2) && (i_n <= 0.4)  && (j_n >= 0.1) && (j_n <= 0.3)  ){
                phi(i,j) = 2.;
                boundary_indices_x[k] = i; boundary_indices_y[k] = j;
            }
            k++;
        }
    }
    return phi;
}

//Writing a CSV file to plot a heatmap
void plot_heatmap(Field phi){
    
    //Destination filename
    std::ofstream outfile ("heatmap.csv");
    //Writes headings
    outfile << "i" <<","<< "j"<<","<< "f" << std::endl;
    for (int j = 0; j<=ny; j++){
        for (int i = 0; i<=nx; i++){
            //Writing to the CSV file
            outfile << i <<","<< j<<","<< phi(i,j) << std::endl;
        }
        cout<<std::endl;
    }
    outfile.close();
}

double field_solver(Field phi, double omega, int itterations){
    
    for( int count=0; count<itterations; count++){
        int    k = 0; //Keep track of boundary indices
        double a = 0;
        for( int j = 1; j < ny; j++ ){
            for( int i = 1; i < nx; i++ ){
                
                //Solve the laplace eqn. using SOR method
                if( i != boundary_indices_x[k] && j!= boundary_indices_y[k] ) {
                    a = phi(i,j);
                    phi(i,j) = omega/4.*( phi(i+1,j) + phi(i-1,j) + phi(i,j-1) + phi(i,j+1)) + (1 - omega) * a;
                }
                k++; //Count for boundary indices
            }
        }
    }
    
    return y_derivative(phi, 0.2*nx, 0.7*ny);
}

int main(){
    
    int itterations = 100;
    double omega    = 0.; // SOR method constant


    //Finding best value of omega using that the system converges to phi'(0.2,0.7) = -1.08617
    
    double best_omega = 5.;
    double best_error = 5.;
    std::ofstream outfile ("error_vs_omega.csv");
    outfile << "omega" <<","<< "dphi" << std::endl;

    for(int i = 0; i<1999; i++){
        omega+=0.001;
        // Initialising the field class
        Field   phi = field_initialise();
        double test = field_solver(phi, omega, itterations);
        double error = abs( -1.08617 - test );
    

        if(error < best_error){
            best_omega = omega;
            best_error = error;
        }

        outfile << omega << ","<< error << std::endl;
        //cout<<"\n"<<"phi'(0.2,0.7)="<<field_solver(phi, omega, itterations)<<"\n";
        //cout<<"\n"<<"omega="<<omega<<"\n";


    }
    outfile.close();

    cout<<"\n"<<"best omega="<<best_omega<<"\n";

    Field   Phi = field_initialise();
    cout<<"\n"<<"phi'(0.2,0.7)="<<field_solver(Phi, best_omega, itterations = 500)<<"\n";

    return 0;
}
