#include <iostream>
#include <fstream>
using std::cout;
#include "Header.h"


int main(){
    
    int           L =  20;
    int           q =   3;
    int itterations = 2000, av_over = 500;
    
    int divisions   = 500;
    double beta_0   = 0.5 , beta_f   = 1.5;
    
    //The division that we itterate over for graphing
    double beta_itts = (beta_f - beta_0)/divisions;

    //Writing CSV file to plot
    std::ofstream outfile ("beta_vs_M.csv");
    outfile << "beta" <<","<< "<M>"<< std::endl;

    //Recording data for plot
    for(double beta = beta_0; beta<=beta_f; beta+=beta_itts){

        //Initialise field
        Field<int> phi(L,L,q,beta);
        //Itterate the field
        phi.itterate(itterations-av_over);

        // Itterate field and record magnetism over 'av_over' itterations
        double mag = 0;
        for(int i = 0; i<av_over;i++){
            mag+=phi.frac_mag();
            phi.itterate(1);
        }
        mag = mag/(av_over+0.);

        // Write values
        outfile << beta << "," << mag << std::endl;
    }
    
    outfile.close();
    return 0;
    
}
