#include <vector>
#include <stdio.h>
#include <stdlib.h>

using std::cout;

//Field class is an LxL grid with integer values.

template <class T>
class Field
{
  private:
    
    int              nx_;
    int              ny_;
    double         beta_;
    T            action_;
    T                 q_;
    std::vector<T> data_;

    //Function which handles periodic indexing on the grid
    int index_(int x, int y) const { return (x+nx_) % nx_ + nx_*((y+ny_) % ny_); }
    
    
    //Applies one update to the Potts-Model
    void update_ij_(int i, int j){

        // Random spin initialised as the current spin
        T random_q = data_[index_(i,j)];
        //Ensure random spin != current spin
        while(random_q==data_[index_(i,j)]) random_q = 1 + drand48() * q_;
        
        T proposed_action = action_ij(              random_q,
                                      data_[index_(i+1,j  )],
                                      data_[index_(i-1,j  )],
                                      data_[index_(i  ,j+1)],
                                      data_[index_(i  ,j-1)] );
        
        T  current_action = action_ij(data_[index_(i  ,j  )],
                                      data_[index_(i+1,j  )],
                                      data_[index_(i-1,j  )],
                                      data_[index_(i  ,j+1)],
                                      data_[index_(i  ,j-1)] );
        
        double dS = proposed_action - current_action; // Change in action
        
        //We accept in if the acceptance probablility > 1
        if(dS < 0 ) data_[ index_(i, j) ] = random_q, action_ += dS;
        
        //Accept with prob exp(-b*dS)
        else if( drand48() < exp(-beta_ * dS)) data_[ index_(i, j) ] = random_q, action_ += dS;
    }
    
    
    //Returns the count of the most frequently-occuring spin
    double f_() {
        
        std::vector<int> count(q_); // Keep track of spin occurrences
        for(int i=0; i<q_; i++) count[i] = 0; //initialise to 0
        
        //For each grid point check what spin (i) it is and add it to 'count[i]'
        for(int j = 0; j < ny_; j++){ for(int i = 0; i < nx_; i++){
                for(int test_q = 1; test_q <= q_; test_q++){
                    if(test_q == data_[index_(i,j)] ){ count[test_q-1]++;}
                }
            }
        }
        
        //Return the max value of the vector
        return *std::max_element(count.begin(), count.end())/(nx_*ny_+0.);
    }
    
    
    
  public:
    
    //Randomly initialising field
    Field(int nx, int ny, int q, double beta) :
    nx_(nx), ny_(ny), data_(nx_*ny_), q_(q), beta_(beta)
    
    {
        for (int i=0;i<nx*ny;i++) data_[i] = 1 + drand48() * q_; //Random initialisation in [1,..,q]

        //Loop over entire grid, and calculate initial action
        for(int j =0; j<ny_; j++) for(int i =0; i<nx_; i++){
                action_+=action_ij(data_[index_(i,  j  )],
                                   data_[index_(i+1,j  )],
                                   data_[index_(i-1,j  )],
                                   data_[index_(i  ,j+1)],
                                   data_[index_(i  ,j-1)]); }
    }
   
    //Calculates the action of phi(i,j); arguements are phi(i,j) followed by its NN's
    T action_ij(T p_ij,T p_ip1j,T p_im1j,T p_ijp1,T p_ijm1){
        T S = 0;
        
        if(p_ij!=p_ip1j)S++;if(p_ij!=p_ijp1)S++;if(p_ij!=p_im1j)S++;if(p_ij!=p_ijm1)S++;
        return S;
    }
    
    
    //Creates an operator which calls the field like phi(a,b)
    T  operator() (int x, int y) const { return data_[index_(x,y)]; }
    T& operator() (int x, int y)       { return data_[index_(x,y)]; }

    
    //Return Field size nx , ny, action
    int     nx() const { return     nx_ ; }
    int     ny() const { return     ny_ ; }
    int action() const { return action_ ; }

    
    //fractional magnetisation of a spin configuration
    double frac_mag() { return (q_*f_() - 1. )/(q_ - 1.); }
    
    // Itterate over its
    void itterate(int its){
        
        //Loop over its and every site i,j
        for(int k =0; k<its; k++) for(int j =0; j<ny_; j++) for(int i =0; i<nx_; i++) update_ij_(i, j);
    }
        

    
    //=========================================================//
    //                 Plotting and animation                  //
    //=========================================================//

    //Writing a CSV file to plot a heatmap
    void plot_heatmap(int itterations) {
        
        //Destination filename
        std::ofstream outfile ("heatmap.csv");
        //Writes headings
        outfile << "i" <<","<< "j"<<","<< "f" << std::endl;
        for(int k =0; k<itterations; k++){
            for(int j =0; j<ny_; j++)for(int i =0; i<nx_; i++){
                    //Writing to the CSV file
                    if(i==nx_-1 && j==ny_-1) outfile << data_[index_(i,j)] << std::flush;
                    else              outfile << data_[index_(i,j)] << "," << std::flush;
                    update_ij_(i, j);
                }

            outfile << std::endl;
        }
        outfile.close();
        
    }
    
    //Printing grid
    void print_grid() const {
        
        cout<<"\n\n";
        for(int j =0; j<ny_; j++){
            for(int i =0; i<nx_; i++) cout<<data_[index_(i,j)]<<" "<<std::flush;
            cout<<"\n";
        }
        cout<<"\n\n";
    }
};


