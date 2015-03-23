#include <iostream>
#include <vector>
#include <cassert>

using namespace std;

const int IMAX = 10;
const int JMAX = IMAX;

const double tau = 2./3.;
const double Re = 100;

const double u_lid = Re*(tau-0.5)/(3.*(IMAX-1.));


// directional weights for f_eq (for D2Q9)
const vector<double> w = {4./9.,1./9.,1./9.,1./9.,1./9.,1./36.,1./36.,1./36.,1./36.};
// particle velocity
const vector<vector<double>> e = {{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};



// initialize f, u
void initialize_variables(vector<vector<vector<double>>>& f,
                          vector<vector<double>>& rho,
                          vector<vector<double>>& u);

//bonunday condition functions
void BottomBC (vector<double>& f);
void LeftBC   (vector<double>& f);
void RightBC  (vector<double>& f);
void TopBC    (vector<double>& f);


int main()
{
    // define variables
    vector<vector<vector<double>>> f       (IMAX,vector<vector<double>>(JMAX,vector<double>(9,0)));
    vector<vector<vector<double>>> f_str   (IMAX,vector<vector<double>>(JMAX,vector<double>(9,0)));
    vector<vector<vector<double>>> f_eq    (IMAX,vector<vector<double>>(JMAX,vector<double>(9,0)));
    vector<vector<double>>         rho     (IMAX,vector<double>(JMAX,0));
    vector<vector<double>>         u       (IMAX,vector<double>(JMAX,0));
    vector<vector<double>>         v       (IMAX,vector<double>(JMAX,0));
    vector<vector<double>>         rho_old (IMAX,vector<double>(JMAX,0));
    vector<vector<double>>         u_old   (IMAX,vector<double>(JMAX,0));
    vector<vector<double>>         v_old   (IMAX,vector<double>(JMAX,0));


    // print setup values
    cout << "Reynolds number = " << Re << endl;
    cout << "Lid Velocity = " << u_lid << endl;
    cout << "Tau = " << tau << endl;
    cout << "Grid size = " << IMAX << " x " << JMAX << endl;

    // initialize
    initialize_variables(f,rho,u);



    cout << "Hello World!" << endl;
    return 0;
}



void initialize_variables(vector<vector<vector<double>>>& f,
                          vector<vector<double>>& rho,
                          vector<vector<double>>& u)
{
    // initialize f
    for(size_t i = 0; i < f.size(); i++)
        for(size_t j = 0; j < f[i].size(); j++)
            for(size_t k = 0; k < 9; k++)
                f[i][j][k] = w[k];

    // initialize rho, u
    for(size_t i = 0; i < u.size(); i++){
        for(size_t j = 0; j < u[i].size(); j++){
            rho[i][j] = 1.0;
            // set lid velocity
            if(j == u[i].size() - 1)
                u[i][j] = u_lid;
        }
    }
}


void BottomBC (vector<double>& f){

    assert(f.size() == 9);
    f[2] = f[4];
    f[5] = f[7];
    f[6] = f[8];
}


void LeftBC   (vector<double>& f){

    assert(f.size() == 9);
    f[1] = f[3];
    f[5] = f[7];
    f[8] = f[6];
}

void RightBC  (vector<double>& f){

    assert(f.size() == 9);
    f[3] = f[1];
    f[6] = f[8];
    f[7] = f[5];
}

void TopBC(vector<double> &f){

    assert(f.size() == 9);

    double rho = f[0] + f[1] + f[3] + 2 * (f[2] + f[5] + f[6]);

    f[4] = f[2];
    f[7] = 0.5*(f[1] - f[3] + 2*f[5] - rho*u_lid);
    f[8] = f[6] + 0.5*(f[3] - f[1] + rho*u_lid);
}
