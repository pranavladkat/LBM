#include <iostream>
#include <vector>

using namespace std;

const int IMAX = 100;
const int JMAX = IMAX;

const double tau = 2./3.;
const double Re = 100;

const double u_lid = Re*(tau-0.5)/(3.*(IMAX-1.));


int main()
{
    // directional weights for f_eq (for D2Q9)
    const vector<double> weight = {4./9.,1./9.,1./9.,1./9.,1./9.,1./36.,1./36.,1./36.,1./36.};

    // particle velocity
    const vector<vector<double>> e = {{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};

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


    cout << "Hello World!" << endl;
    return 0;
}
