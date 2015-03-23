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
void apply_boundary_conditions(vector<vector<vector<double>>>& f);
void BottomBC (vector<double>& f);
void LeftBC   (vector<double>& f);
void RightBC  (vector<double>& f);
void TopBC    (vector<double>& f);

// streaming function
void streaming_step(vector<vector<vector<double>>>& f, vector<vector<vector<double>>>& f_str);

// calculate rho, u, v
void compute_macroscopic_variables(const vector<vector<vector<double>>>& f,
                                   vector<vector<double>>& rho,
                                   vector<vector<double>>& u,
                                   vector<vector<double>>& v);

// compute equilibrium distribution function
void compute_f_equilibrium(const vector<vector<vector<double>>>& f,
                           const vector<vector<double>> rho,
                           const vector<vector<double>> u,
                           const vector<vector<double>> v,
                           vector<vector<vector<double>>>& f_eq);

void collision(const vector<vector<vector<double>>>& f_str,
               const vector<vector<vector<double>>>& f_eq,
               vector<vector<vector<double>>>& f);


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

    // stream
    streaming_step(f,f_str);

    compute_macroscopic_variables(f_str,rho,u,v);

    compute_f_equilibrium(f,rho,u,v,f_eq);

    collision(f_str,f_eq,f);

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

    apply_boundary_conditions(f);
}


void apply_boundary_conditions(vector<vector<vector<double>>>& f){

    //apply bottom bc
    for(size_t i = 0; i < f.size(); i++){
        vector<double>& boundary_node = f[i][0];
        BottomBC(boundary_node);
    }

    //apply left bc
    for(size_t i = 0; i <f[0].size(); i++){
        vector<double>& boundary_node = f[0][i];
        LeftBC(boundary_node);
    }

    //apply right bc
    for(size_t i = 0; i < f[0].size(); i++){
        vector<double>& boundary_node = f[f.size()-1][i];
        RightBC(boundary_node);
    }

    //apply top bc
    for(size_t i = 0; i < f.size(); i++){
        vector<double>& boundary_node = f[i][f[0].size()-1];
        TopBC(boundary_node);
    }

//    for(size_t i = 0; i < f.size(); i++){
//        for(size_t j = 0; j < f[i].size(); j++){
//            for(size_t k = 0; k < 9; k++)
//            cout << f[i][j][k] << "  ";
//            cout << endl;
//        }
//    }

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


void streaming_step(vector<vector<vector<double>>>& f, vector<vector<vector<double>>>& f_str){

    for(size_t i = 0; i < f.size(); i++){
        for(size_t j = 0; j < f[i].size(); j++){

            f_str[i][j][0] = f[i][j][0];

            if(i > 0)
                f_str[i][j][1] = f[i-1][j][1];

            if(j > 0)
                f_str[i][j][2] = f[i][j-1][2];

            if(i < f.size()-1)
                f_str[i][j][3] = f[i+1][j][3];

            if(j < f[i].size()-1)
                f_str[i][j][4] = f[i][j+1][4];

            if(i > 0 && j > 0)
                f_str[i][j][5] = f[i-1][j-1][5];

            if(i < f.size()-1 && j > 0)
                f_str[i][j][6] = f[i+1][j-1][6];

            if(i < f.size()-1 && j < f[i].size()-1)
                f_str[i][j][7] = f[i+1][j+1][7];

            if(i > 0 &&  j < f[i].size()-1)
                f_str[i][j][8] = f[i-1][j+1][8];
        }
    }

    apply_boundary_conditions(f_str);
}




void compute_macroscopic_variables(const vector<vector<vector<double>>>& f,
                                   vector<vector<double>>& rho,
                                   vector<vector<double>>& u,
                                   vector<vector<double>>& v)
{
    for(size_t i = 0; i < f.size(); i++){
        for(size_t j = 0; j < f[i].size(); j++){

            rho[i][j] = 0;
            u  [i][j] = 0;
            v  [i][j] = 0;

            for(size_t k = 0; k < 9; k++){
                rho[i][j] += f[i][j][k];
                u  [i][j] += f[i][j][k] * e[k][0];
                v  [i][j] += f[i][j][k] * e[k][1];
            }

            u[i][j] /= rho[i][j];
            v[i][j] /= rho[i][j];
        }
    }
}



void compute_f_equilibrium(const vector<vector<vector<double>>>& f,
                           const vector<vector<double>> rho,
                           const vector<vector<double>> u,
                           const vector<vector<double>> v,
                           vector<vector<vector<double>>>& f_eq)
{

    for(size_t i = 0; i < f.size(); i++){
        for(size_t j = 0; j < f[i].size(); j++){

            const double& u_eq = u[i][j];
            const double& v_eq = v[i][j];

            double U_square = u_eq*u_eq + v_eq*v_eq;

            f_eq[i][j][0] = w[0] * rho[i][j] * ( 1.0 - 1.5*U_square );

            for(size_t k = 1; k < 9; k++){

                double e_dot_u        = e[k][0] * u_eq + e[k][1] * v_eq;
                double e_dot_u_square = e_dot_u * e_dot_u;

                f_eq[i][j][k] = w[k] * rho[i][j] * (1 + 3*e_dot_u + 4.5*e_dot_u_square - 1.5*U_square);

            }

        }
    }

}



void collision(const vector<vector<vector<double>>>& f_str,
               const vector<vector<vector<double>>>& f_eq,
               vector<vector<vector<double>>>& f)
{

    for(size_t i = 0; i < f.size(); i++){
        for(size_t j = 0; j < f[i].size(); j++){
            for(size_t k = 0; k < 9; k++){

                f[i][j][k] = f_str[i][j][k] - (f_str[i][j][k] - f_eq[i][j][k]) / tau;
            }
        }
    }

}

