
#include <math.h>
#include "mex.h"
#include "matrix.h"



/* Define all constants required for simulation */


const double A[2] = {27.671988	,60.30274};
const double B[2] = {51.14898,1.034566};
const double C[2] = {-30.64454,	-0.192997};
const double D[2] = {6.847911,	0.012540};
const double E[2] = {-0.157906,	-6.860254};

const double AN[3] = {	28.98641,	19.50583,	35.51872};
const double BN[3] = {	1.853978,	19.88705,	1.128728};
const double CN[3] = {	-9.647459,	-8.598535,	-0.196103};
const double DN[3] = {	16.63537,	1.369784,	0.014662};
const double EN[3] = {0.000117,	0.527601,	-4.553760};
const double FN[3] = {-8.671914,	-4.935202,	-18.97091};
const double GN[3] = {226.4168,	212.3900,	224.9810};

const double AO[3] = {	31.32234,	30.03235,	20.91111};
const double BO[3] = {	-20.23531,	8.772972,	10.72071};
const double CO[3] = {	57.86644,	-3.988133,	-2.020498};
const double DO[3] = {	-36.50624,	0.788313,	0.146449};
const double EO[3] = {  -0.007374,	-0.741599,	9.245722};
const double FO[3] = {  -8.903471,	-11.32468,	5.337651}; 
const double GO[3] = {	246.7945,	236.1663,	237.6185};





static void cpf(
        
          
        double      x_nodes,
        double    MW_prod,
        double    MW, 
        double    rho_N2O,
		int	index[],
 		int	indexO[],
        int    indexN[],
        double    Mol_density_prod[],
        double    Mol_density_N2O[],
        double    Mol_density_Total[],
        double     X_N2O[],
        double    X_Prod[],
        double       gas[],
        double      dens[],
        double        cp[]
		   )
{
    int i;
           for (i = 0; i < x_nodes; i++){
               
               Mol_density_prod[i]   =  dens[i]/(MW_prod);
               Mol_density_N2O[i]    = ((rho_N2O)-dens[i])/(MW); 
               Mol_density_Total[i]  = Mol_density_prod[i] + Mol_density_N2O[i];
               
               X_N2O[i]  = Mol_density_N2O[i]/Mol_density_Total[i];
               X_Prod[i] = Mol_density_prod[i]/Mol_density_Total[i];
               
               // Temperature based restrictions
               
               if  (gas[i] < 1400){
                   index[i] = 0;
               }
               else {
                   index[i] = 1;
               }
               
               
               
               if  (gas[i] < 500){
                   indexN[i] = 0;
               }
               else if (gas[i] > 500 && gas[i] < 2000) {
                   indexN[i] = 1; 
               }
               else {
                   indexN[i] = 2; 
               }
               
               if  (gas[i] < 700){
                   indexO[i] = 0;
               }
               else if (gas[i] > 700 && gas[i] < 2000) {
                   indexO[i] = 1; 
               }
               else {
                   indexO[i] = 2; 
               }
               
               double gas2 = (gas[i]/1000.0f);
               double gas3 = gas2*gas2; 
               double gas4 = gas2*gas2*gas2; 
               double cp_prod = ((2.0f/3.0f)*AN[indexN[i]] + (1.0f/3.0f)*(AO[indexO[i]])) + gas2*
                        ((2.0f/3.0f)*BN[indexN[i]] + (1.0f/3.0f)*BO[indexO[i]]) + gas3*((2.0f/3.0f)*
                       CN[indexN[i]] + (1.0f/3.0f)*CO[indexO[i]]) + gas4*((2.0f/3.0f)*DN[indexN[i]] + 
                       (1.0f/3.0f)*DO[indexO[i]]) + (1/gas3)*((2.0f/3.0f)*EN[indexN[i]] + (1.0f/3.0f)*
                       EO[indexO[i]]);
               
               double Average_CP = X_N2O[i]*(A[index[i]] + B[index[i]]*gas2 + C[index[i]]*gas3 + D[index[i]]*gas4 + 
                      (E[index[i]]/gas3)) + X_Prod[i]*cp_prod; 
               
               double Average_MW = X_N2O[i] * MW + X_Prod[i]*MW_prod; 
              // printf("%d\n", indexN[i]);

                  
               
               cp[i] = Average_CP/Average_MW; 
                 
           }
 

}


static void main_iter(

        double      t_steps,        
        double      x_nodes,
        int      plot_step,
        int      count[], 
        int      limit, 
        double    Ru,
        double    cp_b, 
        double  rho_b,
        double      kb,
        double    Q_decomp, 
        double    MW_prod,
        double    MW, 
        double    rho_N2O,
        double       fac,
        double       fac2, 
        double          k,
        double         dx,
        double         dt, 
        double         m_t,
        double       Area, 
        double          h,
        double          u,
        double      sigma, 
		int	     index[],
 		int	     indexO[],
        int      indexN[],
        double    Mol_density_prod[],
        double    Mol_density_N2O[],
        double    Mol_density_Total[],
        double    X_N2O[],
        double    X_Prod[],
        double       gas[],
        double     block[],
        double      dens[],
        double     gas_p[],
        double   block_p[],
        double    dens_p[],
        double        cp[],
        double      beta[],
        double   del_2_term[],
        double   hg_b[],
        double   convec[],
        double   decomp[],
        double   del_2_termb[],
        double  hb_g[], 
        double del_rho[],
        double rho_change[],
        double gas_pl[],
        double block_pl[],
        double dens_pl[],
        int    num_vec,
        double time_to_change, 
        double u_vec[],
        double main_u_vec[], 
        int   var_u
        
        
		   )
{
    
    int n;
    int i; 
    int ind;
    int u_count = 1; 
    int time_before_u = 0; 
    bool stop        = false; 
    
     int found_it   = 0; // Location of decomposition 
     int speed_iter = 0; // # of times gas is solved within dt
   
    
          for (n = 1; n < t_steps && !stop; n++){

              
               // Progess of Simulation
               
              double ratio_t = ((double)(n)/(double)t_steps); 
    
               
               if (n % plot_step == 0) {
                   printf("\n");
                   printf("%f of simulation has been completed\n",ratio_t); 
               }
                                  
          
               // Changing velocity 
              
              if (found_it == 1 && u_count < num_vec && var_u == 1){
                time_before_u = time_before_u +1;
                    if (time_before_u  == time_to_change){
                        printf("\n");
                        printf("On %d of %d iterations", u_count, num_vec);
                        printf("\n");
                        
                        time_before_u = 0; 
                        u = u_vec[u_count];
                        u_count++;
                    }
                
                  
                  
              }
              
               if (speed_iter == limit){
                     speed_iter = 0;
               }

                cpf(x_nodes, MW_prod,MW, rho_N2O, index, indexO, indexN, Mol_density_prod, Mol_density_N2O, 
                 Mol_density_Total, X_N2O, X_Prod, gas, dens, cp);
       
                
                 while(speed_iter<limit){
                   
                    
                     
                        for (i = 0; i < x_nodes; i++) {
                                               
                            if (gas[i] != gas[i]){ // Checking for NaN in the middle, random
                                printf("NaN's have been found within the gas solution\n");
                                stop = true; 
                            
                            }
                            
                            beta[i] = 1.3e11*exp(-59200/(Ru*gas_p[i]));
                           // printf("%f", beta[i]);
                            
                            if (i == 0) {
                                if (n > 1){
                                    gas[i]   = (10.0*block_p[i] + (k/dx)*300.0 + (m_t/Area)*cp[i]*300.0 + 
                                            (k/dx)*gas_p[i+1])/((k/dx) + (m_t/Area)*cp[i] + (k/dx) + 10.0); 
                        
                                    
                                }
                       
                         
                            }
                            
                            else  if  (i == x_nodes-1){
                               if (n>1) {
                                    gas[i]   = (10.0*block_p[i] + (k/dx)*(gas_p[i-1]) + 
                                            (m_t/Area)*cp[i]*(gas_p[i-1]))/((k/dx) + (m_t/Area)*cp[i] + 10.0); 
                                  
                                } 
                            
                            }
                            
                            else {
                  
                            
                               del_2_term[i] =  k*(gas_p[i+1] - 2.0f*gas_p[i] + gas_p[i-1])/(dx*dx);
                              
                                hg_b[i]       = h*sigma*(block_p[i] - gas_p[i]);
                                convec[i]     = u*cp[i]* rho_N2O * (gas_p[i+1] - gas_p[i-1])/(2.0*dx);
                                decomp[i]     = beta[i] *Q_decomp*(rho_N2O - dens_p[i]);
                                gas[i]        = fac*(dt/(rho_N2O*cp[i]))*(del_2_term[i] + hg_b[i] - convec[i]
                                        + decomp[i]) + gas_p[i]; 
                                 
                                }
                            
                        }
                        speed_iter++;
                     // printf("%d\n", speed_iter);
                       // printf("%f\n", cp[x_nodes-1]); 
                     
                        // Will need to sum over add the variables into gas_p before while loop
                     for (i = 0; i<x_nodes; i++){
                      gas_p[i]    = gas[i];
                     }
             
                 }
                       
          
                 for (i = 0; i < x_nodes; i++) {
                     
                                          
                            if (block[i] != block[i]){ // Checking for NaN in the middle, random
                                printf("NaN's have been found within the block solution \n");
                                stop = true; 
                            
                            }

                            if (i == 0) {
                                if (n > 1){
                                    block[i]   = ((20.0*300.0) + kb*block_p[i+1]/dx) * ( 1/((kb/dx) + 20.0)); 

                                }
                       
                         
                            }
                            
                            else  if  (i == x_nodes-1){
                               if (n>1) {
                                    block[i]   = (20.0*gas_p[i] + kb*block_p[i-1]/dx) * (1/((kb/dx) + 20.0));
                               
                            
                                } 
                            
                            }
                            
                            else {
                                 del_2_termb[i] = kb*(block_p[i+1] - 2.0*block_p[i] + block_p[i-1])/(dx*dx);
                                 hb_g[i]        = h*sigma*(gas_p[i] - block_p[i]);
                                block[i]        = (fac2*dt/(rho_b*cp_b))*(del_2_termb[i] + hb_g[i]) +  block_p[i];
                              
                            
                                    
                                }
                 }
                
                         
                 for (i = 0; i < x_nodes; i++) {
                     
                                          
                            if (dens[i] != dens[i]){ // Checking for NaN in the middle, random
                                printf("NaN's have been found within the dens solution \n");
                                stop = true; 
                            
                            }

                            if (i == 0) {
                                    dens[i]   = 0;
                            }
                            
                            else {
                                   rho_change[i] = beta[i]*(rho_N2O - dens_p[i]); 
                                    del_rho[i]    = u*(dens_p[i] - dens_p[i-1])/dx; 
                                    dens[i]       = dt*(-del_rho[i] + rho_change[i]) + dens_p[i]; 
                                } 
                 }
                            
 
                 for (i = 0; i<x_nodes; i++){
                      block_p[i] = block[i];
                      dens_p[i]  = dens[i];
                 }
            
                      if (n % plot_step == 0) {
                          main_u_vec[count[0]]  = u;
                         for (i = 0; i<x_nodes; i++){
                             ind = (count[0]*x_nodes) + i; 
                            gas_pl[ind]    = gas[i];
                            block_pl[ind]  = block[i];
                            dens_pl[ind]   = dens[i];
                         }
                         count[0]++;
                      
                      
               }
            
             if (found_it < 1){
                 int index =  (int)round((x_nodes*0.05/10.0)*10);
                if (dens[index]>0.99*rho_N2O == 1){
                    printf("\n");
                    printf("Decomposition has been achieved within this system (at ~1/20 of system length) for this iteration"); 
                    printf("\n");
                    found_it = 1;
                }
                 
                 
             }
                
         } 
}


void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{
    
     // RUN THIS WHEN TRYING TO RUN THE STUFF 
   // MEX_MATLAB(single(rho_N2O), single(u), single(k), single(sigma), single(h), single(kb), single(rho_b), single(cp_b), 
    //single(L), single(dx), single(cp_b), int32(x_nodes), single(t_steps), single(x_mesh), single(fac), single(limit), single(fac2), single(gas_pl), single(block_pl), single(dens_pl), single(gas), single(block), single(dens), single(cp), int32(0),int32(0),int32(0),int32(0))

   
    /* loading all values needed for this ridiculously long ass code */ 
    
    int i = 0;
    
    // Gas Constants loaded in
    double rho_N2O          = (double)  mxGetScalar(prhs[i]); i++;
    double u                = (double)  mxGetScalar(prhs[i]); i++;
    double k                = (double)  mxGetScalar(prhs[i]); i++;
    
    // Material Properties 
    double sigma            = (double)  mxGetScalar(prhs[i]); i++;
    double h                = (double)  mxGetScalar(prhs[i]); i++;
    double kb               = (double)  mxGetScalar(prhs[i]); i++;
    double rho_b            = (double)  mxGetScalar(prhs[i]); i++;
    double cp_b             = (double)  mxGetScalar(prhs[i]); i++;
    double *cp              = (double*) mxGetData  (prhs[i]); i++;
    
    // System Parameters
    double L                  = (double)  mxGetScalar(prhs[i]); i++;
    double dx                 = (double) mxGetScalar(prhs[i]); i++;
    double dt                 = (double) mxGetScalar(prhs[i]); i++;
    double x_nodes            = (double)    mxGetScalar(prhs[i]); i++;
    double t_steps            = (double)    mxGetScalar(prhs[i]); i++;
    double *x_mesh            = (double*) mxGetData  (prhs[i]); i++;
    int    plot_step          =  (int)   mxGetScalar(prhs[i]); i++;
    int  *count               =  (int*)   mxGetData(prhs[i]); i++;
             
    // Speed of Simulation Parameters
    double fac             = (double)  mxGetScalar(prhs[i]); i++;
    double limit           = (double)  mxGetScalar(prhs[i]); i++;
    double fac2            = (double)  mxGetScalar(prhs[i]); i++;
    
            
    /* Pointers for all arrays that need to change */
    // Pointers to arrays for plotting        
    double  *gas_pl         = (double*) mxGetData  (prhs[i]); i++;
    double  *block_pl       = (double*) mxGetData  (prhs[i]); i++;
    double  *dens_pl        = (double*) mxGetData  (prhs[i]); i++;
    
    // Pointers for all arrays where each iterating-calculation occurs
    double  *gas            = (double*) mxGetData  (prhs[i]); i++;
    double  *block          = (double*) mxGetData  (prhs[i]); i++;
    double  *dens           = (double*) mxGetData  (prhs[i]); i++;   
      
   
    
    // Loading in parameters if alternating speed 
    int var_u               = (int)    mxGetScalar(prhs[i]); i++;
    double min_u            = (double)  mxGetScalar(prhs[i]); i++;
    double max_u            = (double)  mxGetScalar(prhs[i]); i++;
    int time_to_change      = (int)    mxGetScalar (prhs[i]); i++;
    double *u_vec            = (double*)  mxGetData(prhs[i]); i++;
    double *main_u_vec      = (double*) mxGetData (prhs[i]); i++; 
    int num_vec          = (int)  mxGetScalar(prhs[i]); i++; 
    
    // Initializing arrays for storing previous time-step solutions 
    double *gas_p;
    gas_p   = (double*) malloc(x_nodes*sizeof(double));
    double *block_p;
    block_p = (double*) malloc(x_nodes*sizeof(double));
    double *dens_p;
    dens_p  = (double*) malloc(x_nodes*sizeof(double));
    
    // Initializing arrays for storing each component of gas heat transfer
    double *del_2_term;
    del_2_term   = (double*) malloc(x_nodes*sizeof(double));
    double *hg_b;
    hg_b = (double*) malloc(x_nodes*sizeof(double));
    double *convec;
    convec  = (double*) malloc(x_nodes*sizeof(double));
    double *decomp;
    decomp   = (double*) malloc(x_nodes*sizeof(double));
    double *beta;
    beta = (double*) malloc(x_nodes*sizeof(double));
   
    
    // Initializing arrays for storing each component of block heat transfer
    double *del_2_termb;
    del_2_termb   = (double*) malloc(x_nodes*sizeof(double));
    double *hb_g;
    hb_g = (double*) malloc(x_nodes*sizeof(double));
    
    // Initializing arrays for density terms 
    double *del_rho;
    del_rho   = (double*) malloc(x_nodes*sizeof(double));
    double *rho_change;
    rho_change = (double*) malloc(x_nodes*sizeof(double));

    // Initilizing arrays needed for calculating cp
    int *index;
    index   = (int*) malloc(x_nodes*sizeof(int));
    int *indexN;
    indexN = (int*) malloc(x_nodes*sizeof(int));
    int *indexO;
    indexO  = (int*) malloc(x_nodes*sizeof(int));

    double *Mol_density_prod;
    Mol_density_prod  = (double*) malloc(x_nodes*sizeof(double));
    double *Mol_density_N2O;
    Mol_density_N2O = (double*) malloc(x_nodes*sizeof(double));
    double *Mol_density_Total;
    Mol_density_Total  = (double*) malloc(x_nodes*sizeof(double));
    double * X_N2O;
     X_N2O  = (double*) malloc(x_nodes*sizeof(double));
    double *X_Prod;
    X_Prod  = (double*) malloc(x_nodes*sizeof(double));


    
    
        
    // Area constant doesn't really matter
    const double Area     = 0.0154; 
    
    // Calculating rest of the variables 
    const double MW       = 44.013f/1000.0f;
    const double Q_decomp = ((82.0f)*(1000.0f))/MW;
    const double m_t      =  u*rho_N2O*Area;
    const double MW_prod  = (2.0f*28.01340f + 31.999f)/(3.0f*1000.0f);

    
    /* Constants that aren't passed in */
    const double Ru       = 1.987; 

    
    // Time to copy all initial values to the previous matrixes. 
    
    int ind;
    for (ind = 0; ind<x_nodes; ind++){
        gas_p[ind]   = gas[ind]; 
      //  printf("%f\n", gas[ind]);
        block_p[ind] = block[ind];
        dens_p[ind]  = dens[ind]; 
    }
    
  
 //printf("%d\n",t_steps);
 //printf("%f\n",x_mesh[x_nodes-1]);
// printf("%f\n",gas[200]);
 //printf("%f\n", cp[10]);
 //printf("%lf\n",dx);

  

        main_iter(t_steps, x_nodes, plot_step, count, limit, Ru,cp_b, rho_b,kb, Q_decomp, MW_prod, MW, rho_N2O, fac,fac2,  k, dx, dt, m_t, Area,
        h, u, sigma, index, indexO, indexN, Mol_density_prod, Mol_density_N2O, Mol_density_Total, X_N2O, X_Prod, gas, 
        block, dens, gas_p, block_p, dens_p, cp, beta, del_2_term, hg_b, convec, decomp, del_2_termb, hb_g,del_rho, rho_change, gas_pl, block_pl, dens_pl,num_vec,
        time_to_change, u_vec, main_u_vec,var_u);

    
    
   // Will need to code in variables for changing speed 
    
    
    
}
