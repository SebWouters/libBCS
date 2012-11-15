#include <iostream>
#include <cstdlib>
#include <cmath>
#include <time.h>

using std::cout;
using std::endl;

#include "BCS.h"

// Sebastian Wouters, November 14, 2012.

extern "C" { //both are ARPACK functions

   void dsaupd_(int *ido, char *bmat, int *n, char *which, int *nev, double *tol, double *resid, int *ncv, double *v, int *ldv, int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int*info);
   void dseupd_(bool *rvec,char *howmny,bool *select,double *d,double *Z,int *ldz,double *sigma,char *bmat,int *n,char *which,int *nev,double *tol,double *resid,int *ncv,double *v,int *ldv,int *iparam,int *ipntr,double *workd,double *workl, int *lworkl,int *info);

}

double BCS::Solve(){

   int MAGICNUMBER = 32;
   double tol = 0.0; // tolerance for the calculation   -->   0.0=machine precision

   int ido = 0;  // internal communication flag
   char bmat = 'I';  // standard eigenvalue problem : A*x = lambda*x
   int n = pointer[NN];  // dimension of the eigenvalue problem

   char * which = new char[2];  // which eigenvalues to compute : the nev smallest ones
   which[0] = 'S';
   which[1] = 'A';
   int nev = 1; // only the smallest one

   double * resid = new double[n]; //residual vector; ini guess = previous A
   for (int counter=0; counter < n; counter++){
      resid[counter] = ((double) rand())/RAND_MAX;
   }
   int info = 0; // we use a random initial guess for resid --> info=0 means random guess
   
   int ncv = (n>MAGICNUMBER)?MAGICNUMBER:n; //number of Arnoldi vectors; capped to be tractable
   double * v = new double [n*ncv]; //space for the arnoldi vectors
   int ldv = n; // leading dimension of the v-mx

   int * iparam = new int[11]; // some options passed to the program
   iparam[0] = 1; // use exact shifts
   iparam[2] = 1000; // maximum number of iterations
   iparam[6] = 1; // standard eigenvalue problem A*x = l*x

   int * ipntr = new int[11]; // pointer to memory
   double * workd = new double[3*n]; // work space for the package, not to be used by users during the optim.
   int lworkl = ncv*ncv + 8*ncv; // size of work space nr. 2
   double * workl = new double[lworkl]; // work space nr. 2

   dsaupd_(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);

   while (ido != 99) // ido 99 means convergence has been reached
   {
   
      HtimesVec(workd + ipntr[1] - 1, workd + ipntr[0] - 1);
      dsaupd_(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);

   }

   int rvec = 0; // compute eigenvectors
   bool * rvec2 = (bool *) &rvec; //fortran <-> C++ mismatch of bool size
   *rvec2 = true;
   char howmny = 'A'; // calculate all of the nev eigenvectors
   int select[ncv]; // workspace for reordering the eigenvalues
   bool * select2 = (bool *) select; //fortran <-> C++ mismatch of bool size

   int ldz = n; //leading dimension Z-array = n x nev array containing the ritz vectors --> store it in A
   double sigma; //unreferenced pm

   dseupd_(rvec2,&howmny,select2,&diagEnergy,diagresult,&ldz,&sigma,&bmat,&n,which,&nev,&tol,resid,&ncv,v,&ldv,iparam,ipntr,workd,workl,&lworkl,&info);

   delete [] which;
   delete [] resid;
   delete [] v;
   delete [] iparam;
   delete [] ipntr;
   delete [] workd;
   delete [] workl;

   if(info != 0)
   {
      cout << "dseupd has exited with an error, program aborted" << endl;
      cout << "check the dseupd manual for error code: " << info << endl;
      exit(666);
   }
   
   CalculateRDMs(); //Calculate the RDMs once the problem is solved.
   
   return diagEnergy;

}

void BCS::HtimesVec(double * result, double * vec){

   //As result is always accessed in the particle number sector defined by Nvals[cnt], the following loop can be parallellized
   #pragma omp parallel for schedule(dynamic) default(none) shared(result,vec)
   for (int cnt=0; cnt<NN; cnt++){
   
      int Nup = (Nvals[cnt] + TwoSz)/2;
      int Ndown = (Nvals[cnt] - TwoSz)/2;
      
      //Set result to the chemical potential part
      double prefact = - (mu_up * Nup + mu_down * Ndown);
      for (unsigned long cnt2=pointer[cnt]; cnt2<pointer[cnt+1]; cnt2++){
         result[cnt2] = prefact * vec[cnt2];
      }
      
      int * state_alpha = new int[L]; //State has length L and contains zeros and ones corresponding to resp. empty and full sites
      int * state_work2 = new int[L];
      for (unsigned int cnt_alpha=0; cnt_alpha<gBinomial(L,Nup); cnt_alpha++){
         unsigned int State_Alpha_IntForm = CntToState[Nup][cnt_alpha];
         ConvertType3(state_alpha,State_Alpha_IntForm);
         for (unsigned int cnt_beta=0; cnt_beta<gBinomial(L,Ndown); cnt_beta++){
            unsigned int State_Beta_IntForm = CntToState[Ndown][cnt_beta];
            
            unsigned long rslt_ptr = pointer[cnt]+cnt_alpha+gBinomial(L,Nup)*cnt_beta;
            prefact = 0.0;
            
            //Do a bitwise comparison to get the doubly occupied sites for the repulsion energy, which are added to prefact
            unsigned int State_Double = State_Alpha_IntForm & State_Beta_IntForm;
            ConvertType3(state_work2,State_Double);
            for (int cnt2=0; cnt2<L; cnt2++){
               if (state_work2[cnt2]==1){
                  prefact += Uelem[cnt2];
               }
            }

            //Kinetic energy [ALPHA el.]: Diag terms of the alpha electrons are added to prefact, off-diagonal terms are added immediately
            for (int cnt2=0; cnt2<L; cnt2++){
               if (state_alpha[cnt2]==1){
                  prefact += Telem[cnt2*(L+1)];
                  int Nswaps = 0;
                  for (int cnt3=cnt2+1; cnt3<L; cnt3++){
                     if (state_alpha[cnt3]==0){
                        state_alpha[cnt2] = 0;
                        state_alpha[cnt3] = 1;
                        unsigned int cnt_alpha_bis = StateToCnt[ConvertType2(state_alpha)];
                        result[rslt_ptr] += ((Nswaps%2==0)?1:-1) * Telem[cnt2 + L*cnt3] * vec[pointer[cnt]+cnt_alpha_bis+gBinomial(L,Nup)*cnt_beta];
                        state_alpha[cnt3] = 0;
                        state_alpha[cnt2] = 1;
                     } else {
                        Nswaps++;
                     }
                  }
                  
                  Nswaps = 0;
                  for (int cnt3=cnt2-1; cnt3>=0; cnt3--){
                     if (state_alpha[cnt3]==0){
                        state_alpha[cnt2] = 0;
                        state_alpha[cnt3] = 1;
                        unsigned int cnt_alpha_bis = StateToCnt[ConvertType2(state_alpha)];
                        result[rslt_ptr] += ((Nswaps%2==0)?1:-1) * Telem[cnt2 + L*cnt3] * vec[pointer[cnt]+cnt_alpha_bis+gBinomial(L,Nup)*cnt_beta];
                        state_alpha[cnt3] = 0;
                        state_alpha[cnt2] = 1;
                     } else {
                        Nswaps++;
                     }
                  }
               }
            }
            
            //Kinetic energy [BETA el.]: Diag terms of the beta electrons are added to prefact, off-diagonal terms are added immediately
            ConvertType3(state_work2,State_Beta_IntForm);
            for (int cnt2=0; cnt2<L; cnt2++){
               if (state_work2[cnt2]==1){
                  prefact += Telem[cnt2*(L+1)];
                  int Nswaps = 0;
                  for (int cnt3=cnt2+1; cnt3<L; cnt3++){
                     if (state_work2[cnt3]==0){
                        state_work2[cnt2] = 0;
                        state_work2[cnt3] = 1;
                        unsigned int cnt_beta_bis = StateToCnt[ConvertType2(state_work2)];
                        result[rslt_ptr] += ((Nswaps%2==0)?1:-1) * Telem[cnt2 + L*cnt3] * vec[pointer[cnt]+cnt_alpha+gBinomial(L,Nup)*cnt_beta_bis];
                        state_work2[cnt3] = 0;
                        state_work2[cnt2] = 1;
                     } else {
                        Nswaps++;
                     }
                  }
                  
                  Nswaps = 0;
                  for (int cnt3=cnt2-1; cnt3>=0; cnt3--){
                     if (state_work2[cnt3]==0){
                        state_work2[cnt2] = 0;
                        state_work2[cnt3] = 1;
                        unsigned int cnt_beta_bis = StateToCnt[ConvertType2(state_work2)];
                        result[rslt_ptr] += ((Nswaps%2==0)?1:-1) * Telem[cnt2 + L*cnt3] * vec[pointer[cnt]+cnt_alpha+gBinomial(L,Nup)*cnt_beta_bis];
                        state_work2[cnt3] = 0;
                        state_work2[cnt2] = 1;
                     } else {
                        Nswaps++;
                     }
                  }
               }
            }
            
            //Alpha & beta electrons diag. kinetic energy terms + on-site repulsion terms.
            result[rslt_ptr] += prefact * vec[rslt_ptr];
            
            //Do the particle non-conserving part: raising.
            if (cnt>0){ //So the result state can be created by adding two particles to a lower particle number state state.
               int Nswaps_alpha = 0;
               for (int cnt2=L-1; cnt2>=0; cnt2--){
                  if (state_alpha[cnt2]==1){
                     state_alpha[cnt2] = 0;
                     unsigned int cnt_alpha_bis = StateToCnt[ConvertType2(state_alpha)];
                     int Nswaps = Nswaps_alpha;
                     for (int cnt3=0; cnt3<L; cnt3++){
                        if (state_work2[cnt3]==1){
                           state_work2[cnt3] = 0;
                           unsigned int cnt_beta_bis = StateToCnt[ConvertType2(state_work2)];
                           result[rslt_ptr] += ((Nswaps%2==0)?1:-1) * Delta[cnt2+L*cnt3] * vec[pointer[cnt-1]+cnt_alpha_bis+gBinomial(L,Nup-1)*cnt_beta_bis];
                           state_work2[cnt3] = 1;
                           Nswaps++;
                        }
                     }
                     state_alpha[cnt2] = 1;
                     Nswaps_alpha++;
                  }
               }
            }
            
            //Do the particle non-conserving part: lowering.
            if (cnt<NN-1){ //So the result state can be created by removing two particles from a higher particle number state.
               int Nswaps_alpha = 0;
               for (int cnt2=L-1; cnt2>=0; cnt2--){
                  if (state_alpha[cnt2]==0){
                     state_alpha[cnt2] = 1;
                     unsigned int cnt_alpha_bis = StateToCnt[ConvertType2(state_alpha)];
                     int Nswaps = Nswaps_alpha;
                     for (int cnt3=0; cnt3<L; cnt3++){
                        if (state_work2[cnt3]==0){
                           state_work2[cnt3] = 1;
                           unsigned int cnt_beta_bis = StateToCnt[ConvertType2(state_work2)];
                           result[rslt_ptr] += ((Nswaps%2==0)?1:-1) * Delta[cnt2+L*cnt3] * vec[pointer[cnt+1]+cnt_alpha_bis+gBinomial(L,Nup+1)*cnt_beta_bis];
                           state_work2[cnt3] = 0;
                        } else {
                           Nswaps++;
                        }
                     }
                     state_alpha[cnt2] = 0;
                  } else {
                     Nswaps_alpha++;
                  }
               }
            }
            
         }
      }
      delete [] state_alpha;
      delete [] state_work2;
   }
   
}


//Calculate the RDMs once the problem is solved (called by solver)
void BCS::CalculateRDMs(){

   for (int cnt=0; cnt<L*L; cnt++){
      RDMup[cnt] = 0.0;
      RDMdown[cnt] = 0.0;
      AddTwoParticles[cnt] = 0.0;
   }
   
   for (int cnt=0; cnt<NN; cnt++){
   
      int Nup = (Nvals[cnt] + TwoSz)/2;
      int Ndown = (Nvals[cnt] - TwoSz)/2;
      
      int * state_alpha = new int[L]; //State has length L and contains zeros and ones corresponding to resp. empty and full sites
      int * state_beta = new int[L];
      for (unsigned int cnt_alpha=0; cnt_alpha<gBinomial(L,Nup); cnt_alpha++){
         unsigned int State_Alpha_IntForm = CntToState[Nup][cnt_alpha];
         ConvertType3(state_alpha,State_Alpha_IntForm);
         for (unsigned int cnt_beta=0; cnt_beta<gBinomial(L,Ndown); cnt_beta++){
            unsigned int State_Beta_IntForm = CntToState[Ndown][cnt_beta];
            ConvertType3(state_beta,State_Beta_IntForm);

            double bra_coeff = diagresult[pointer[cnt]+cnt_alpha+gBinomial(L,Nup)*cnt_beta];

            //RDMup
            for (int cnt2=0; cnt2<L; cnt2++){
               if (state_alpha[cnt2]==1){
                  RDMup[cnt2*(L+1)] += bra_coeff * bra_coeff; //diag
                  int Nswaps = 0;
                  for (int cnt3=cnt2+1; cnt3<L; cnt3++){ //upper diag
                     if (state_alpha[cnt3]==0){
                        state_alpha[cnt2] = 0;
                        state_alpha[cnt3] = 1;
                        unsigned int cnt_alpha_bis = StateToCnt[ConvertType2(state_alpha)];
                        RDMup[cnt2+L*cnt3] += ((Nswaps%2==0)?1:-1) * bra_coeff * diagresult[pointer[cnt]+cnt_alpha_bis+gBinomial(L,Nup)*cnt_beta];
                        state_alpha[cnt3] = 0;
                        state_alpha[cnt2] = 1;
                     } else {
                        Nswaps++;
                     }
                  }
               }
            }
            
            //RDMdown
            for (int cnt2=0; cnt2<L; cnt2++){
               if (state_beta[cnt2]==1){
                  RDMdown[cnt2*(L+1)] += bra_coeff * bra_coeff; //diag
                  int Nswaps = 0;
                  for (int cnt3=cnt2+1; cnt3<L; cnt3++){ //upper diag
                     if (state_beta[cnt3]==0){
                        state_beta[cnt2] = 0;
                        state_beta[cnt3] = 1;
                        unsigned int cnt_beta_bis = StateToCnt[ConvertType2(state_beta)];
                        RDMdown[cnt2+L*cnt3] += ((Nswaps%2==0)?1:-1) * bra_coeff * diagresult[pointer[cnt]+cnt_alpha+gBinomial(L,Nup)*cnt_beta_bis];
                        state_beta[cnt3] = 0;
                        state_beta[cnt2] = 1;
                     } else {
                        Nswaps++;
                     }
                  }
               }
            }
            
            //AddTwoParticles
            if (cnt>0){ //So the bra state can be created by adding two particles to a lower particle number state state.
               int Nswaps_alpha = 0;
               for (int cnt2=L-1; cnt2>=0; cnt2--){
                  if (state_alpha[cnt2]==1){
                     state_alpha[cnt2] = 0;
                     unsigned int cnt_alpha_bis = StateToCnt[ConvertType2(state_alpha)];
                     int Nswaps = Nswaps_alpha;
                     for (int cnt3=0; cnt3<L; cnt3++){
                        if (state_beta[cnt3]==1){
                           state_beta[cnt3] = 0;
                           unsigned int cnt_beta_bis = StateToCnt[ConvertType2(state_beta)];
                           AddTwoParticles[cnt2+L*cnt3] += ((Nswaps%2==0)?1:-1)*bra_coeff*diagresult[pointer[cnt-1]+cnt_alpha_bis+gBinomial(L,Nup-1)*cnt_beta_bis];
                           state_beta[cnt3] = 1;
                           Nswaps++;
                        }
                     }
                     state_alpha[cnt2] = 1;
                     Nswaps_alpha++;
                  }
               }
            }

         }
      }
      delete [] state_alpha;
      delete [] state_beta;
   }
   
   //Copy to lower diagonal
   for (int irow=0; irow<L; irow++){
      for (int icol=irow+1; icol<L; icol++){
         RDMup[icol+L*irow] = RDMup[irow+L*icol];
         RDMdown[icol+L*irow] = RDMdown[irow+L*icol];
      }
   }

}

