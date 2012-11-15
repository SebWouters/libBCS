#include <iostream>
#include <cstdlib>
#include <time.h>

using std::cout;
using std::endl;

#include "BCS.h"

// Sebastian Wouters, November 14, 2012.

BCS::BCS(int L, int TwoSz){

   srand(time(NULL));

   int maxL = gMaxL();
   if (L>maxL){
      cout << "Maximum value of L that is supported on this machine = " << maxL << endl;
      cout << "Integer overflow possible. Edit code!" << endl;
      exit(666);
   }

   //System size stuff: Depending on Sz and L, there are NN possible particle number sectors, with possible values in Nvals
   this->L = L;
   this->TwoSz = TwoSz;
   int Nmin = abs(TwoSz);
   NN = L - Nmin + 1;
   Nvals = new int[NN];
   for (int cnt=0; cnt<NN; cnt++){
      Nvals[cnt] = Nmin + 2*cnt;
   }
   
   //Allocate memory for the matrix elem stuff
   Telem = new double[L*L];
   Delta = new double[L*L];
   Uelem = new double[L];
   
   //Find out how much space each particle space needs: for Nvals[cnt], the space needed is pointer[cnt+1]-pointer[cnt]
   Binomial = new unsigned int[(L+1)*(L+1)];
   FillBinomial();
   pointer = new unsigned long[NN+1];
   pointer[0] = 0;
   for (int cnt=0; cnt<NN; cnt++){
      pointer[cnt+1] = pointer[cnt] + gBinomial(L,(Nvals[cnt]+TwoSz)/2)*gBinomial(L,(Nvals[cnt]-TwoSz)/2);
   }
   
   //Allocate and fill the state to counter conversions
   CntToState = new unsigned int*[L+1];
   for (int cnt=0; cnt<=L; cnt++){
      CntToState[cnt] = new unsigned int[gBinomial(L,cnt)];
   }
   StateToCnt = new unsigned int[(1 << L)]; // (1<<L) = 2^L
   FillStateCounterConversions();
   
   //Allocate result memory
   diagresult = new double[pointer[NN]];
   RDMup = new double[L*L];
   RDMdown = new double[L*L];
   AddTwoParticles = new double[L*L];

}

BCS::~BCS(){

   delete [] Nvals;

   delete [] Telem;
   delete [] Delta;
   delete [] Uelem;
   
   delete [] Binomial;
   delete [] pointer;
   
   delete [] StateToCnt;
   for (int cnt=0; cnt<=L; cnt++){
      delete [] CntToState[cnt];
   }
   delete [] CntToState;
   
   delete [] diagresult;
   delete [] RDMup;
   delete [] RDMdown;
   delete [] AddTwoParticles;

}

void BCS::FillBinomial(){ //From wikipedia :-)

   for (int cnt=0; cnt<=L; cnt++){
      Binomial[(L+1)*cnt] = 0;
      Binomial[cnt] = 1;
   }
   
   for (int n=1; n<=L; n++){
      for (int k=1; k<=n; k++){
         Binomial[n + (L+1)*k] = gBinomial(n-1,k-1) + gBinomial(n-1,k);
      }
   }

}

unsigned int BCS::gBinomial(int n, int k){

   if ((n<0) || (k<0) || (k>n) || (n>L)){
      return 0;
   }
   
   return Binomial[n + (L+1)*k];

}

void BCS::FillStateCounterConversions(){

   //0 particles --> only Cnt=0, State=0
   CntToState[0][0] = 0;
   StateToCnt[0] = 0;
   
   //L particles --> only Cnt=0, State=2^L-1
   unsigned int TwoToPowerLMinusOne = (1<<L)-1;
   CntToState[L][0] = TwoToPowerLMinusOne;
   StateToCnt[TwoToPowerLMinusOne] = 0;
   
   //Npart particles
   for (int Npart=1; Npart<=L-1; Npart++){
   
      //Start with Counter = 0 and all particles to the left; state contains the ordered positions of the Npart particles.
      unsigned int Counter = 0;
      int * state = new int[Npart];
      for (int cnt=0; cnt<Npart; cnt++){
         state[cnt] = cnt;
      }
      
      //Loop over the different Npart states
      bool stop = false;
      do{
         unsigned int ConvertedState = 0;
         for (int cnt=0; cnt<Npart; cnt++){
            ConvertedState += (1<<state[cnt]); //(1<<n) = 2^n 
         }
         CntToState[Npart][Counter] = ConvertedState;
         StateToCnt[ConvertedState] = Counter;
         stop = getNextState(state,Npart);
         Counter++;
      } while (!stop);
      
      if (Counter!=gBinomial(L,Npart)){
         cout << "Error at FillStateCounterConversions: Nstates = " << Counter << " vs. C(" << L << "," << Npart << ")" << endl;
      }
   
      delete [] state;
   }

}

//State is of length Npart and contains the Npart particle positions in an ordered way, with pos = 0 ... L-1.
//For 5 particles for example the state [1,2,3,6,7] would generate [0,1,4,6,7] as the next state.
//The bool is false if a new state was generated and true if no new states can be generated.
bool BCS::getNextState(int * state, int Npart){

   int index = 0;
   
   while (index<Npart-1){
      if (state[index]==state[index+1]-1){
         index += 1;
      } else {
         state[index] += 1;
         for (int cnt=0; cnt<index; cnt++){
            state[cnt] = cnt;
         }
         return false; //don't stop yet
      }
   }
   
   //now index==Npart-1
   if (state[index]<L-1){
      state[index] += 1;
      for (int cnt=0; cnt<index; cnt++){
         state[cnt] = cnt;
      }
      return false; //don't stop yet
   }
   
   //no updates possible anymore
   return true; //stop

}

//State has length L and contains zeros and ones corresponding to resp. empty and full sites, thereturn is the corresponding integer: sum_i state[i] * 2^i.
unsigned int BCS::ConvertType2(int * state){

   unsigned int thereturn = 0;
   for (int cnt=0; cnt<L; cnt++){
      if (state[cnt]==1){
         thereturn += (1<<cnt);
      }
   }
   return thereturn;

}

//Convert theinteger = sum_i state[i] * 2^i into state (opposite action of ConvertType2)
void BCS::ConvertType3(int * state, unsigned int theinteger){

   for (int cnt=0; cnt<L; cnt++){
      state[cnt] = (theinteger & (1<<cnt))?1:0;
   }

}

double * BCS::gTelem(){ return Telem; }
double * BCS::gUelem(){ return Uelem; }
double * BCS::gDelta(){ return Delta; }
double * BCS::gmu_up(){ return &mu_up; }
double * BCS::gmu_down(){ return &mu_down; }

int BCS::gMaxL(){

   int maxL = 1;
   while ((1<<maxL) < (1<<(maxL+1))){
      maxL ++;
   }
   return maxL;

}

double * BCS::gRDMup(){ return RDMup; }
double * BCS::gRDMdown(){ return RDMdown; }
double * BCS::gAddTwoParticles(){ return AddTwoParticles; }


