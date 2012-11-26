#include "BCS.h"

extern "C"{

   BCS * BCS_create( int L, int TwoSz ){ return new BCS(L, TwoSz); }
  
   void BCS_destroy( BCS * prob ){ delete prob; }
  
   void BCS_ClearHam( BCS * prob, int L ){
      for (int cnt=0; cnt<L; cnt++){
         for (int cnt2=0; cnt2<L; cnt2++){
            prob->gTelem()[cnt+L*cnt2] = 0.0;
            prob->gDelta()[cnt+L*cnt2] = 0.0;
         }
         prob->gUelem()[cnt] = 0.0;
      }
      *(prob->gmu_up()) = 0.0;
      *(prob->gmu_down()) = 0.0;
   }
  
   void BCS_setTelem( BCS * prob, int L, int i, int j, double val){
      prob->gTelem()[i+L*j] = val;
   }
  
   void BCS_setDelta( BCS * prob, int L, int i, int j, double val){
      prob->gDelta()[i+L*j] = val;
   }
  
   void BCS_setUelem( BCS * prob, int i, double val){
      prob->gUelem()[i] = val;
   }
  
   void BCS_setmu_up( BCS * prob, double val){
      *(prob->gmu_up()) = val;
   }
  
   void BCS_setmu_down( BCS * prob, double val){
      *(prob->gmu_down()) = val;
   }
   
   double BCS_Solve( BCS * prob){
      return prob->Solve();
   }
  
   double BCS_getRDMup( BCS * prob, int L, int i, int j){
      return prob->gRDMup()[i+L*j];
   }
  
   double BCS_getRDMdown( BCS * prob, int L, int i, int j){
      return prob->gRDMdown()[i+L*j];
   }
  
   double BCS_getAddTwoParticles( BCS * prob, int L, int i, int j){
      return prob->gAddTwoParticles()[i+L*j];
   }
   
   double BCS_getDoubleOcc( BCS * prob, int L, int i){
      return prob->gDoubleOcc()[i];
   }
 
}  

