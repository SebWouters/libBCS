#ifndef BCS_H
#define BCS_H

/** BCS class.
    \author Sebastian Wouters <sebastian.wouters@ugent.be>
    \date November 14, 2012
  
    \section secBCSClass_1 Goal of the BCS class
 
    The goal of this class is to find the exact groundstate of the non-particle conserving hermitian BCS Hamiltonian in the site basis: \n
    \f$ \hat{H} = \sum\limits_{ij} t_{ij} \sum\limits_{\sigma} \hat{a}_{i \sigma}^{\dagger} \hat{a}_{j \sigma} + \sum\limits_{i} U_i \hat{n}_{i \uparrow} \hat{n}_{i \downarrow} - \mu_{\uparrow} \sum\limits_{i} \hat{n}_{i \uparrow} - \mu_{\downarrow} \sum\limits_{i} \hat{n}_{i \downarrow} + \sum\limits_{ij} \Delta_{ij} \left( \hat{a}_{i \uparrow}^{\dagger} \hat{a}_{j \downarrow}^{\dagger} + \hat{a}_{j \downarrow} \hat{a}_{i \uparrow} \right) \f$ \n
    where the latin letters denote site-indices \f$ \left(i = 0, 1 , ... , L-1\right) \f$ and the greek letters spin projections \f$ \left(\sigma = \uparrow, \downarrow\right) \f$. All Hamiltonian parameters should be real. Hermiticity requires the t-matrix to be symmetric. There are no requirements on the \f$ \Delta \f$-matrix. This Hamiltonian preserves spin projection, but not particle number. */
class BCS {

   public:
   
      //! Constructor
      /** \param L The number of sites in the model
          \param TwoSz Two times the targetted \f$ s_z \f$ value */      
      BCS(int L, int TwoSz);
      
      //! Destructor
      virtual ~BCS();
      
      //! Getter of the pointer to the t-matrix. See detailed description for conventions.
      /**  \return The pointer to the t-matrix. This matrix is stored in column major: \f$ t_{ij} = Telem[i + L*j] \f$. (This however doesn't matter as the t-matrix should be symmetric.) The indices i and j run from 0 to L-1. See Hamiltonian in \ref secBCSClass_1 "BCS class information". */
      double * gTelem();
      
      //! Getter of the pointer to the U-array. See detailed description for conventions.
      /** \return The pointer to the U-array. The index runs from 0 to L-1. See Hamiltonian in \ref secBCSClass_1 "BCS class information". */
      double * gUelem();
      
      //! Getter of the pointer to the \f$ \Delta \f$-matrix. See detailed description for conventions.
      /** \return The pointer to the \f$ \Delta \f$-matrix. This matrix is stored in column major: \f$ \Delta_{ij} = Delta[i + L*j] \f$. The indices i and j run from 0 to L-1. See Hamiltonian in \ref secBCSClass_1 "BCS class information". */
      double * gDelta();
      
      //! Getter of the pointer to the \f$ \mu_{\uparrow} \f$-variable.
      /** \return The pointer to the \f$ \mu_{\uparrow} \f$-variable. See Hamiltonian in \ref secBCSClass_1 "BCS class information". */
      double * gmu_up();
      
      //! Getter of the pointer to the \f$ \mu_{\downarrow} \f$-variable.
      /** \return The pointer to the \f$ \mu_{\downarrow} \f$-variable. See Hamiltonian in \ref secBCSClass_1 "BCS class information". */
      double * gmu_down();
      
      //! Solve the Hamiltonian. Returns the ground state energy. Please set ALL Hamiltonian parameters first by getting the pointers and filling in the values at the correct places.
      /** \return The ground state energy of the defined Hamiltonian (L, \f$ t \f$, \f$ U \f$, \f$ \Delta \f$, \f$ \mu_{\uparrow} \f$, \f$ \mu_{\downarrow} \f$) in the \f$ s_z \f$ symmetry sector, see Hamiltonian in \ref secBCSClass_1 "BCS class information". Please set ALL Hamiltonian parameters first by getting the pointers and filling in the values at the correct places. This also includes zero values. */
      double Solve();
      
      //! Getter of the pointer to the reduced density matrix of alpha electrons, see detailed description. Solve() must be called first.
      /** \return Pointer to the RDM of alpha electrons: \f$ \braket{\hat{a}_{i \uparrow}^{\dagger} \hat{a}_{j \uparrow}}\f$. This matrix is stored in column major: \f$ \braket{\hat{a}_{i \uparrow}^{\dagger} \hat{a}_{j \uparrow}} = array[i + L*j] \f$. (This however doesn't matter as the 1-RDM is symmetric.) The indices i and j run from 0 to L-1. */
      double * gRDMup();
      
      //! Getter of the pointer to the reduced density matrix of beta electrons, see detailed description. Solve() must be called first.
      /** \return Pointer to the RDM of beta electrons: \f$ \braket{\hat{a}_{i \downarrow}^{\dagger} \hat{a}_{j \downarrow}}\f$. This matrix is stored in column major: \f$ \braket{\hat{a}_{i \downarrow}^{\dagger} \hat{a}_{j \downarrow}} = array[i + L*j] \f$. (This however doesn't matter as the 1-RDM is symmetric.) The indices i and j run from 0 to L-1. */
      double * gRDMdown();
      
      //! Getter of the pointer to the matrix of expectation values for adding two particles, see detailed description. Solve() must be called first.
      /** \return Pointer to the matrix of expectation values for adding two particles: \f$ \braket{\hat{a}_{i \uparrow}^{\dagger} \hat{a}_{j \downarrow}^{\dagger}}\f$. This matrix is stored in column major: \f$  \braket{\hat{a}_{i \uparrow}^{\dagger} \hat{a}_{j \downarrow}^{\dagger}} = array[i + L*j] \f$. The indices i and j run from 0 to L-1. */
      double * gAddTwoParticles();
      
      //! Getter of the pointer to the array of local double occupancy expectation values, see detailed description. Solve() must be called first.
      /** \return Pointer to the array of local double occupancy expectation values: \f$ \braket{\hat{n}_{i\uparrow}\hat{n}_{i\downarrow}} = array[i] \f$. The index i runs from 0 to L-1. */
      double * gDoubleOcc();
   
   private:
   
      //Problem size parameters
      int L;
      int TwoSz;
      int NN;
      int * Nvals;
      
      //Function that determines the max. L-value for which no integer overflows occur on the current machine.
      int gMaxL();
      
      //Hamiltonian parameters --> have to be set to zero by user when zero (column major (fortran) for Telem and Delta).
      double * Telem;
      double * Uelem;
      double * Delta;
      double mu_up;
      double mu_down;
      
      //Binomial coefficient stuff.
      unsigned int * Binomial;
      void FillBinomial();
      unsigned int gBinomial(int, int);
      
      //Memory helper
      unsigned long * pointer;
      
      //State to counter conversions: what matters for the general part of the program
      //See comments at ConvertType2/3 for information on the states used.
      unsigned int ** CntToState;
      unsigned int * StateToCnt;
      unsigned int ConvertType2(int *);
      void ConvertType3(int *, unsigned int);
      
      //State to counter conversions: at the beginning of the program
      void FillStateCounterConversions();
      bool getNextState(int *, int);
      
      //Hamiltonian times a FCI vector
      void HtimesVec(double *, double *);
      
      //The result of the diagonalization
      double * diagresult;
      double diagEnergy;
      
      //Density matrix stuff
      void CalculateRDMs();
      double * RDMup;
      double * RDMdown;
      double * AddTwoParticles;
      double * DoubleOcc;
      

};

#endif
