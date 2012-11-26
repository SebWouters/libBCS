import numpy as np
from ctypes import *
libBCS = cdll.LoadLibrary('./libBCS.so')
#Set the argtypes
libBCS.BCS_create.argtypes = [c_int,c_int]
libBCS.BCS_destroy.argtypes = [c_void_p]
libBCS.BCS_ClearHam.argtypes = [c_void_p, c_int]
libBCS.BCS_setTelem.argtypes = [c_void_p,c_int,c_int,c_int,c_double]
libBCS.BCS_setDelta.argtypes = [c_void_p,c_int,c_int,c_int,c_double]
libBCS.BCS_setUelem.argtypes = [c_void_p,c_int,c_double]
libBCS.BCS_setmu_up.argtypes = [c_void_p,c_double]
libBCS.BCS_setmu_down.argtypes = [c_void_p,c_double]
libBCS.BCS_Solve.argtypes = [c_void_p]
libBCS.BCS_getRDMup.argtypes = [c_void_p,c_int,c_int,c_int]
libBCS.BCS_getRDMdown.argtypes = [c_void_p,c_int,c_int,c_int]
libBCS.BCS_getAddTwoParticles.argtypes = [c_void_p,c_int,c_int,c_int]
#Set the restype
libBCS.BCS_Solve.restype = c_double
libBCS.BCS_getRDMup.restype = c_double
libBCS.BCS_getRDMdown.restype = c_double
libBCS.BCS_getAddTwoParticles.restype = c_double
libBCS.BCS_getDoubleOcc.restype = c_double

## PyBCS class
#
#  \author Sebastian Wouters <sebastian.wouters@ugent.be>
#  \date November 15, 2012
#    
#  \section secPyBCSgoal Goal of the PyBCS class
#  
#  The C++ class BCS, together with the C-bindings in source/BCS_cfunctions.cpp, can be compiled into the library libBCS.so. This class provides an easy interface to the library.
#
#  \section secPyBCSham The Hamiltonian
#
#  For the sake of completeness, we repeat here the Hamiltonian for which libBCS finds the exact groundstate: the non-particle conserving hermitian BCS Hamiltonian in the site basis: \n
#  \f$ \hat{H} = \sum\limits_{ij} t_{ij} \sum\limits_{\sigma} \hat{a}_{i \sigma}^{\dagger} \hat{a}_{j \sigma} + \sum\limits_{i} U_i \hat{n}_{i \uparrow} \hat{n}_{i \downarrow} - \mu_{\uparrow} \sum\limits_{i} \hat{n}_{i \uparrow} - \mu_{\downarrow} \sum\limits_{i} \hat{n}_{i \downarrow} + \sum\limits_{ij} \Delta_{ij} \left( \hat{a}_{i \uparrow}^{\dagger} \hat{a}_{j \downarrow}^{\dagger} + \hat{a}_{j \downarrow} \hat{a}_{i \uparrow} \right) \f$ \n
#  where the latin letters denote site-indices \f$ \left(i = 0, 1 , ... , L-1\right) \f$ and the greek letters spin projections \f$ \left(\sigma = \uparrow, \downarrow\right) \f$. All Hamiltonian parameters should be real. Hermiticity requires the t-matrix to be symmetric. There are no requirements on the \f$ \Delta \f$-matrix. This Hamiltonian preserves spin projection, but not particle number.
class PyBCS(object):

   ## Constructor of a PyBCS object.
   #  The constructor creates a C++ BCS object.
   #  \param L The number of sites in the model
   #  \param TwoSz Two times the targetted \f$ s_z \f$ value
   def __init__(self,L,TwoSz):
      ## The number of sites in the model
      self.L = L
      ## Two times the targetted \f$ s_z \f$ value
      self.TwoSz = TwoSz
      ## Pointer to the C++ BCS object
      self.obj = libBCS.BCS_create(L,TwoSz)
   
   ## Destructor of the PyBCS object.
   #  Very important to destruct the BCS object located at self.obj!
   def __del__(self):
      libBCS.BCS_destroy(self.obj)

   ## Clear all Hamiltonian parameters in the BCS object.
   def ClearHam(self):
      libBCS.BCS_ClearHam(self.obj, self.L)
   
   ## Set the hopping terms.
   #  \param Tmat Tmat should be a 2D numpy array of size L by L, containing the hopping terms between the different sites: \f$ t_{ij} = Tmat[i,j] \f$.
   def SetTmat(self, Tmat):
      for irow in range(0,self.L):
         for icol in range(0,self.L):
            libBCS.BCS_setTelem(self.obj, self.L, irow, icol, Tmat[irow,icol])
            
   ## Set the \f$\Delta\f$-matrix.
   #  \param Delta Delta should be a 2D numpy array of size L by L: \f$\Delta_{ij} = Delta\left[i,j\right]\f$.
   def SetDelta(self, Delta):
      for irow in range(0,self.L):
         for icol in range(0,self.L):
            libBCS.BCS_setDelta(self.obj, self.L, irow, icol, Delta[irow,icol])
   
   ## Set the U-array.
   #  \param Uarray Uarray should be a 1D numpy array of size L.
   def SetUarray(self, Uarray):
      for irow in range(0,self.L):
         libBCS.BCS_setUelem(self.obj, irow, Uarray[irow])
   
   ## Set \f$\mu_{\uparrow}\f$.
   #  \param mu_up The value for \f$\mu_{\uparrow}\f$.
   def SetMuUp(self, mu_up):
      libBCS.BCS_setmu_up(self.obj, mu_up)
 
   ## Set \f$\mu_{\downarrow}\f$.
   #  \param mu_down The value for \f$\mu_{\downarrow}\f$.
   def SetMuDown(self, mu_down):   
      libBCS.BCS_setmu_down(self.obj, mu_down)
   
   ## Find the ground state.
   #  Note that the Hamiltonian should be set first.
   #  \return The ground state energy.
   def Solve(self):
      return libBCS.BCS_Solve(self.obj)
   
   ## Get the 1-RDM for the alpha electrons.
   #  \return The 1-RDM of the alpha electrons \f$\braket{\hat{a}_{i\uparrow}^{\dagger}\hat{a}_{j\uparrow}}\f$.
   #  Note that the ground state should be computed first by Solve(). It is resturned as a 2D numpy array of size L by L: \f$\braket{\hat{a}_{i\uparrow}^{\dagger}\hat{a}_{j\uparrow}} = array[i,j]\f$.
   def GetRDMup(self):
      A = np.zeros([self.L,self.L])
      for irow in range(0,self.L):
         for icol in range(0, self.L):
            A[irow,icol] = libBCS.BCS_getRDMup(self.obj, self.L, irow, icol)
      return A
      
   ## Get the 1-RDM for the beta electrons.
   #  \return The 1-RDM of the beta electrons \f$\braket{\hat{a}_{i\downarrow}^{\dagger}\hat{a}_{j\downarrow}}\f$.
   #  Note that the ground state should be computed first by Solve(). It is resturned as a 2D numpy array of size L by L: \f$\braket{\hat{a}_{i\downarrow}^{\dagger}\hat{a}_{j\downarrow}} = array[i,j]\f$.
   def GetRDMdown(self):
      A = np.zeros([self.L,self.L])
      for irow in range(0,self.L):
         for icol in range(0, self.L):
            A[irow,icol] = libBCS.BCS_getRDMdown(self.obj, self.L, irow, icol)
      return A
   
   ## Get the two-particle addition probability.
   #  \return The two-particle addition probability \f$\braket{\hat{a}_{i\uparrow}^{\dagger}\hat{a}_{j\downarrow}^{\dagger}}\f$.
   #  Note that the ground state should be computed first by Solve(). It is resturned as a 2D numpy array of size L by L: \f$\braket{\hat{a}_{i\uparrow}^{\dagger}\hat{a}_{j\downarrow}^{\dagger}} = array[i,j]\f$.
   def GetAddTwoParticles(self):
      A = np.zeros([self.L,self.L])
      for irow in range(0,self.L):
         for icol in range(0, self.L):
            A[irow,icol] = libBCS.BCS_getAddTwoParticles(self.obj, self.L, irow, icol)
      return A
      
   ## Get the local double occupancy expectation value.
   #  \return The local double occupancy expectation value \f$\braket{\hat{n}_{i\uparrow}\hat{n}_{i\downarrow}}\f$.
   #  Note that the ground state should be computed first by Solve(). It is resturned as a 1D numpy array of size L: \f$\braket{\hat{n}_{i\uparrow}\hat{n}_{i\downarrow}} = array[i]\f$.
   def GetDoubleOcc(self):
      A = np.zeros([self.L])
      for icnt in range(0,self.L):
         A[icnt] = libBCS.BCS_getDoubleOcc(self.obj, self.L, icnt)
      return A



