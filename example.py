import PyBCS as BCS
import numpy as np

def doHubbard(L,TwoSz,U,mu_up,mu_down):

   #First set the Hamiltonian parameters
   Tmat = np.zeros([L,L])
   for cnt in range(0,L-1): #Hence OBC
      Tmat[cnt,cnt+1] = -1.0
      Tmat[cnt+1,cnt] = -1.0
   Delta = np.zeros([L,L])
   Uarray = U * np.ones([L])

   #Let the BCS library do its work
   Problem = BCS.PyBCS(L,TwoSz)
   Problem.ClearHam() #Not strictly necessary as ALL Hamiltonian parameters are set in the following
   Problem.SetTmat(Tmat)
   Problem.SetDelta(Delta)
   Problem.SetUarray(Uarray)
   Problem.SetMuUp(mu_up)
   Problem.SetMuDown(mu_down)
   Energy = Problem.Solve()
   RDMup = Problem.GetRDMup()
   RDMdown = Problem.GetRDMdown()
   AddTwoPart = Problem.GetAddTwoParticles()

   #Print what there is to print
   print "###  Hubbard  ###"
   print "Min. eigenvalue                                                 = ", Energy
   Nup = np.trace(RDMup)
   Ndown = np.trace(RDMdown)
   NormAdd2Part = np.linalg.norm(AddTwoPart)
   print "Norm of matrix with expectation values for adding two electrons = ",NormAdd2Part
   print "Number of electrons                                             = ",Nup + Ndown
   print "Energy (with chemical potential removed)                        = ",Energy+mu_up*Nup+mu_down*Ndown


doHubbard(8, 2, 9.0, 3.0, 3.0)

def doRandomSingleParticle(L):

   #First set the Hamiltonian parameters
   Tmat = -np.random.random([L,L])
   for cnt in range(0,L):
      for cnt2 in range(cnt+1,L):
         Tmat[cnt,cnt2] = Tmat[cnt2,cnt]
   Delta = np.random.random([L,L])
   mu_up = np.random.random([1])[0]
   mu_down = np.random.random([1])[0]
   
   #Let the BCS library do its work
   TwoSz = 0 #For the subsequent Bogoliubov part to remain valid.
   Problem = BCS.PyBCS(L,TwoSz)
   Problem.ClearHam()
   Problem.SetTmat(Tmat)
   Problem.SetDelta(Delta)
   Problem.SetMuUp(mu_up)
   Problem.SetMuDown(mu_down)
   BCS_Energy = Problem.Solve()
   RDMup = Problem.GetRDMup()
   RDMdown = Problem.GetRDMdown()
   AddTwoPart = Problem.GetAddTwoParticles()
   
   #Check the consistency of the solution.
   print "###  Single particle BCS  ###"
   print "Energy with BCS solver                                                  = ",BCS_Energy
   RDM_Energy = np.sum(Tmat * (RDMup + RDMdown)) - mu_up * np.trace(RDMup) - mu_down * np.trace(RDMdown) + 2.0 * np.sum(Delta * AddTwoPart)
   print "Energy recalculated based on the 1-RDMs and <a+ a+>                     = ",RDM_Energy
   
   #Check by diagonalizing the single particle Hamiltonian.
   Offset = np.trace(Tmat) - 0.5*(mu_up + mu_down)*L
   temp = np.zeros([2*L,2*L])
   temp[0:L,0:L] = Tmat - mu_up * np.diag(np.ones([L]))
   temp[L:2*L,L:2*L] = - Tmat + mu_down * np.diag(np.ones([L]))
   temp[0:L,L:2*L] = Delta
   temp[L:2*L,0:L] = np.array(np.mat(Delta).T)
   eigs, eigvecs = np.linalg.eigh(temp)
   SP_Energy = Offset + 0.5*np.sum(eigs[0:L]) - 0.5*np.sum(eigs[L:2*L])
   print "Energy by a Bogoliubov transormation of the single particle Hamiltonian = ",SP_Energy

   
doRandomSingleParticle(8)

