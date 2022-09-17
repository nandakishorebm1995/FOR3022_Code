# Damage Identification in Fiber Metal Laminates using Bayesian Analysis with Model Order Reduction
 
### Description of the directories: 

* EnKF: 
     + Implements the ensemble Kalman filter method as a sequential data assimilation algorithm. It requires all the data from the TU BS cloud storage along with all the MATLAB files within this directory and to successfully conduct the test for damage identification.
 
* MCMC: 
    + Implements the Markov chain Monte Carlo method with Metropolis-Hastings algorithm in a Bayesian framework. It requires all the data from the TU BS cloud storage along with all the MATLAB files within this directory and to successfully conduct the test for damage identification.
    
* PMOR: 
    + Implements adaptive parametric model order reduction method that reduces the computation time by a factor of approximately 33 compared to that of the high-fidelity simulation using COMSOL-Multiphysics software. It requires the COMSOL model ("Stahl_laminat_032021_mesh.mph") from the TU BS cloud storage along with all the MATLAB files within this directory and to successfully generate the required number of global reduced order bases (ROBs).
    
### Description of the files:

* FML_FEM_HDM.m:
    + It takes damage parameters as inputs, generates an FEM model on COMSOL, solves using a time-dependent "Generalized alpha method" or Newmark-time integration method and outputs the high-fidelity solution of the system

* FML_FEM_ROM.m:
    + This function reduces the HiFi model of the underlying system. It consumes system matrices, simulation parameters and generated global ROBs and outputs the system response. 
    + It calls the "NewmarkIntegration3.m" function which actually performs the Newmark-time integration to solve the system. 
    
* FML_FEM_MATRIX_EXTRACT.m: 
    + This function takes damage parameters and excitation signal as the inputs, generates an FEM model on COMSOL and outputs the system matrices (mass matrix, stiffness matrix and load matrix) without solving the system. 
    
* NewmarkIntegration3.m: 
    + This function implements the implicit Newmark-beta time integration method. It is an explicit integration scheme. Here, alpha and beta are the constants that are set to 0.25 and 0.5 respectively, such that the scheme remains unconditionally stable.
