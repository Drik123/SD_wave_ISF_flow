# SD_wave_ISF_flow
Repository for the manuscript titled "Ionic waves drive fluid flow by inducing spatiotemporal porosity variations in the brain interstitial space" undergoing peer review (preprint available on bioRxiv: https://www.biorxiv.org/content/10.64898/2026.11.18.766825v1)

**(1) SD wave code:**

  **File type:** Fortran90
  
  **Description:** Simulates a physiologically realistic SD wave with initiation near the top right corner of an x-y planar domain of brain cortex. Instructions to change any of the parameters are commented in the code itself.
  
  **Output:** Possible outputs are  (i) intracellular K+ concentration field and (ii) intracellular Na+ concentration field, (iii) extracellular K+ concentration field and (iv) extracellular Na+ concentration field, (v) membrane action potential field, (vi) porosity field, (vii) temporal derivative of porosity, (viii) pressure field, (ix) x-velocity field, and (x) y-velocity field.


**(2) Theoretical wave:**

**File type:** Fortran90

  **Description:** Simulates a prescribed porosity wave in an x-y planar domain of brain cortex. Instructions to change any of the parameters are commented in the code itself.
  
  **Output:** Possible outputs are (i) porosity field, (ii) temporal derivative of porosity, (iii) pressure field, (iv) x-velocity field, and (v) y-velocity field,

**(3) Compile and run instructions:**
    Linear algebra packages llapack and lblas are needed. The following lines need to be executed to compile and run the .f90 code.

    > gfortran SD_wave_code.f90 -o -llapack -lblas -mcmodel=large
    > ./a.out

**(4) Citation:** 

**Description:**  Please cite the preprint available on bioRxiv (https://www.biorxiv.org/content/10.64898/2026.11.18.766825v1).
