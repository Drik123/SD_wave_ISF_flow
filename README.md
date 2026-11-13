# SD_wave_ISF_flow
Repository for the paper titled "Ionic waves drive fluid flow by inducing spatiotemporal variations in porosity of the brain interstitial space"

**(1) SD wave code:**

  **File type:** Fortran90
  
  **Description:** Simulates a real SD wave with initiation at the top right corner of an x-y planar domain of brain cortex. Instructions to change any of the parameters are commented in the code itself.
  
  **Output:** Possible outputs are intracellular (i) K+ conc. and (ii) Na+ conc., extracellular (iii) K+ conc. and (iv) Na+ conc., (v) membrane action potential, (vi) porosity variation, (vii) temporal derivative of porosity change, (viii) pressure variation, (ix) x-velocity, and (x) y-velocity.


**(2) Theoretical wave:**

**File type:** Fortran90

  **Description:** Simulates a theoretical porosity wave which varies sinusoidally and travels diagonally from the bottom left to the top right in an x-y planar domain of brain cortex. Instructions to change any of the parameters are commented in the code itself.
  
  **Output:** Possible outputs are (i) porosity variation, (ii) temporal derivative of porosity change, (iii) pressure variation, (iv) x-velocity, and (v) y-velocity,
