!Control file for nrattle (data lines are free format)
! (note that text lines preceded by "!" are ignored)
!Revision of program involving a change in the control file on this date:
   09/09/11
!Stem name for output files (the file containing model information and complex spectra will be named 
! <stem>.nrattle_complex_spectra.out 
! and the one containing the absolute spectral amplitudes used in plotting will be named 
! <stem>.nrattle_amps4plot.out
!  nrattle_bj97_vs30_1100_ai_45
  test_nrattle_02mar11
!Number of frequencies (including f = 0.0) and highest freq. (nfreq, fhigh):
   46977   100
!Depth at which response is computed (0.0 for free surface): 
    0
!Model information:
!For each layer- thickness, Vs, density, 1/Q (if > 1.0, the program assumes
!that the value represents Q; note that damping in percent = 100*0.5*1/Q):
   1   10   75   1.2   12.5
!  Halfspace Vs and density (1/Q is automatically set to 0.0):
    435  2.7 
!Layer corresponding to halfspace (can be less than actual halfspace)   
! (a large number means the source is in the halfspace) and    
! angle of incidence (< 0.0 means obtain theta from site_amp file):   
! angle of incidence (< 0.0 means obtain theta from site_amp file):    
   81  0
