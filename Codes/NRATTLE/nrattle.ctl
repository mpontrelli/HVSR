!Control file for nrattle (data lines are free format)
!Revision of program involving a change in the control file on this date:
   09/09/11
!Stem name for output files (the file containing model information and complex spectra will be named 
! <stem>.nrattle_complex_spectra.out 
! and the one containing the absolute spectral amplitudes used in plotting will be named 
! <stem>.nrattle_amps4plot.out
!  nrattle_bj97_vs30_1100_ai_45
  test_nrattle_02mar11
!Number of frequencies (including f = 0.0) and highest freq. (nfreq, fhigh):   
   10000    20
!Depth at which response is computed (0.0 for free surface):   
    0
!Model information:
!For each layer- thickness, Vs, density, 1/Q (if > 1.0, the program assumes
!that the value represents Q; note that damping in percent = 100*0.5*1/Q):
    1   1   151.5445   2.7   20
    2   1   155.9295   2.7   20
    3   1   160.3145   2.7   20
    4   1   164.6995   2.7   20
    5   1   169.0845   2.7   20
    6   1   173.4695   2.7   20
    7   1   177.8545   2.7   20
    8   1   182.2395   2.7   20
    9   1   186.6245   2.7   20
    10   1   191.0095   2.7   20
    11   1   195.3945   2.7   20
    12   1   199.7795   2.7   20
    13   1   204.1645   2.7   20
    14   1   208.5495   2.7   20
    15   1   212.9345   2.7   20
    16   1   217.3195   2.7   20
    17   1   221.7045   2.7   20
    18   1   226.0895   2.7   20
    19   1   230.4745   2.7   20
    20   1   234.8595   2.7   20
    21   1   239.2445   2.7   20
    22   1   243.6295   2.7   20
    23   1   248.0145   2.7   20
    24   1   252.3995   2.7   20
    25   1   256.7845   2.7   20
    26   1   261.1695   2.7   20
    27   1   265.5545   2.7   20
    28   1   269.9395   2.7   20
    29   1   274.3245   2.7   20
    30   1   278.7095   2.7   20
    31   1   283.0945   2.7   20
    32   1   287.4795   2.7   20
    33   1   291.8645   2.7   20
    34   1   296.2495   2.7   20
    35   1   300.6345   2.7   20
    36   1   305.0195   2.7   20
    37   1   309.4045   2.7   20
    38   1   313.7895   2.7   20
    39   1   318.1745   2.7   20
    40   1   322.5595   2.7   20
    41   1   326.9445   2.7   20
!  Halfspace Vs and density (1/Q is automatically set to 0.0):
    2000   2.7
!Layer corresponding to halfspace (can be less than actual halfspace)   
! (a large number means the source is in the halfspace) and 
! angle of incidence (< 0.0 means obtain theta from site_amp file):    
   81  0
