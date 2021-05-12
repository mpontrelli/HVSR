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
   8688    10
!Depth at which response is computed (0.0 for free surface):   
    0
!Model information:
!For each layer- thickness, Vs, density, 1/Q (if > 1.0, the program assumes
!that the value represents Q; note that damping in percent = 100*0.5*1/Q):
    1   1   150.586   1.7   0.033
    2   1   161.9635   1.7   0.033
    3   1   165.783   1.7   0.033
    4   1   178.7812   1.7   0.033
    5   1   181.7867   1.7   0.033
    6   1   196.1409   1.7   0.033
    7   1   196.1409   1.7   0.033
    8   1   205.4657   1.7   0.033
    9   1   207.9657   1.7   0.033
    10   1   207.9657   1.7   0.033
    11   1   205.2977   1.7   0.033
    12   1   202.1097   1.7   0.033
    13   1   201.8936   1.7   0.033
    14   1   199.9706   1.7   0.033
    15   1   199.5859   1.7   0.033
    16   1   207.7685   1.7   0.033
    17   1   208.5377   1.7   0.033
    18   1   206.9992   1.7   0.033
    19   1   205.4608   1.7   0.033
    20   1   206.6146   1.7   0.033
    21   1   228.8198   1.7   0.033
    22   1   229.2044   1.7   0.033
    23   1   230.1422   1.7   0.033
    24   1   230.9114   1.7   0.033
    25   1   230.1422   1.7   0.033
    26   1   255.6313   1.7   0.033
    27   1   254.6082   1.7   0.033
    28   1   252.1082   1.7   0.033
    29   1   255.0249   1.7   0.033
    30   1   254.1915   1.7   0.033
    31   1   273.3364   1.7   0.033
    32   1   270.8985   1.7   0.033
    33   1   270.8985   1.7   0.033
    34   1   283.8144   1.7   0.033
    35   1   287.3925   1.7   0.033
    36   1   325.927   1.7   0.033
    37   1   341.0594   1.7   0.033
    38   1   339.6308   1.7   0.033
    39   1   349.6308   1.7   0.033
    40   1   360.4026   1.7   0.033
    41   1   374.2417   1.7   0.033
    42   1   346.3844   1.7   0.033
    43   1   361.7305   1.7   0.033
    44   1   361.7305   1.7   0.033
    45   1   361.7305   1.7   0.033
    46   1   399.4549   1.7   0.033
    47   1   324.6939   1.7   0.033
    48   1   324.6939   1.7   0.033
    49   1   324.6939   1.7   0.033
    50   1   324.6939   1.7   0.033
!  Halfspace Vs and density (1/Q is automatically set to 0.0):
    3000   2.75
!Layer corresponding to halfspace (can be less than actual halfspace)   
! (a large number means the source is in the halfspace) and 
! angle of incidence (< 0.0 means obtain theta from site_amp file):    
   81  0
