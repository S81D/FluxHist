# FluxHist
ANNIE Flux plots for CC and NC XS.

Will output Flux histogram plot as well as a root file with histograms for each flavor that can be used in combination with the GENIE XS spline files to calculate predicted cross sections from GENIEv3.

Note the spline files are in energy bins of 50 MeV. Therefore the default histogram binning here is recommended for 50 MeV. Modify the configurations (start of the code) accordingly. The default units are [neutrinos / 50 MeV / cm^2 / POT]; these units can be modified within the code.

### Configuration
(Modify in `FluxRates.py` header)
```
plot_name = 'example.png'        # name of output Flux Histogram plot
use_full_detector = True         # whether to use full detector (True) volume or a reduced, FV (False) --> FV is defined below the header and can be modified accordingly
generate_rootfile = False        # generate an output root file with flux histograms for each flavor
rootfile_name = 'example.root'   # name of the output root file
generate_summary_txt = True      # generate a summary txt with flux statistics
txt_name = 'example.txt'         # name of the output txt file
TEST_RUN = True                  # will only run over 20 gsimple files (used for fast debugging and testing)
gsimple_directory = 'filepath/'  # path for the gsimple flux files
```

### Usage
`python3 FluxRates.py`

### Output
1. `root` file with the flux histograms per flavor
2. Printed flux statistics
3. Output txt file with flux statistics
