# FluxHist
ANNIE Flux plots for CC and NC XS.

Will output Flux histogram plot as well as a root file with histograms for each flavor that can be used in combination with the GENIE XS spline files to calculate predicted cross sections from GENIEv3.

Note the spline files are in energy bins of 50 MeV. Therefore the default histogram binning here is recommended for 50 MeV. Modify the configurations (start of the code) accordingly.

Example plots are included with different units.

### Usage
`python3 FluxRates.py`
