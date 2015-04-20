# kexcalc
A ZZ-exchange global fitting program
Owen Williams

1. Introduction

Kexcalc is a Matlab function written to analyze NMR data from a ZZ-exchange experiment. It uses time series data for the
auto and cross peaks and calculates the exchange rate that produces the best fit to the experimental data. The fitting
equations were taken from Palmer et al. (1) and Kleckner and Foster (2).

2. Requirements

Kexcalc requires a working installation of Matlab. It should work on any Windows, Mac, or Linux system as it does not use
any OS-specific functions; however, it has only been tested on Windows 7 (as that's what I have).

3. Installation

Download kexcalc.m and place it in a directory on your Matlab path.

4. Usage

a. Import your peak intensities and mixing times into Matlab. Make sure these are all the same length!
b. Call kexcalc and optionally assign its output to a variable, e.g.
    results = kexcalc(Iaa, Ibb, Iab, I,ba, times, [guesses]);
  where Iaa and Ibb are vectors with intensity data for the auto peaks, Iab and Iba are vectors of intensity data for
  the cross peaks, times is a vector of the corresponding mixing times, and guesses is a structure with non-default
  initial guesses for the exchange parameters.
c. The fitting results are automatically plotted and printed to the console, but are also saved in a result structure
  for later use.
  
5. Structures

a. guesses: starting point for the minimization function. Override the defaults here.
  guesses.kex: exchange rate
  guesses.pA: initial population of peak A
  guesses.R1A: longitudinal relaxation rate of peak A
  guesses.R1B: longitudinal relaxation rate of peak B
  
b. results:
  results.kex: exchange rate
  results.pA: population of peak A
  results.pB: population of peak B
  results.R1A: longitudinal relaxation rate of peak A
  results.R1B: longitudinal relaxation rate of peak B
  results.datafit: matrix of normalized intensity data and fitted points, in the format
    [dataAA dataBB dataAB dataBA fitAA fitBB fitAB fitBA]
  results.fitIAA: global fit to the A auto peak intensites
  results.fitIBB: global fit to the B auto peak intensites
  results.fitIAB: global fit to the AB cross peak intensites
  results.fitIBA: global fit to the BA cross peak intensites
  results.dataIAA: normalized intensites for auto peak A
  results.dataIBB: normalized intensites for auto peak B
  results.dataIAB: normalized intensites for cross peak AB
  results.dataIBA: normalized intensites for cross peak BA
  results.residual: Sum over all time points of the squared difference between the experimental data and the fit.

6. References

(1) Palmer, A. G.; Kroenke, C. D.; Loria, J. P. "Nuclear magnetic resonance methods for quantifying 
microsecond-to-millisecond motions in biological molecules." Meth. Enzymology. 2001, 339, 204-238.

(2) Kleckner, I. R.; Foster, M. P. "An introduction to NMR-based approaches for measuring protein dynamics." 
Biochimica et Biophysica Acta. 2011, 1814, 942-968.
  
7. Contact Information

Questions, comments, or suggestions? Feel free to email me at olwilliams256@gmail.com.
