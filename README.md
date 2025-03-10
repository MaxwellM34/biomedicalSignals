# biomedicalSignals

This repository contains MATLAB code that processes EEG data from two sources, EEGrelaxed.csv and Grandmal.txt, to compare normal EEG activity with seizure activity and to implement a basic automated seizure detection algorithm. The code reads EEGrelaxed.csv, which includes elapsed time and channels FP1, FP2, O1, and O2, and plots its time-domain and frequency-domain features using FFT analysis to characterize the baseline EEG. It then loads Grandmal.txt, a single-channel recording of seizure activity, and segments the data into one-second windows with a 50 percent overlap. For each window, the code computes the short-time energy (the sum of squared amplitudes) and establishes a detection threshold set to the mean energy plus two standard deviations. Windows that exceed this threshold are flagged as containing seizure activity. The script produces plots that display the energy levels for each window as well as a binary detection indicator over time, providing visual justification for the seizure detection approach. This code requires MATLAB with the Signal Processing Toolbox and has been tested with MATLAB R2018b and later.
