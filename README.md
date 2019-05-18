# HRV detrend
Matlab implementation of Heart Rate Variability (HRV) detrending method by MP Tarvainen et al., 2002

## Usage
* detrendFast divides the signal in *order* windows and compute the detrend accordingly
* detrendSample divides the signal according to the window length parameter and detrends the last window separately. Better option for signals with different or unknown number of samples
