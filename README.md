This tool was designed to analyze the ciliary beating frequency from videos recorded with a high-speed camera (~1000 frame/s).  
  
It allows users to set multiple region-of-interest (ROIs) following the lines. The size and spacing of ROIs can be adjusted by the user, and the number of ROIs is determined accordingly. Intensity changes across time were extracted from every pixel in individual ROIs and used to compute the Fast Fourier Transform to analyze the frequency components. The spectrum of mean frequency components was used to compute the representative frequency of each ROI's intensity change. These values were stored and plotted for validation.  
  
This code is written in MATLAB.  
  
Please refer to Kang et al., 2024, "Loss of Katnal2 leads to ependymal ciliary hyperfunction and autism-related phenotypes in mice", Plos Biology, for further information.  
  
# Example
The video for example was recorded ex vivo for the movements of cilia located in the lateral ventricular ependymal cell layer under oxygenated and pH-adjusted conditions. We provide 4 lines across the ependymal cell linings for ROI generation. 33 ROIs of size 11 x 11 pixels, spaced 30 pixels apart, were generated and used in the example.
