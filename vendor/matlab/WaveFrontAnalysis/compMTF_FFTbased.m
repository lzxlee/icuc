% The MTF Is the modulus of the OPTICAL TRANSFER FUNCTION computed from the point spread function
% OTF is the 2-D FFT of the point spread function

OTF = fft2(PSF);

% The MTF is the real part of the Optical Transfer Function
MTF = fftshift(abs(OTF)); % FFT shift sets the origin to the center of the MTF plot
MTF = MTF./max(max(MTF)); % scale the MTF so that the peak (at the origin) is equal to 1

% Select the subsection within the cutoff frequency for plotting the MTF

cutofffull = (PARAMS.PupilFieldSize)/(PARAMS.ImagingWavelength/1000)/57.3;
cutoffaper = (PARAMS.PupilSize)/(PARAMS.ImagingWavelength/1000)/57.3;
halfsize = ceil(0.5*(cutoffaper/cutofffull)*PARAMS.PixelDimension);

%if(halfsize<PARAMS(1))
%   MTFaper = MTF(PARAMS(1)/2-halfsize:PARAMS(1)/2+halfsize,PARAMS(1)/2-halfsize:PARAMS(1)/2+halfsize);
%	cutoff=cutoffaper;   
%else
   MTFaper=MTF;
   cutoff=cutofffull/2;
%end

% Define the axes of the MTF plot in cyc/degree
n = size(MTFaper);
axisMTF = -cutoff:2*cutoff/n(1):(cutoff-(2*cutoff/n(1)));

figure('Position',[1200 500 400 400]);
surf(axisMTF,axisMTF,MTFaper);
title('Modulation Transfer Function');
xlabel('c/deg');
ylabel('c/deg');
zlabel('modulation');
xlim([-cutoff cutoff]);
ylim([-cutoff cutoff]);
colormap(hsv)
clear OTF MTF;

