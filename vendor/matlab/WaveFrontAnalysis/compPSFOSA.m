function PSF = CompPSFOSA(c,PARAMS);


% Zphase_Mahajan generates the complex pupil function for the Fourier transform

[pupilfunc PixelsInPupil] = Zphase_MahajanOSA(c,PARAMS,PARAMS.PupilFieldSize);

pupilfunc=transpose(pupilfunc);

% The amplitude of the point spread function is the Fourier transform
% of the wavefront aberration (when is it expressed in terms of phase)
Hamp=fft2(pupilfunc);

% The intensity of the PSF is the square of the amplitude
% the complex conjugate is a way of multiplying out the complex part of the function

Hint=(Hamp .* conj(Hamp));

% Define the size of the PSF plot in arcmin.
% NOTE: The dimension of a single pixel in the PSF in radians is the wavelength
% divided by the size of the pupil field.

plotdimension = 60*PARAMS.PixelDimension*(180*60/3.1416)*PARAMS.ImagingWavelength*.001/PARAMS.PupilFieldSize;

%fprintf('The size of the point spread image is %g seconds of arc\n',plotdimension);

PSF = real(fftshift(Hint)); % this comment reorients the PSF so the origin is at the center of the image
PSF = PSF./(PixelsInPupil^2); % scale the PSF so that peak represents the Strehl ratio
clear Hint;
%fprintf('The Strehl ratio is: %d\n',max(max(PSF)));


% Compute the dimension of one side of the PSF
% NOTE: The dimension of a single pixel in the PSF image (in radians) is the wavelength
% divided by the size of the pupil field.

axisPSF=-plotdimension/2:plotdimension/PARAMS.PixelDimension:(plotdimension/2)-plotdimension/PARAMS.PixelDimension;

scrsz = get(0,'ScreenSize');
% Set parameters and display the wavefront aberration as a grayscale image
iptsetpref('ImshowAxesVisible','on');;
figure('Position',[400 500 400 400]);
%imshow(PSF, [0 max(max(PSF))],'XData',axisPSF,'YData',axisPSF);
imagesc(axisPSF,axisPSF,PSF);
%surfc(axisPSF,axisPSF,PSF);
title('Point Spread Function');
xlabel('arcsec');
ylabel('arcsec');
axis square;
zlabel('height (relative to diffraction limited peak)');
xlim([-plotdimension/2 plotdimension/2]);
ylim([-plotdimension/2 plotdimension/2]);
colormap(gray);

