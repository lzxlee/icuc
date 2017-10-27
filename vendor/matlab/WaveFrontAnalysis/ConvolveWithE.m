function [TumbE convolution]=ConvolveWithE(LetterSize,PARAMS,PSF)

plotdimension = 60*PARAMS.PixelDimension*(180*60/3.1416)*PARAMS.ImagingWavelength*.001/PARAMS.PupilFieldSize;
LetterSizePixels = PARAMS.PixelDimension*LetterSize*60/plotdimension;

boxsize = PARAMS.PixelDimension;

TumbE=ones(boxsize,boxsize);
imcenter=size(TumbE)./2;

gapsize=round(LetterSizePixels/5);

TumbE=ones(boxsize,boxsize);
TumbE(imcenter(1)-2*gapsize:imcenter(1)+3*gapsize-1,imcenter(2)-2*gapsize:imcenter(2)+3*gapsize-1)=0;
TumbE(imcenter(1)-1*gapsize:imcenter(1)-0*gapsize-1,imcenter(2)-1*gapsize:imcenter(2)+3*gapsize-1)=1;
TumbE(imcenter(1)+1*gapsize:imcenter(1)+2*gapsize-1,imcenter(2)-1*gapsize:imcenter(2)+3*gapsize-1)=1;

product =fftshift(fft2(double(TumbE))) .* fftshift(fft2(double(PSF)));
convolution = abs(fftshift(ifft2(double(product))));
clear product;

convolution = convolution./max(max(convolution));

axisCONV=-plotdimension/2:plotdimension/PARAMS.PixelDimension:(plotdimension/2)-plotdimension/PARAMS.PixelDimension;

scrsz = get(0,'ScreenSize');
figure('Position',[800 500 400 400]);
colormap(gray);
imagesc(axisCONV,axisCONV,convolution);
title('Convolution');
xlabel('arcsec');
ylabel('arcsec');
axis square;

if(0)
    if isunix %MAC
        fname=['images/convE.tif'];
        fname=['images/E.tif'];
    else %windows
      fname=['images\E.tif'];
      fname=['images\convE.tif'];
    end
    imwrite(TumbE,fname,'tif');
    imwrite(convolution,fname,'tif');
end