% The MTF Is the modulus of the OPTICAL TRANSFER FUNCTION computed from the point spread function
% OTF is the 2-D FFT of the point spread function

[pupilfunc PixelsInPupil] = Zphase_MahajanOSA(c, PARAMS, PARAMS.PupilSize); %Use the entire pupil for pupilfieldsize

OTF = conv2(pupilfunc,conj(fliplr(flipud((pupilfunc)))));%autocorrelation

% The MTF is the magnitude of the Optical Transfer Function
MTF = abs(OTF)'; 
MTF = MTF./max(max(MTF)); % scale the MTF so that the peak (at the origin) is equal to 1


% The PTF is the phase of the Optical Transfer Function
PTF = angle(OTF)';

PTF = PTF*180/pi; %convert PTF into degrees

% Select the subsection within the cutoff frequency for plotting the MTF

cutofffull = 2*(PARAMS.PupilSize)/(PARAMS.ImagingWavelength/1000)/57.3;
cutoffaper = 2*(PARAMS.PupilSize)/(PARAMS.ImagingWavelength/1000)/57.3;
halfsize = ceil(0.5*(cutoffaper/cutofffull)*PARAMS.PixelDimension);

MTFaper=MTF;
cutoff=cutofffull/2;

% Define the axes of the MTF plot in cyc/degree
n = size(MTFaper);
axisMTF = -cutoff:2*cutoff/n(1):(cutoff-(2*cutoff/n(1)));

figure('Position',[1200 500 400 400]);
mesh(axisMTF,axisMTF, MTF);colorbar;
title('Modulation Transfer Function');
xlabel('c/deg');
ylabel('c/deg');
zlabel('modulation');
xlim([-cutoff cutoff]);
ylim([-cutoff cutoff]);
colormap(hot)

if(1)%Save the MTF to a BMP and CVS file

    if isunix %return 'true' if you are using a MAC (use forward slash for subdirectories)
        name=['images/' PARAMS.FileName(1:length(PARAMS.FileName)-4) '_' num2str(DefocusInDiopters) '_astig_axis_' num2str(PupilSize) '_MTF.tif'];
    else %if you are using Windows (use backward slash for subdirectories)
        name=['images\' PARAMS.FileName(1:length(PARAMS.FileName)-4) '_' num2str(DefocusInDiopters) '_astig_axis_' num2str(PupilSize) '_MTF.tif'];
    end
    saveas(gcf,[name(1:length(name)-4) '.bmp']);
    csvwrite([name(1:length(name)-4) '.csv'],MTF);
end

if(1)
v=-195:30:195;
figure('Position',[0 0 400 400]);
imagesc(axisMTF,axisMTF,PTF);colorbar;
title('Phase Transfer Function');
xlabel('c/deg');
ylabel('c/deg');
zlabel('modulation');
axis square;
xlim([-cutoff cutoff]);
ylim([-cutoff cutoff]);
colormap(hsv);
end

if(1)    %Save the PTF to a BMP and CSV file

    if isunix %return 'true' if you are using a MAC (use forward slash for subdirectories)
        name=['images/' PARAMS.FileName(1:length(PARAMS.FileName)-4) '_' num2str(DefocusInDiopters) '_astig_axis_' num2str(PupilSize) '_PTF.tif'];
    else %if you are using Windows (use backward slash for subdirectories)
        name=['images\' PARAMS.FileName(1:length(PARAMS.FileName)-4) '_' num2str(DefocusInDiopters) '_astig_axis_' num2str(PupilSize) '_PTF.tif'];
    end
    saveas(gcf,[name(1:length(name)-4) '.bmp']);
    csvwrite([name(1:length(name)-4) '.csv'],PTF);
end

%This script computes and displays the average radial and horizontal and vertical MTFs and plots them on the same figure
%This script can only be used after running compMTF.m

[sizMTF junk] = size(MTFaper);
[junk peaky]=max(max(MTFaper));
[junk peakx]=max(max(transpose(MTFaper)));

maxradius=sqrt(2*(sizMTF/2)^2);

%MTFarray = zeros(sizMTF,sizMTF);
numbins=100;%sizMTF/4; %number of bins in the radial MTF
radialMTFarray = zeros(numbins,2);

for (i=1:sizMTF)
   for (j=1:sizMTF)
      xpos=i-peakx;
      ypos=j-peaky;
      radius=ceil(numbins*sqrt(xpos^2+ypos^2)/maxradius);
      if (radius==0)
      else
         radialMTFarray(radius,1)=radialMTFarray(radius,1)+MTFaper(i,j);
         radialMTFarray(radius,2)=radialMTFarray(radius,2)+1;
      end
   end
end


radialMTF=radialMTFarray(:,1)./radialMTFarray(:,2);

radialMTF=cat(1,1,radialMTF);
axisradMTF=0:sqrt(2)*max(abs(axisMTF))/numbins:sqrt(2)*max(abs(axisMTF));
figure('Position',[400 0 400 400]);
plot(axisradMTF,radialMTF,axisMTF(peakx+1:n),MTFaper(peakx+1:n,peaky),axisMTF(peaky+1:n),MTFaper(peakx,peaky+1:n));
legend('radial avg', 'horizontal', 'vertical');
title('MTF');
xlabel('c/deg');
ylabel('modulation');

if(1)%Save the Horizontal, vertical and radial average MTF plots to a BMP file. 

    if isunix %return 'true' if you are using a MAC (use forward slash for subdirectories)
        name=['images/' PARAMS.FileName(1:length(PARAMS.FileName)-4) '_' num2str(DefocusInDiopters) '_astig_axis_' num2str(PupilSize) '_HVRMTF.tif'];
    else %if you are using Windows (use backward slash for subdirectories)
        name=['images\' PARAMS.FileName(1:length(PARAMS.FileName)-4) '_' num2str(DefocusInDiopters) '_astig_axis_' num2str(PupilSize) '_HVRMTF.tif'];
    end
    saveas(gcf,[name(1:length(name)-4) '.bmp']);
end