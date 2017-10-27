% Optical performance.m
% Austin Roorda, June 20, 1999 (with many updates)
% This program is used to generate plots of the 
% wavefront, PSF, PTF and MTF from a set of preset coefficients or a file of Zernike coefficients

% The ordering of the Zernike polynomials is based on OSA standard ordering;

% *******************************************************************************************
% ************ PART 1: SET UP THE PROGRAM WITH APPROPRIATE PARAMETERS ***********************
% ************** and decide what to run  ********************************
% *******************************************************************************************
close all; clear all;

PARAMS.PixelDimension = 256;% size of pupil aperture field in pixels (this defines the resolution of the calculation)
PARAMS.PupilSize = 5;
PARAMS.PupilFitSize = 5;
PARAMS.PupilFieldSize = 15;
PARAMS.ImagingWavelength = 0.550;% imaging wavelength in microns
PARAMS.WavefrontResolution = 20;% increase to enhance the display of the wavefront (doesn't affect calculation)

saveWF = 1;  %set to zero to skip this part of the code
savePSF = 1; %set to zero to skip this part of the code
runCONV = 1; %set to zero to skip this part of the code
runMTFPTF = 1; %set to zero to skip this part of the code
saveMTF = 0;

% *******************************************************************************************
% ************ PART 2: GET THE ZERNIKE COEFFICIENTS *****************************************
% *******************************************************************************************

choice = inputdlg('1: Preset coeffs 2:*.zer file 3: random eye','Choose 1,2 or 3');
choice = str2num(choice{1});%-- 8/13/09  9:35 AM --%

if (choice==1) %use preset coefficients   
    % Zernike coeffs
    % Preset all coefficients to 0 and edit each value as desired

    % tilt
    c(1)=0; 	c(2)=0.0;
    %defocus and astigmatism
    c(3)=0; 	c(4)=0;	c(5)=0;
    %coma like
    c(6)=0.0;	c(7)=0.0;	c(8)=0.0;	c(9)=0.0;	
    %spherical aberration like
    c(10)=0;	c(11)=0.0;	c(12)=0.0;	c(13)=0.0;	c(14)=0;	
    %higher order (each row is a new radial order)
    c(15)=0;	c(16)=0;	c(17)=0.1;	c(18)=0;	c(19)=0;	c(20)=0;	
    c(21)=0;	c(22)=0;	c(23)=0;	c(24)=0;	c(25)=0.0;	c(26)=0;	c(27)=0;	
    c(28)=0;	c(29)=0;	c(30)=0;	c(31)=0;	c(32)=0;	c(33)=0;	c(34)=0;	c(35)=0;	
    c(36)=0;	c(37)=0;	c(38)=0;	c(39)=0;	c(40)=0;	c(41)=0;	c(42)=0;	c(43)=0;	c(44)=0;	
    c(45)=0;	c(46)=0;	c(47)=0;	c(48)=0;	c(49)=0;	c(50)=0;	c(51)=0;	c(52)=0;	c(53)=0;	c(54)=0;	
    c(55)=0;	c(56)=0;	c(57)=0;	c(58)=0;	c(59)=0;	c(60)=0;	c(61)=0;	c(62)=0;	c(63)=0;	c(64)=0;	c(65)=0;
    PARAMS.PupilSize = 5;% size of pupil in mm for which PSF and MTF is to be calculated
    PARAMS.PupilFitSize = 5;% size of pupil in mm that Zernike coefficients define. NOTE: You can define the aberrations
                    % for any pupil size and calculate their effects for any smaller aperture.

    % Because of the reciprocal relationship between the pupil function and its Fourier transform,
    % it helps to define a pupil that makes up only a small central region of the pupil aperture.
    % Otherwise the point spread function is too small.
    PARAMS.PupilFieldSize = 20;	param4orig = PARAMS.PupilFieldSize; 	% size of pupil field in mm (use a large field to magnify the PSF)

    %boxsize = 145.6; %set parameter 4 to get a specific PSF image dimension;
    %PARAMS.PupilFieldSize = 60*PARAMS.PixelDimension*(180/3.1416)*PARAMS.ImagingWavelength*.001/boxsize;
    PARAMS.FileName = 'preset';    
    fprintf('\nYou are working with manually set coefficients\n');    
    
elseif (choice ==2) % Read a list of Zernike coefficients
      
    [fname,pname] = uigetfile('*.zer','Open Coefficient File');
    fid = fopen([pname fname],'r');
    version = fgetl(fid);% fscanf(fid,'%s',1);
    instrument = fgetl(fid);% fscanf(fid,'%s',1);
    manuf = fgetl(fid);% fscanf(fid,'%s',1);
    oper = fgetl(fid);% fscanf(fid,'%s',1);
    pupoff = fgetl(fid);% fscanf(fid,'%s',1);
    geooff = fgetl(fid);% fscanf(fid,'%s',1);
    datatype = fgetl(fid);% fscanf(fid,'%s',1);
    Rfit = fscanf(fid,'%s %f\n',[1 2]);% fscanf(fid,'%s',1);
    Rfit=Rfit(length(Rfit));
    Rmax = fgetl(fid);% fscanf(fid,'%s',1);
    waverms = fgetl(fid);% fscanf(fid,'%s',1);
    order = fgetl(fid);% fscanf(fid,'%s',1);
    strehl = fgetl(fid);% fscanf(fid,'%s',1);
    refent = fgetl(fid);% fscanf(fid,'%s',1);
    refcor = fgetl(fid);% fscanf(fid,'%s',1);
    resspec = fgetl(fid);% fscanf(fid,'%s',1);
    data = fgetl(fid);% fscanf(fid,'%s',1);
    c = fscanf(fid,'%i %i %g',[3 inf]); %read the first line
    fclose(fid);
    c=c(3,2:66); %ignore the piston term (check this carefully!!!!!)

    %set these parameter based on what is contained in the *.zer file
    PARAMS.PupilSize=2*Rfit; %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
    PARAMS.PupilFitSize=2*Rfit; 
    PARAMS.PupilFieldSize=20; %automatically compute the field size
    fprintf('\nThe current *.zer filename is %s\n',fname);
    PARAMS.FileName = fname;
elseif (choice == 3)
    randcoeffs=VirtualEyes(1,6); %make it for a 6 mm pupil
    c=zeros(1,65);
    c(1,1:length(randcoeffs)-1) = randcoeffs(2:36)';
    PARAMS.FileName = 'Random';
    PARAMS.PupilSize = 6;
    PARAMS.PupilFitSize = 6;
end

PARAMS % display the parameters to the screen

% *******************************************************************************************
% ************ PART 3: COMPUTE THE SPECTACLE CORRECTION *************************************
% This section will calculate the best spectacle correction for the eye, considering all terms
% *******************************************************************************************
   
Prescription = spectacleFuncOSA(PARAMS,c); % call the function spectacleFuncOSA.m

fprintf('The best spectacle correction is %3.5g DS X %3.5g DC : AXIS %4g deg\n',Prescription.Sphere,Prescription.Cylinder,Prescription.Axis);
 
% *******************************************************************************************
% ************* PART 4: Calculate the RMS of the wave aberration*****************************
% In this section, you can calculate the RMS for all terms or any subest of terms. For example,
% to calculate the rms for defocus, simply caculate the sqrt of the sum of the squares for terms 
% 3,4 and 5. To calculate the RMS for 3rd order terms and higher, calculate the sqrt of the 
% sum of the squares for terms 6 - 66. 
% NOTE: If you change the values of any of the coefficients at this stage, those changes will
% be applied to all the next calculations. So, if you want to calculate the PSF and MTF for 
% high order aberrations only, then set all low order terms to zero at this point.
% *******************************************************************************************

% Remove the useless terms
c(1)=0;    % this term is tilt and does not have any effect on image quality
c(2)=0;	   % this term is tilt and does not have any effect on image quality

% remove other terms as desired

if(1) % set to one to set astigmatism to zero
    c(3)=0;    % remove astigmatism
    c(5)=0;    % remove astigmatism
    fprintf('Astigmatism has been set to zero for all calculations\n\n');
end

c(4) = 0;  % remove defocus

fprintf('Defocus\tPupil\tRMS\tStrehl\n');

corig=c;%write the original coefficients into a new variable

PupilSize = PARAMS.PupilSize; %change this if you want to analyze the wavefront for a new pupil size. You can also add this to the for loop below

%**************************************************************************
%** Notes on using the For loop********************************************
% The 'for' loop can be used to cycle through a range of variables. The two
% most common variable to cycle through are DefocusInDiopters and Pupilsize. 
% examples
%              DefocusInDiopters = -1:.25:1; << this runs a series of values of defocus from -1 D to 1 D in 0.25 D steps. 
%              PupilSize = 1:6; << this runs thr program through a series of pupil size values from 1 to 6 in 1 mm steps
%
% If you do not want to run through a series of values, then just set a
% fixed variable after the for command in one or both of the lines. 
% eg DefocusInDiopters = 0;

% Caution: never set PARAMS.PupilSize to have a value larger than the pupil
% size that the Zernike Coefficients define (PARAMS.PupilFitSize)
%*************************************************************************

for DefocusInDiopters=-1:1:1; % or use something like DefocusInDiopters = -1:0.25:1 to specify a range of defocus values. 

for PupilSize = PARAMS.PupilSize; % or use something like PupilSize = 1:1:6 to specify a range of pupil sizes. 
    
    PARAMS.PupilSize=PupilSize;
    
    if (PARAMS.PupilSize < PARAMS.PupilFitSize)
        c = TransformC(corig,PARAMS.PupilFitSize,PARAMS.PupilSize,0,0,0); %compute new coeffficients if pupil size has changed. 
        c(1:5) = 0; %set tilt, defocus and astigmatism to zero    
    else
        c=corig;
    end
    
    c(4)=(1e6/(4*sqrt(3)))*DefocusInDiopters*((PARAMS.PupilSize/2000)^2); % convert DefocusInDiopters into a Zernike coefficient

    % calculate the RMS
    rms=sqrt(sum(c(1:65).^2)); % compute the sum of the sqaure of all the coefficients in the list that you define. 

    % print the result to the screen
    fprintf('%g\t%g\t%g\t',DefocusInDiopters, PupilSize, rms);

    % *******************************************************************************************
    % ************ PART 5: GENERATE AND DISPLAY THE WAVEFRONT ABERRATION ************************
    % *******************************************************************************************
    
    waveabermap = compWaveOSA(PARAMS,c); % call the function compWaveOSA.m
    
    if(saveWF)%Save the wavefront map into a BMP and CSV file
        if isunix %return 'true' if you are using a MAC (use forward slash for subdirectories)
            name=['images/' PARAMS.FileName(1:length(PARAMS.FileName)-4) '_' num2str(DefocusInDiopters) '_astig_axis_' num2str(PupilSize) '_WF.tif'];
        else %if you are using Windows (use backward slash for subdirectories)
            name=['images\' PARAMS.FileName(1:length(PARAMS.FileName)-4) '_' num2str(DefocusInDiopters) '_astig_axis_' num2str(PupilSize) '_WF.tif'];
        end
        saveas(gcf,[name(1:length(name)-4) '.bmp']);
        csvwrite([name(1:length(name)-4) '.csv'],waveabermap);
    end
 
    % *******************************************************************************************
    % ************ PART 6: GENERATE AND DISPLAY THE POINT SPREAD FUNCTION ***********************
    % *******************************************************************************************
    
    PSF = compPSFOSA(c,PARAMS); % call the function compPSFOSA.m

    % write some parameters to the screen
    fprintf('%g\n',max(max(PSF))); %print the strehl ratio to the screen

    if(savePSF) %Save the point spread function into a TIF, BMP and CVS file. 
        if isunix %return 'true' if you are using a MAC (use forward slash for subdirectories)
            name=['images/' PARAMS.FileName(1:length(PARAMS.FileName)-4) '_' num2str(DefocusInDiopters) '_astig_axis_' num2str(PupilSize) '_PSF.tif'];
        else %if you are using Windows (use backward slash for subdirectories)
            name=['images\' PARAMS.FileName(1:length(PARAMS.FileName)-4) '_' num2str(DefocusInDiopters) '_astig_axis_' num2str(PupilSize) '_PSF.tif'];
        end          
        PSF=PSF/(max(max(PSF)));
        imwrite(PSF,name,'tif');
        saveas(gcf,[name(1:length(name)-4) '.bmp']);
        csvwrite([name(1:length(name)-4) '.csv'],PSF);
    end
    
    % *******************************************************************************************
    % ************ PART 7: GENERATE AND DISPLAY A CONVOLUTION OF AN E WITH THE PSF ****************
    % *******************************************************************************************
   
    if(runCONV) %convolve the current PSF with a letter E (size LetterSize) and save as BMP and CSV file. 
        LetterSize=5; %This is the lettersize in minutes of arc: examples --> 20/10=2.5; 20/20 = 5; 20/40 = 10 ...
        [TumbE ConvE] = ConvolveWithE(LetterSize,PARAMS,PSF);%/(max(max(PSF))));
        if isunix %return 'true' if you are using a MAC (use forward slash for subdirectories)
            name=['images/' PARAMS.FileName(1:length(PARAMS.FileName)-4) '_' num2str(DefocusInDiopters) '_astig_axis_' num2str(PupilSize) '_CONV.tif'];
        else %if you are using Windows (use backward slash for subdirectories)
            name=['images\' PARAMS.FileName(1:length(PARAMS.FileName)-4) '_' num2str(DefocusInDiopters) '_astig_axis_' num2str(PupilSize) '_CONV.tif'];
        end
        saveas(gcf,[name(1:length(name)-4) '.bmp']);
        csvwrite([name(1:length(name)-4) '.csv'],ConvE);
    end
    
    % *******************************************************************************************
    % ************ PART 8: GENERATE AND DISPLAY THE MODULATION TRANSFER FUNCTION ****************
    % *******************************************************************************************
    
    if(runMTFPTF)
        
        %compMTF_FFTbased;
        compMTFPTF_radhorvert;

        if(saveMTF)
            if isunix %return 'true' if you are using a MAC (use forward slash for subdirectories)
                namerad=['images/' PARAMS.FileName(1:length(PARAMS.FileName)-4) '_' num2str(DefocusInDiopters) '_astig_axis_' num2str(PupilSize) '_radMTF.csv'];
                nameHV=['images/' PARAMS.FileName(1:length(PARAMS.FileName)-4) '_' num2str(DefocusInDiopters) '_astig_axis_' num2str(PupilSize) '_HVMTF.csv'];   
            else %if you are using Windows (use backward slash for subdirectories)
                namerad=['images/' PARAMS.FileName(1:length(PARAMS.FileName)-4) '_' num2str(DefocusInDiopters) '_astig_axis_' num2str(PupilSize) '_radMTF.csv'];
                nameHV=['images/' PARAMS.FileName(1:length(PARAMS.FileName)-4) '_' num2str(DefocusInDiopters) '_astig_axis_' num2str(PupilSize) '_HVMTF.csv'];
            end
        
        radMTF = cat(2,axisradMTF',radialMTF);
        csvwrite(namerad,radMTF);
        
        HVMTF = cat(2,axisMTF(peaky+1:n)',MTFaper(peakx+1:n,peaky),MTFaper(peakx,peaky+1:n)');
        csvwrite(nameHV,HVMTF);
        end
        
    end
    drawnow;
end
end
