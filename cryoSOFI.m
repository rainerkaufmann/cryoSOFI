% cryoSOFI(d,'align','roi','quick')
%
%   INPUT:
%   d           raw data stack
%   'align'     Do drift correction.
%   'roi'       Usees ROI defined by user for drift determination.
%               --> click on top left corner, then bottom right corner of
%                   desired ROI, then press ENTER
%   'quick'     Uses only one set of taus [tau1, tau2, tau3, tau4].
%               Without 'quick' a permutation of all possibilities is used
%               and averaged.
%
%   EXAMPLE FOR USAGE: X = cryoSOFI(d,'align','quick');
%       --> Data will be corrected for drift based on whole field of view
%           before calculating SOFI images in the 'quick' mode.
%
%   OUTPUT:
%   should be mostly self explanatory...
%   AC: auto-correlation
%   XC: cross-correlation
%   XG4: 4th order cross-correlation that also contains lower order
%        information (while XC4 is the 4th order cumulant)
%   deconv: deconvolved
%
%
% Copyright 2016-2019 by Felipe Moser & Rainer Kaufmann
% All rights reserved.
%
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

function [ output ] = cryoSOFI( varargin )
d = varargin{1};
if any( strcmpi( varargin, 'align'))
    if any( strcmpi( varargin, 'roi'))
        roi = 1;
    else
        roi = 0;
    end
        d = DriftCorrection( d, roi );
end
d = single(d);
delta = bsxfun(@minus, d, mean(d,3));
average = mean(d,3);
tau1 = 2;
tau2 = 3;
tau3 = 4;
tau4 = 5;
quick = 0;


if nargin >1
    if any(strcmpi(varargin,'tau'))
        tauVarPos=find(strcmpi(varargin,'tau'))+1;
        tauVarIn=varargin{tauVarPos};
        if numel(tauVarIn)~=4
            msgbox('[ERROR] Tau has to be a vector with 4 elements');
            return;
        else
            tau1 = tauVarIn(1);
            tau2 = tauVarIn(2);
            tau3 = tauVarIn(3);
            tau4 = tauVarIn(4);
        end
    end
    if any( strcmpi( varargin, 'quick'))
        quick = 1;
    end
    type = varargin{2};
    if strcmpi(type, 'psf')
        [ measuredPSF, measuredPSFstd, horCor, vertCor, diagCor1, diagCor2 ] = pixelDistanceCorrection( delta, average, tau1, tau2, tau3, tau4, quick );
        output = struct('measuredPSF', measuredPSF, 'measuredPSFstd', measuredPSFstd,...
            'horCor', horCor, 'vertCor', vertCor, 'diagCor1', diagCor1, 'diagCor2', diagCor2);
        
    elseif strcmpi(type, 'auto')
        [  AC2, AC3, AG4, AC4 ] = AC_cryoSOFI( delta, tau1, tau2, tau3, tau4, quick);
        output = struct('average',average, 'AC2', AC2, 'AC3', AC3, 'AG4', AG4, 'AC4', AC4);
    elseif strcmpi(type, 'cross')
        [ ~, measuredPSFstd, horCor, vertCor, diagCor1, diagCor2 ] = pixelDistanceCorrection( delta, average, tau1, tau2, tau3, tau4, quick );
        [ XC2, XC3, XG4, XC4 ] = XC_cryoSOFI(  delta, tau1, tau2, tau3, tau4, quick);
        [ XC2, XC3, XG4, XC4 ] = XC_pixelDistanceCorrection (XC2, XC3, XG4, XC4,  horCor, vertCor, diagCor1, diagCor2);
        [ XC2deconv, XC3deconv, XG4deconv, XC4deconv, averageDeconv ] = XC_deconvolution (XC2, XC3, XG4, XC4, average, measuredPSFstd);
        output = struct('average', average, 'averageDeconv', averageDeconv, 'XC2', XC2, 'XC3', XC3, 'XG4', XG4, 'XC4', XC4,...
            'XC2deconv', XC2deconv, 'XC3deconv', XC3deconv,'XG4deconv', XG4deconv,'XC4deconv',XC4deconv, 'measuredPSFstd', measuredPSFstd);
    else
        [ AC2, AC3, AG4, AC4 ] = AC_cryoSOFI( delta, tau1, tau2, tau3, tau4, quick);
        [ ~, measuredPSFstd, horCor, vertCor, diagCor1, diagCor2 ] = pixelDistanceCorrection( delta, average, tau1, tau2, tau3, tau4, quick );
        [ XC2, XC3, XG4, XC4 ] = XC_cryoSOFI(  delta, tau1, tau2, tau3, tau4, quick);
        [ XC2, XC3, XG4, XC4 ] = XC_pixelDistanceCorrection (XC2, XC3, XG4, XC4,  horCor, vertCor, diagCor1, diagCor2);
        [ XC2deconv, XC3deconv, XG4deconv, XC4deconv, averageDeconv ] = XC_deconvolution (XC2, XC3, XG4, XC4, average, measuredPSFstd);
        output = struct('average', average, 'averageDeconv', averageDeconv, 'AC2', AC2, 'AC3', AC3, 'AG4', AG4, 'AC4', AC4,...
            'XC2', XC2, 'XC3', XC3, 'XG4', XG4, 'XC4', XC4,...
            'XC2deconv', XC2deconv, 'XC3deconv', XC3deconv,'XG4deconv', XG4deconv,'XC4deconv',XC4deconv, 'measuredPSFstd', measuredPSFstd);
    end
    
else
    [  AC2, AC3, AG4, AC4 ] = AC_cryoSOFI( delta, tau1, tau2, tau3, tau4, quick);
    [ ~, measuredPSFstd, horCor, vertCor, diagCor1, diagCor2 ] = pixelDistanceCorrection( delta, average, tau1, tau2, tau3, tau4, quick );
    [ XC2, XC3, XG4, XC4 ] = XC_cryoSOFI(  delta, tau1, tau2, tau3, tau4, quick);
    [ XC2, XC3, XG4, XC4 ] = XC_pixelDistanceCorrection (XC2, XC3, XG4, XC4,  horCor, vertCor, diagCor1, diagCor2);
    [ XC2deconv, XC3deconv, XG4deconv, XC4deconv, averageDeconv ] = XC_deconvolution (XC2, XC3, XG4, XC4, average, measuredPSFstd);
    output = struct('average', average, 'averageDeconv', averageDeconv, 'AC2', AC2, 'AC3', AC3, 'AG4', AG4, 'AC4', AC4,...
        'XC2', XC2, 'XC3', XC3, 'XG4', XG4, 'XC4', XC4,...
        'XC2deconv', XC2deconv, 'XC3deconv', XC3deconv,'XG4deconv', XG4deconv,'XC4deconv',XC4deconv, 'measuredPSFstd', measuredPSFstd);
end

end

function [ AC2, AC3, AG4, AC4 ] = AC_cryoSOFI ( delta, tau1, tau2, tau3, tau4, quick )
xpix = size(delta,1);
ypix = size(delta,2);
nimages = size(delta,3);
if quick == 1
    tau = [tau1, tau2, tau3, tau4];
else
    tau = perms( [tau1, tau2, tau3, tau4] );
end
tauMax = max( [tau1, tau2, tau3, tau4] );
i = 1:xpix;
j = 1:ypix;
t = 1:( nimages-tauMax );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%
%%%%%%%%%%% Calculate auto-correlation SOFI images
%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Calculating auto-correlation SOFI images...');



%%%%%%%%%%% Auto-correlations
AC2 = zeros(xpix, ypix);
AC3 = zeros(xpix, ypix);
AG4 = zeros(xpix, ypix);
AG4_12 = zeros(xpix, ypix);
AG4_13 = zeros(xpix, ypix);
AG4_14 = zeros(xpix, ypix);
AG4_23 = zeros(xpix, ypix);
AG4_24 = zeros(xpix, ypix);
AG4_34 = zeros(xpix, ypix);

for k = 1:size(tau,1)
    AC2 = AC2 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,2))),3)/size(tau,1);
    AC3 = AC3 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,2)).*delta(i,j,t+tau(k,3))),3)/size(tau,1);
    AG4 = AG4 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,2)).*delta(i,j,t+tau(k,3)).*delta(i,j,t+tau(k,4))),3)/size(tau,1);
    AG4_12 = AG4_12 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,2))),3)/size(tau,1);
    AG4_13 = AG4_13 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,3))),3)/size(tau,1);
    AG4_14 = AG4_14 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,4))),3)/size(tau,1);
    AG4_23 = AG4_23 + mean(abs(delta(i,j,t+tau(k,2)).*delta(i,j,t+tau(k,3))),3)/size(tau,1);
    AG4_24 = AG4_24 + mean(abs(delta(i,j,t+tau(k,2)).*delta(i,j,t+tau(k,4))),3)/size(tau,1);
    AG4_34 = AG4_34 + mean(abs(delta(i,j,t+tau(k,3)).*delta(i,j,t+tau(k,4))),3)/size(tau,1);
end

AC4 = abs(AG4-AG4_12.*AG4_34-AG4_13.*AG4_24-AG4_14.*AG4_23);


fprintf('Finished. \n');
end

function [ measuredPSF, measuredPSFstd, horCor, vertCor, diagCor1, diagCor2, sigma ] = pixelDistanceCorrection ( delta, average, tau1, tau2, tau3, tau4, quick)
xpix = size( delta,1 );
ypix = size( delta,2 );
nimages = size( delta,3 );
if quick == 1
    tau = [tau1, tau2, tau3, tau4];
else
    tau = perms( [tau1, tau2, tau3, tau4] );
end
tauMax = max( [tau1, tau2, tau3, tau4] );
i = 2:xpix-1;
j = 2:ypix-1;
t = 1:( nimages-tauMax );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%
%%%%%%%%%%% Calculate cross-correlation pixel distance factors
%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf( 'Calculating cross-correlation pixel distance factors...' );



%%%%%%%%%%% 2nd-order auto- and cross-correlations
AC2 = zeros(xpix-2, ypix-2);
XC2hor = zeros(xpix-2, ypix-2);
XC2vert = zeros(xpix-2, ypix-2);
XC2diag1 = zeros(xpix-2, ypix-2);
XC2diag2 = zeros(xpix-2, ypix-2);

for k=1:size(tau,1)
    AC2 = AC2 + mean( abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,2))),3)/size(tau,1);
    XC2hor = XC2hor + mean(abs(delta(i,j-1,t+tau(k,1)).*delta(i,j+1,t+tau(k,2))),3)/size(tau,1);
    XC2vert = XC2vert +mean(abs(delta(i-1,j,t+tau(k,1)).*delta(i+1,j,t+tau(k,2))),3)/size(tau,1);
    XC2diag1 = XC2diag1+ mean(abs(delta(i-1,j-1,t+tau(k,1)).*delta(i+1,j+1,t+tau(k,2))),3)/size(tau,1);
    XC2diag2 = XC2diag2 +mean(abs(delta(i-1,j+1,t+tau(k,1)).*delta(i+1,j-1,t+tau(k,2))),3)/size(tau,1);
end

%%%%%%%%%%% Filter to use only optimal parts of the image for calculating
%%%%%%%%%%% the pixel distance corrections.
averageBlur = imgaussfilt( average, 2);
averageBlur = averageBlur-min(averageBlur(:));
averageBlur = mat2gray(averageBlur);
filterGradient = bwconvhull(abs(gradient(averageBlur))>0.01, 'objects');
filterIntensity = (averageBlur>0.25 & averageBlur<0.9);
filter = logical(filterGradient.*filterIntensity);
filter = filter( 2:xpix-1, 2:ypix-1);

%%%%%%%%%%% Pixel distance correction factors
horCor=XC2hor./AC2;
horCor(horCor>1 | horCor<0 | isnan(horCor))=0;
horCor=mean(horCor(filter));
vertCor=XC2vert./AC2;
vertCor(vertCor>1 | vertCor<0 | isnan(vertCor))=0;
vertCor=mean(vertCor(filter));
diagCor1=XC2diag1./AC2;
diagCor1(diagCor1>1 | diagCor1<0 | isnan(diagCor1))=0;
diagCor2=XC2diag2./AC2;
diagCor2(diagCor2>1 | diagCor2<0 | isnan(diagCor2))=0;
diagCor1=mean(diagCor1(filter));
diagCor2=mean(diagCor2(filter));

%%%%%%%%%%% Standard deviation of the PSF from pixel distance correction factors
sigmaX = sqrt(-1/(log(horCor)));
sigmaY = sqrt(-1/(log(vertCor)));
sigma = (sigmaX + sigmaY)/2;

%%%%%%%%%%% Measured PSF of the original images calculated from the distance correction factors
measuredPSF = fspecial('gaussian', [xpix, ypix] , double(sigma));
measuredPSFstd = sigma;

fprintf('Finished. \n');
end

function [ XC2, XC3, XG4, XC4 ] = XC_cryoSOFI ( delta, tau1, tau2, tau3, tau4, quick)
delta = double(delta);
xpix = size(delta,1);
ypix = size(delta,2);
nimages = size(delta,3);
if quick == 1
    tau = [tau1, tau2, tau3, tau4];
else
    tau = perms( [tau1, tau2, tau3, tau4] );
end
tauMax = max( [tau1, tau2, tau3, tau4] );
i = 2:xpix-1;
j = 2:ypix-1;
t = 1:( nimages-tauMax );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%
%%%%%%%%%%% Calculate cross-correlation SOFI images
%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Calculating cross-correlation SOFI images...\n');


%%%%%%%%%%% 2nd-order cross-correlation
fprintf('Calculating 2nd-order cross-correlation SOFI image...');

XC2Central = zeros(xpix-2, ypix-2);
XC2Horizontal = zeros(xpix-2, ypix-2);
XC2Vertical = zeros(xpix-2, ypix-2);
XC2Diagonal1 = zeros(xpix-2, ypix-2);
XC2Diagonal2 = zeros(xpix-2, ypix-2);
for k = 1:size(tau,1)
    XC2Central = XC2Central + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,2))),3)/size(tau,1);
    XC2Horizontal = XC2Horizontal + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j+1,t+tau(k,2))),3)/size(tau,1);
    XC2Vertical = XC2Vertical + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j,t+tau(k,2))),3)/size(tau,1);
    XC2Diagonal1 = XC2Diagonal1 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j+1,t+tau(k,2))),3)/size(tau,1);
    XC2Diagonal2 = XC2Diagonal2 + mean(abs(delta(i,j+1,t+tau(k,1)).*delta(i+1,j,t+tau(k,2))),3)/size(tau,1);
end
XC2Diagonal = (XC2Diagonal1+XC2Diagonal2)/2;
XC2xpixOUT = 2*(xpix-2)-1;
XC2ypixOUT = 2*(ypix-2)-1;
XC2 = zeros(XC2xpixOUT, XC2ypixOUT);
XC2(1:2:end,1:2:end) = XC2Central(1:end,1:end);
XC2(1:2:end,2:2:end) = XC2Horizontal(1:end,1:end-1);
XC2(2:2:end,1:2:end) = XC2Vertical(1:end-1,1:end);
XC2(2:2:end,2:2:end) = XC2Diagonal(1:end-1,1:end-1);
clear XC2Central XC2Horizontal XC2Vertical...
    XC2Diagonal1 XC2Diagonal2 XC2Diagonal;
fprintf('Finished.\n');


%%%%%%%%%%% 3rd-order cross-correlation
fprintf('Calculating 3rd-order cross-correlation SOFI image...');

XC3Central = zeros(xpix-2, ypix-2);
XC3Horizontal1 = zeros(xpix-2, ypix-2);
XC3Horizontal2 = zeros(xpix-2, ypix-2);
XC3Vertical1 = zeros(xpix-2, ypix-2);
XC3Vertical2 = zeros(xpix-2, ypix-2);
XC3Diagonal1 = zeros(xpix-2, ypix-2);
XC3Diagonal2 = zeros(xpix-2, ypix-2);
XC3Diagonal3 = zeros(xpix-2, ypix-2);
XC3Diagonal4 = zeros(xpix-2, ypix-2);
for k = 1:size(tau,1)
    XC3Central = XC3Central + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,2)).*delta(i,j,t+tau(k,3))),3)/size(tau,1);
    XC3Horizontal1 = XC3Horizontal1 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,2)).*delta(i,j+1,t+tau(k,3))),3)/size(tau,1);
    XC3Horizontal2 = XC3Horizontal2 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j+1,t+tau(k,2)).*delta(i,j+1,t+tau(k,3))),3)/size(tau,1);
    XC3Vertical1 = XC3Vertical1 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,2)).*delta(i+1,j,t+tau(k,3))),3)/size(tau,1);
    XC3Vertical2 =  XC3Vertical2 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j,t+tau(k,2)).*delta(i+1,j,t+tau(k,3))),3)/size(tau,1);
    XC3Diagonal1 = XC3Diagonal1 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,2)).*delta(i+1,j+1,t+tau(k,3))),3)/size(tau,1);
    XC3Diagonal2 = XC3Diagonal2 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j+1,t+tau(k,2)).*delta(i+1,j+1,t+tau(k,3))),3)/size(tau,1);
    XC3Diagonal3 = XC3Diagonal3 + mean(abs(delta(i,j+1,t+tau(k,1)).*delta(i,j+1,t+tau(k,2)).*delta(i+1,j,t+tau(k,3))),3)/size(tau,1);
    XC3Diagonal4 = XC3Diagonal4 + mean(abs(delta(i,j+1,t+tau(k,1)).*delta(i+1,j,t+tau(k,2)).*delta(i+1,j,t+tau(k,3))),3)/size(tau,1);
end
XC3xpixOUT = 3*(xpix-2)-2;
XC3ypixOUT = 3*(ypix-2)-2;
XC3 = zeros(XC3xpixOUT, XC3ypixOUT);
XC3(1:3:end,1:3:end) = XC3Central(1:end,1:end);
XC3(1:3:end,2:3:end) = XC3Horizontal1(1:end,1:end-1);
XC3(1:3:end,3:3:end) = XC3Horizontal2(1:end,1:end-1);
XC3(2:3:end,1:3:end) = XC3Vertical1(1:end-1,1:end);
XC3(3:3:end,1:3:end) = XC3Vertical2(1:end-1,1:end);
XC3(2:3:end,2:3:end) = XC3Diagonal1(1:end-1,1:end-1);
XC3(3:3:end,3:3:end) = XC3Diagonal2(1:end-1,1:end-1);
XC3(2:3:end,3:3:end) = XC3Diagonal3(1:end-1,1:end-1);
XC3(3:3:end,2:3:end) = XC3Diagonal4(1:end-1,1:end-1);

fprintf('Finished.\n');
clear XC3Central XC3Horizontal1 XC3Horizontal2...
    XC3Vertical1 XC3Vertical2 XC3Diagonal1...
    XC3Diagonal2 XC3Diagonal3 XC3Diagonal4;


%%%%%%%%%%% 4th-order cross-correlation
fprintf('Calculating 4th-order cross-correlation SOFI image (without cumulants)...');

XG4Central = zeros(xpix-2, ypix-2);
XG4Horizontal1 = zeros(xpix-2, ypix-2);
XG4Horizontal2 = zeros(xpix-2, ypix-2);
XG4Horizontal3 = zeros(xpix-2, ypix-2);
XG4Vertical1 = zeros(xpix-2, ypix-2);
XG4Vertical2 = zeros(xpix-2, ypix-2);
XG4Vertical3 = zeros(xpix-2, ypix-2);
XG4Diagonal1 = zeros(xpix-2, ypix-2);
XG4Diagonal2 = zeros(xpix-2, ypix-2);
XG4Diagonal3 = zeros(xpix-2, ypix-2);
XG4Diagonal4 = zeros(xpix-2, ypix-2);
XG4Diagonal5 = zeros(xpix-2, ypix-2);
XG4Diagonal6 = zeros(xpix-2, ypix-2);
XG4Diagonal7 = zeros(xpix-2, ypix-2);
XG4Diagonal8 = zeros(xpix-2, ypix-2);
XG4Diagonal9 = zeros(xpix-2, ypix-2);
for k = 1:size(tau,1)
    XG4Central = XG4Central + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,2)).*delta(i,j,t+tau(k,3)).*delta(i,j+1,t+tau(k,4))),3)/size(tau,1);
    XG4Horizontal1 = XG4Horizontal1 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,2)).*delta(i,j,t+tau(k,3)).*delta(i,j+1,t+tau(k,4))),3)/size(tau,1);
    XG4Horizontal2 = XG4Horizontal2 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,2)).*delta(i,j+1,t+tau(k,3)).*delta(i,j+1,t+tau(k,4))),3)/size(tau,1);
    XG4Horizontal3 = XG4Horizontal3 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j+1,t+tau(k,2)).*delta(i,j+1,t+tau(k,3)).*delta(i,j+1,t+tau(k,4))),3)/size(tau,1);
    XG4Vertical1 = XG4Vertical1 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,2)).*delta(i,j,t+tau(k,3)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4Vertical2 = XG4Vertical2 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,2)).*delta(i+1,j,t+tau(k,3)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4Vertical3 = XG4Vertical3 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j,t+tau(k,2)).*delta(i+1,j,t+tau(k,3)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4Diagonal1 = XG4Diagonal1 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,2)).*delta(i,j+1,t+tau(k,3)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4Diagonal2 = XG4Diagonal2 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j+1,t+tau(k,2)).*delta(i,j+1,t+tau(k,3)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4Diagonal3 = XG4Diagonal3 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j+1,t+tau(k,2)).*delta(i,j+1,t+tau(k,3)).*delta(i+1,j+1,t+tau(k,4))),3)/size(tau,1);
    XG4Diagonal4 = XG4Diagonal4 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j+1,t+tau(k,2)).*delta(i+1,j,t+tau(k,3)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4Diagonal5 = XG4Diagonal5 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j+1,t+tau(k,2)).*delta(i+1,j+1,t+tau(k,3)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4Diagonal6 = XG4Diagonal6 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j+1,t+tau(k,2)).*delta(i+1,j+1,t+tau(k,3)).*delta(i+1,j+1,t+tau(k,4))),3)/size(tau,1);
    XG4Diagonal7 = XG4Diagonal7 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j+1,t+tau(k,2)).*delta(i+1,j,t+tau(k,3)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4Diagonal8 = XG4Diagonal8 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j+1,t+tau(k,2)).*delta(i+1,j+1,t+tau(k,3)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4Diagonal9 = XG4Diagonal9 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j+1,t+tau(k,2)).*delta(i+1,j+1,t+tau(k,3)).*delta(i+1,j+1,t+tau(k,4))),3)/size(tau,1);
end
XG4xpixOUT = 4*(xpix-2)-3;
XG4ypixOUT = 4*(ypix-2)-3;
XG4 = zeros(XG4xpixOUT, XG4ypixOUT);
XG4(1:4:end,1:4:end) = XG4Central(1:end,1:end);
XG4(1:4:end,2:4:end) = XG4Horizontal1(1:end,1:end-1);
XG4(1:4:end,3:4:end) = XG4Horizontal2(1:end,1:end-1);
XG4(1:4:end,4:4:end) = XG4Horizontal3(1:end,1:end-1);
XG4(2:4:end,1:4:end) = XG4Vertical1(1:end-1,1:end);
XG4(3:4:end,1:4:end) = XG4Vertical2(1:end-1,1:end);
XG4(4:4:end,1:4:end) = XG4Vertical3(1:end-1,1:end);
XG4(2:4:end,2:4:end) = XG4Diagonal1(1:end-1,1:end-1);
XG4(2:4:end,3:4:end) = XG4Diagonal2(1:end-1,1:end-1);
XG4(2:4:end,4:4:end) = XG4Diagonal3(1:end-1,1:end-1);
XG4(3:4:end,2:4:end) = XG4Diagonal4(1:end-1,1:end-1);
XG4(3:4:end,3:4:end) = XG4Diagonal5(1:end-1,1:end-1);
XG4(3:4:end,4:4:end) = XG4Diagonal6(1:end-1,1:end-1);
XG4(4:4:end,2:4:end) = XG4Diagonal7(1:end-1,1:end-1);
XG4(4:4:end,3:4:end) = XG4Diagonal8(1:end-1,1:end-1);
XG4(4:4:end,4:4:end) = XG4Diagonal9(1:end-1,1:end-1);
clear XG4Central XG4Horizontal1 XG4Horizontal2 XG4Horizontal3...
    XG4Vertical1 XG4Vertical2 XG4Vertical3...
    XG4Diagonal1 XG4Diagonal2 XG4Diagonal3...
    XG4Diagonal4 XG4Diagonal5 XG4Diagonal6...
    XG4Diagonal7 XG4Diagonal8 XG4Diagonal9;
fprintf('Finished.\n');

%%%%%%%%%%% 4th-order cross-correlation cumulant
fprintf('Calculating 4th-order cross-correlation SOFI image (with cumulants)...');
%%%%%%%%%%% XG4_12
XG4_12_Central = zeros(xpix-2, ypix-2);
XG4_12_Horizontal1 = zeros(xpix-2, ypix-2);
XG4_12_Horizontal2 = zeros(xpix-2, ypix-2);
XG4_12_Horizontal3 = zeros(xpix-2, ypix-2);
XG4_12_Vertical1 = zeros(xpix-2, ypix-2);
XG4_12_Vertical2 = zeros(xpix-2, ypix-2);
XG4_12_Vertical3 = zeros(xpix-2, ypix-2);
XG4_12_Diagonal1 = zeros(xpix-2, ypix-2);
XG4_12_Diagonal2 = zeros(xpix-2, ypix-2);
XG4_12_Diagonal3 = zeros(xpix-2, ypix-2);
XG4_12_Diagonal4 = zeros(xpix-2, ypix-2);
XG4_12_Diagonal5 = zeros(xpix-2, ypix-2);
XG4_12_Diagonal6 = zeros(xpix-2, ypix-2);
XG4_12_Diagonal7 = zeros(xpix-2, ypix-2);
XG4_12_Diagonal8 = zeros(xpix-2, ypix-2);
XG4_12_Diagonal9 = zeros(xpix-2, ypix-2);
for k = 1:size(tau,1)
    XG4_12_Central = XG4_12_Central + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,2))),3)/size(tau,1);
    XG4_12_Horizontal1 = XG4_12_Horizontal1 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,2))),3)/size(tau,1);
    XG4_12_Horizontal2 = XG4_12_Horizontal2 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,2))),3)/size(tau,1);
    XG4_12_Horizontal3 = XG4_12_Horizontal3 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j+1,t+tau(k,2))),3)/size(tau,1);
    XG4_12_Vertical1 = XG4_12_Vertical1 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,2))),3)/size(tau,1);
    XG4_12_Vertical2 = XG4_12_Vertical2 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,2))),3)/size(tau,1);
    XG4_12_Vertical3 = XG4_12_Vertical3 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j,t+tau(k,2))),3)/size(tau,1);
    XG4_12_Diagonal1 = XG4_12_Diagonal1 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,2))),3)/size(tau,1);
    XG4_12_Diagonal2 = XG4_12_Diagonal2 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j+1,t+tau(k,2))),3)/size(tau,1);
    XG4_12_Diagonal3 = XG4_12_Diagonal3 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j+1,t+tau(k,2))),3)/size(tau,1);
    XG4_12_Diagonal4 = XG4_12_Diagonal4 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j+1,t+tau(k,2))),3)/size(tau,1);
    XG4_12_Diagonal5 = XG4_12_Diagonal5 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j+1,t+tau(k,2))),3)/size(tau,1);
    XG4_12_Diagonal6 = XG4_12_Diagonal6 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j+1,t+tau(k,2))),3)/size(tau,1);
    XG4_12_Diagonal7 = XG4_12_Diagonal7 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j+1,t+tau(k,2))),3)/size(tau,1);
    XG4_12_Diagonal8 = XG4_12_Diagonal8 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j+1,t+tau(k,2))),3)/size(tau,1);
    XG4_12_Diagonal9 = XG4_12_Diagonal9 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j+1,t+tau(k,2))),3)/size(tau,1);
end
XG4_12 = zeros(XG4xpixOUT, XG4ypixOUT);
XG4_12(1:4:end,1:4:end) = XG4_12_Central(1:end,1:end);
XG4_12(1:4:end,2:4:end) = XG4_12_Horizontal1(1:end,1:end-1);
XG4_12(1:4:end,3:4:end) = XG4_12_Horizontal2(1:end,1:end-1);
XG4_12(1:4:end,4:4:end) = XG4_12_Horizontal3(1:end,1:end-1);
XG4_12(2:4:end,1:4:end) = XG4_12_Vertical1(1:end-1,1:end);
XG4_12(3:4:end,1:4:end) = XG4_12_Vertical2(1:end-1,1:end);
XG4_12(4:4:end,1:4:end) = XG4_12_Vertical3(1:end-1,1:end);
XG4_12(2:4:end,2:4:end) = XG4_12_Diagonal1(1:end-1,1:end-1);
XG4_12(2:4:end,3:4:end) = XG4_12_Diagonal2(1:end-1,1:end-1);
XG4_12(2:4:end,4:4:end) = XG4_12_Diagonal3(1:end-1,1:end-1);
XG4_12(3:4:end,2:4:end) = XG4_12_Diagonal4(1:end-1,1:end-1);
XG4_12(3:4:end,3:4:end) = XG4_12_Diagonal5(1:end-1,1:end-1);
XG4_12(3:4:end,4:4:end) = XG4_12_Diagonal6(1:end-1,1:end-1);
XG4_12(4:4:end,2:4:end) = XG4_12_Diagonal7(1:end-1,1:end-1);
XG4_12(4:4:end,3:4:end) = XG4_12_Diagonal8(1:end-1,1:end-1);
XG4_12(4:4:end,4:4:end) = XG4_12_Diagonal9(1:end-1,1:end-1);
clear XG4_12_Central XG4_12_Horizontal1 XG4_12_Horizontal2 XG4_12_Horizontal3...
    XG4_12_Vertical1 XG4_12_Vertical2 XG4_12_Vertical3...
    XG4_12_Diagonal1 XG4_12_Diagonal2 XG4_12_Diagonal3...
    XG4_12_Diagonal4 XG4_12_Diagonal5 XG4_12_Diagonal6...
    XG4_12_Diagonal7 XG4_12_Diagonal8 XG4_12_Diagonal9;
%%%%%%%%%%% XG4_13
XG4_13_Central = zeros(xpix-2, ypix-2);
XG4_13_Horizontal1 = zeros(xpix-2, ypix-2);
XG4_13_Horizontal2 = zeros(xpix-2, ypix-2);
XG4_13_Horizontal3 = zeros(xpix-2, ypix-2);
XG4_13_Vertical1 = zeros(xpix-2, ypix-2);
XG4_13_Vertical2 = zeros(xpix-2, ypix-2);
XG4_13_Vertical3 = zeros(xpix-2, ypix-2);
XG4_13_Diagonal1 = zeros(xpix-2, ypix-2);
XG4_13_Diagonal2 = zeros(xpix-2, ypix-2);
XG4_13_Diagonal3 = zeros(xpix-2, ypix-2);
XG4_13_Diagonal4 = zeros(xpix-2, ypix-2);
XG4_13_Diagonal5 = zeros(xpix-2, ypix-2);
XG4_13_Diagonal6 = zeros(xpix-2, ypix-2);
XG4_13_Diagonal7 = zeros(xpix-2, ypix-2);
XG4_13_Diagonal8 = zeros(xpix-2, ypix-2);
XG4_13_Diagonal9 = zeros(xpix-2, ypix-2);
for k = 1:size(tau,1)
    XG4_13_Central = XG4_13_Central + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,3))),3)/size(tau,1);
    XG4_13_Horizontal1 = XG4_13_Horizontal1 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,3))),3)/size(tau,1);
    XG4_13_Horizontal2 = XG4_13_Horizontal2 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j+1,t+tau(k,3))),3)/size(tau,1);
    XG4_13_Horizontal3 = XG4_13_Horizontal3 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j+1,t+tau(k,3))),3)/size(tau,1);
    XG4_13_Vertical1 = XG4_13_Vertical1 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j,t+tau(k,3))),3)/size(tau,1);
    XG4_13_Vertical2 = XG4_13_Vertical2 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j,t+tau(k,3))),3)/size(tau,1);
    XG4_13_Vertical3 = XG4_13_Vertical3 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j,t+tau(k,3))),3)/size(tau,1);
    XG4_13_Diagonal1 = XG4_13_Diagonal1 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j+1,t+tau(k,3))),3)/size(tau,1);
    XG4_13_Diagonal2 = XG4_13_Diagonal2 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j+1,t+tau(k,3))),3)/size(tau,1);
    XG4_13_Diagonal3 = XG4_13_Diagonal3 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j+1,t+tau(k,3))),3)/size(tau,1);
    XG4_13_Diagonal4 = XG4_13_Diagonal4 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j,t+tau(k,3))),3)/size(tau,1);
    XG4_13_Diagonal5 = XG4_13_Diagonal5 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j+1,t+tau(k,3))),3)/size(tau,1);
    XG4_13_Diagonal6 = XG4_13_Diagonal6 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j+1,t+tau(k,3))),3)/size(tau,1);
    XG4_13_Diagonal7 = XG4_13_Diagonal7 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j,t+tau(k,3))),3)/size(tau,1);
    XG4_13_Diagonal8 = XG4_13_Diagonal8 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j+1,t+tau(k,3))),3)/size(tau,1);
    XG4_13_Diagonal9 = XG4_13_Diagonal9 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j+1,t+tau(k,3))),3)/size(tau,1);
end
XG4_13 = zeros(XG4xpixOUT, XG4ypixOUT);
XG4_13(1:4:end,1:4:end) = XG4_13_Central(1:end,1:end);
XG4_13(1:4:end,2:4:end) = XG4_13_Horizontal1(1:end,1:end-1);
XG4_13(1:4:end,3:4:end) = XG4_13_Horizontal2(1:end,1:end-1);
XG4_13(1:4:end,4:4:end) = XG4_13_Horizontal3(1:end,1:end-1);
XG4_13(2:4:end,1:4:end) = XG4_13_Vertical1(1:end-1,1:end);
XG4_13(3:4:end,1:4:end) = XG4_13_Vertical2(1:end-1,1:end);
XG4_13(4:4:end,1:4:end) = XG4_13_Vertical3(1:end-1,1:end);
XG4_13(2:4:end,2:4:end) = XG4_13_Diagonal1(1:end-1,1:end-1);
XG4_13(2:4:end,3:4:end) = XG4_13_Diagonal2(1:end-1,1:end-1);
XG4_13(2:4:end,4:4:end) = XG4_13_Diagonal3(1:end-1,1:end-1);
XG4_13(3:4:end,2:4:end) = XG4_13_Diagonal4(1:end-1,1:end-1);
XG4_13(3:4:end,3:4:end) = XG4_13_Diagonal5(1:end-1,1:end-1);
XG4_13(3:4:end,4:4:end) = XG4_13_Diagonal6(1:end-1,1:end-1);
XG4_13(4:4:end,2:4:end) = XG4_13_Diagonal7(1:end-1,1:end-1);
XG4_13(4:4:end,3:4:end) = XG4_13_Diagonal8(1:end-1,1:end-1);
XG4_13(4:4:end,4:4:end) = XG4_13_Diagonal9(1:end-1,1:end-1);
clear XG4_13_Central XG4_13_Horizontal1 XG4_13_Horizontal2 XG4_13_Horizontal3...
    XG4_13_Vertical1 XG4_13_Vertical2 XG4_13_Vertical3...
    XG4_13_Diagonal1 XG4_13_Diagonal2 XG4_13_Diagonal3...
    XG4_13_Diagonal4 XG4_13_Diagonal5 XG4_13_Diagonal6...
    XG4_13_Diagonal7 XG4_13_Diagonal8 XG4_13_Diagonal9;
%%%%%%%%%%% XG4_14
XG4_14_Central = zeros(xpix-2, ypix-2);
XG4_14_Horizontal1 = zeros(xpix-2, ypix-2);
XG4_14_Horizontal2 = zeros(xpix-2, ypix-2);
XG4_14_Horizontal3 = zeros(xpix-2, ypix-2);
XG4_14_Vertical1 = zeros(xpix-2, ypix-2);
XG4_14_Vertical2 = zeros(xpix-2, ypix-2);
XG4_14_Vertical3 = zeros(xpix-2, ypix-2);
XG4_14_Diagonal1 = zeros(xpix-2, ypix-2);
XG4_14_Diagonal2 = zeros(xpix-2, ypix-2);
XG4_14_Diagonal3 = zeros(xpix-2, ypix-2);
XG4_14_Diagonal4 = zeros(xpix-2, ypix-2);
XG4_14_Diagonal5 = zeros(xpix-2, ypix-2);
XG4_14_Diagonal6 = zeros(xpix-2, ypix-2);
XG4_14_Diagonal7 = zeros(xpix-2, ypix-2);
XG4_14_Diagonal8 = zeros(xpix-2, ypix-2);
XG4_14_Diagonal9 = zeros(xpix-2, ypix-2);
for k = 1:size(tau,1)
    XG4_14_Central = XG4_14_Central + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j+1,t+tau(k,4))),3)/size(tau,1);
    XG4_14_Horizontal1 = XG4_14_Horizontal1 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j+1,t+tau(k,4))),3)/size(tau,1);
    XG4_14_Horizontal2 = XG4_14_Horizontal2 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j+1,t+tau(k,4))),3)/size(tau,1);
    XG4_14_Horizontal3 = XG4_14_Horizontal3 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i,j+1,t+tau(k,4))),3)/size(tau,1);
    XG4_14_Vertical1 = XG4_14_Vertical1 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_14_Vertical2 = XG4_14_Vertical2 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_14_Vertical3 = XG4_14_Vertical3 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_14_Diagonal1 = XG4_14_Diagonal1 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_14_Diagonal2 = XG4_14_Diagonal2 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_14_Diagonal3 = XG4_14_Diagonal3 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j+1,t+tau(k,4))),3)/size(tau,1);
    XG4_14_Diagonal4 = XG4_14_Diagonal4 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_14_Diagonal5 = XG4_14_Diagonal5 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_14_Diagonal6 = XG4_14_Diagonal6 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j+1,t+tau(k,4))),3)/size(tau,1);
    XG4_14_Diagonal7 = XG4_14_Diagonal7 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_14_Diagonal8 = XG4_14_Diagonal8 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_14_Diagonal9 = XG4_14_Diagonal9 + mean(abs(delta(i,j,t+tau(k,1)).*delta(i+1,j+1,t+tau(k,4))),3)/size(tau,1);
end
XG4_14 = zeros(XG4xpixOUT, XG4ypixOUT);
XG4_14(1:4:end,1:4:end) = XG4_14_Central(1:end,1:end);
XG4_14(1:4:end,2:4:end) = XG4_14_Horizontal1(1:end,1:end-1);
XG4_14(1:4:end,3:4:end) = XG4_14_Horizontal2(1:end,1:end-1);
XG4_14(1:4:end,4:4:end) = XG4_14_Horizontal3(1:end,1:end-1);
XG4_14(2:4:end,1:4:end) = XG4_14_Vertical1(1:end-1,1:end);
XG4_14(3:4:end,1:4:end) = XG4_14_Vertical2(1:end-1,1:end);
XG4_14(4:4:end,1:4:end) = XG4_14_Vertical3(1:end-1,1:end);
XG4_14(2:4:end,2:4:end) = XG4_14_Diagonal1(1:end-1,1:end-1);
XG4_14(2:4:end,3:4:end) = XG4_14_Diagonal2(1:end-1,1:end-1);
XG4_14(2:4:end,4:4:end) = XG4_14_Diagonal3(1:end-1,1:end-1);
XG4_14(3:4:end,2:4:end) = XG4_14_Diagonal4(1:end-1,1:end-1);
XG4_14(3:4:end,3:4:end) = XG4_14_Diagonal5(1:end-1,1:end-1);
XG4_14(3:4:end,4:4:end) = XG4_14_Diagonal6(1:end-1,1:end-1);
XG4_14(4:4:end,2:4:end) = XG4_14_Diagonal7(1:end-1,1:end-1);
XG4_14(4:4:end,3:4:end) = XG4_14_Diagonal8(1:end-1,1:end-1);
XG4_14(4:4:end,4:4:end) = XG4_14_Diagonal9(1:end-1,1:end-1);
clear XG4_14_Central XG4_14_Horizontal1 XG4_14_Horizontal2 XG4_14_Horizontal3...
    XG4_14_Vertical1 XG4_14_Vertical2 XG4_14_Vertical3...
    XG4_14_Diagonal1 XG4_14_Diagonal2 XG4_14_Diagonal3...
    XG4_14_Diagonal4 XG4_14_Diagonal5 XG4_14_Diagonal6...
    XG4_14_Diagonal7 XG4_14_Diagonal8 XG4_14_Diagonal9;
%%%%%%%%%%% XG4_23
XG4_23_Central = zeros(xpix-2, ypix-2);
XG4_23_Horizontal1 = zeros(xpix-2, ypix-2);
XG4_23_Horizontal2 = zeros(xpix-2, ypix-2);
XG4_23_Horizontal3 = zeros(xpix-2, ypix-2);
XG4_23_Vertical1 = zeros(xpix-2, ypix-2);
XG4_23_Vertical2 = zeros(xpix-2, ypix-2);
XG4_23_Vertical3 = zeros(xpix-2, ypix-2);
XG4_23_Diagonal1 = zeros(xpix-2, ypix-2);
XG4_23_Diagonal2 = zeros(xpix-2, ypix-2);
XG4_23_Diagonal3 = zeros(xpix-2, ypix-2);
XG4_23_Diagonal4 = zeros(xpix-2, ypix-2);
XG4_23_Diagonal5 = zeros(xpix-2, ypix-2);
XG4_23_Diagonal6 = zeros(xpix-2, ypix-2);
XG4_23_Diagonal7 = zeros(xpix-2, ypix-2);
XG4_23_Diagonal8 = zeros(xpix-2, ypix-2);
XG4_23_Diagonal9 = zeros(xpix-2, ypix-2);
for k = 1:size(tau,1)
    XG4_23_Central = XG4_23_Central + mean(abs(delta(i,j,t+tau(k,2)).*delta(i,j,t+tau(k,3))),3)/size(tau,1);
    XG4_23_Horizontal1 = XG4_23_Horizontal1 + mean(abs(delta(i,j,t+tau(k,2)).*delta(i,j,t+tau(k,3))),3)/size(tau,1);
    XG4_23_Horizontal2 = XG4_23_Horizontal2 + mean(abs(delta(i,j,t+tau(k,2)).*delta(i,j+1,t+tau(k,3))),3)/size(tau,1);
    XG4_23_Horizontal3 = XG4_23_Horizontal3 + mean(abs(delta(i,j+1,t+tau(k,2)).*delta(i,j+1,t+tau(k,3))),3)/size(tau,1);
    XG4_23_Vertical1 = XG4_23_Vertical1 + mean(abs(delta(i,j,t+tau(k,2)).*delta(i,j,t+tau(k,3))),3)/size(tau,1);
    XG4_23_Vertical2 = XG4_23_Vertical2 + mean(abs(delta(i,j,t+tau(k,2)).*delta(i+1,j,t+tau(k,3))),3)/size(tau,1);
    XG4_23_Vertical3 = XG4_23_Vertical3 + mean(abs(delta(i+1,j,t+tau(k,2)).*delta(i+1,j,t+tau(k,3))),3)/size(tau,1);
    XG4_23_Diagonal1 = XG4_23_Diagonal1 + mean(abs(delta(i,j,t+tau(k,2)).*delta(i,j+1,t+tau(k,3))),3)/size(tau,1);
    XG4_23_Diagonal2 = XG4_23_Diagonal2 + mean(abs(delta(i,j+1,t+tau(k,2)).*delta(i,j+1,t+tau(k,3))),3)/size(tau,1);
    XG4_23_Diagonal3 = XG4_23_Diagonal3 + mean(abs(delta(i,j+1,t+tau(k,2)).*delta(i,j+1,t+tau(k,3))),3)/size(tau,1);
    XG4_23_Diagonal4 = XG4_23_Diagonal4 + mean(abs(delta(i,j+1,t+tau(k,2)).*delta(i+1,j,t+tau(k,3))),3)/size(tau,1);
    XG4_23_Diagonal5 = XG4_23_Diagonal5 + mean(abs(delta(i,j+1,t+tau(k,2)).*delta(i+1,j+1,t+tau(k,3))),3)/size(tau,1);
    XG4_23_Diagonal6 = XG4_23_Diagonal6 + mean(abs(delta(i,j+1,t+tau(k,2)).*delta(i+1,j+1,t+tau(k,3))),3)/size(tau,1);
    XG4_23_Diagonal7 = XG4_23_Diagonal7 + mean(abs(delta(i+1,j+1,t+tau(k,2)).*delta(i+1,j,t+tau(k,3))),3)/size(tau,1);
    XG4_23_Diagonal8 = XG4_23_Diagonal8 + mean(abs(delta(i+1,j+1,t+tau(k,2)).*delta(i+1,j+1,t+tau(k,3))),3)/size(tau,1);
    XG4_23_Diagonal9 = XG4_23_Diagonal9 + mean(abs(delta(i+1,j+1,t+tau(k,2)).*delta(i+1,j+1,t+tau(k,3))),3)/size(tau,1);
end
XG4_23 = zeros(XG4xpixOUT, XG4ypixOUT);
XG4_23(1:4:end,1:4:end) = XG4_23_Central(1:end,1:end);
XG4_23(1:4:end,2:4:end) = XG4_23_Horizontal1(1:end,1:end-1);
XG4_23(1:4:end,3:4:end) = XG4_23_Horizontal2(1:end,1:end-1);
XG4_23(1:4:end,4:4:end) = XG4_23_Horizontal3(1:end,1:end-1);
XG4_23(2:4:end,1:4:end) = XG4_23_Vertical1(1:end-1,1:end);
XG4_23(3:4:end,1:4:end) = XG4_23_Vertical2(1:end-1,1:end);
XG4_23(4:4:end,1:4:end) = XG4_23_Vertical3(1:end-1,1:end);
XG4_23(2:4:end,2:4:end) = XG4_23_Diagonal1(1:end-1,1:end-1);
XG4_23(2:4:end,3:4:end) = XG4_23_Diagonal2(1:end-1,1:end-1);
XG4_23(2:4:end,4:4:end) = XG4_23_Diagonal3(1:end-1,1:end-1);
XG4_23(3:4:end,2:4:end) = XG4_23_Diagonal4(1:end-1,1:end-1);
XG4_23(3:4:end,3:4:end) = XG4_23_Diagonal5(1:end-1,1:end-1);
XG4_23(3:4:end,4:4:end) = XG4_23_Diagonal6(1:end-1,1:end-1);
XG4_23(4:4:end,2:4:end) = XG4_23_Diagonal7(1:end-1,1:end-1);
XG4_23(4:4:end,3:4:end) = XG4_23_Diagonal8(1:end-1,1:end-1);
XG4_23(4:4:end,4:4:end) = XG4_23_Diagonal9(1:end-1,1:end-1);
clear XG4_23_Central XG4_23_Horizontal1 XG4_23_Horizontal2 XG4_23_Horizontal3...
    XG4_23_Vertical1 XG4_23_Vertical2 XG4_23_Vertical3...
    XG4_23_Diagonal1 XG4_23_Diagonal2 XG4_23_Diagonal3...
    XG4_23_Diagonal4 XG4_23_Diagonal5 XG4_23_Diagonal6...
    XG4_23_Diagonal7 XG4_23_Diagonal8 XG4_23_Diagonal9;
%%%%%%%%%%% XG4_24
XG4_24_Central = zeros(xpix-2, ypix-2);
XG4_24_Horizontal1 = zeros(xpix-2, ypix-2);
XG4_24_Horizontal2 = zeros(xpix-2, ypix-2);
XG4_24_Horizontal3 = zeros(xpix-2, ypix-2);
XG4_24_Vertical1 = zeros(xpix-2, ypix-2);
XG4_24_Vertical2 = zeros(xpix-2, ypix-2);
XG4_24_Vertical3 = zeros(xpix-2, ypix-2);
XG4_24_Diagonal1 = zeros(xpix-2, ypix-2);
XG4_24_Diagonal2 = zeros(xpix-2, ypix-2);
XG4_24_Diagonal3 = zeros(xpix-2, ypix-2);
XG4_24_Diagonal4 = zeros(xpix-2, ypix-2);
XG4_24_Diagonal5 = zeros(xpix-2, ypix-2);
XG4_24_Diagonal6 = zeros(xpix-2, ypix-2);
XG4_24_Diagonal7 = zeros(xpix-2, ypix-2);
XG4_24_Diagonal8 = zeros(xpix-2, ypix-2);
XG4_24_Diagonal9 = zeros(xpix-2, ypix-2);
for k = 1:size(tau,1)
    XG4_24_Central = XG4_24_Central + mean(abs(delta(i,j,t+tau(k,2)).*delta(i,j+1,t+tau(k,4))),3)/size(tau,1);
    XG4_24_Horizontal1 = XG4_24_Horizontal1 + mean(abs(delta(i,j,t+tau(k,2)).*delta(i,j+1,t+tau(k,4))),3)/size(tau,1);
    XG4_24_Horizontal2 = XG4_24_Horizontal2 + mean(abs(delta(i,j,t+tau(k,2)).*delta(i,j+1,t+tau(k,4))),3)/size(tau,1);
    XG4_24_Horizontal3 = XG4_24_Horizontal3 + mean(abs(delta(i,j+1,t+tau(k,2)).*delta(i,j+1,t+tau(k,4))),3)/size(tau,1);
    XG4_24_Vertical1 = XG4_24_Vertical1 + mean(abs(delta(i,j,t+tau(k,2)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_24_Vertical2 = XG4_24_Vertical2 + mean(abs(delta(i,j,t+tau(k,2)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_24_Vertical3 = XG4_24_Vertical3 + mean(abs(delta(i+1,j,t+tau(k,2)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_24_Diagonal1 = XG4_24_Diagonal1 + mean(abs(delta(i,j,t+tau(k,2)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_24_Diagonal2 = XG4_24_Diagonal2 + mean(abs(delta(i,j+1,t+tau(k,2)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_24_Diagonal3 = XG4_24_Diagonal3 + mean(abs(delta(i,j+1,t+tau(k,2)).*delta(i+1,j+1,t+tau(k,4))),3)/size(tau,1);
    XG4_24_Diagonal4 = XG4_24_Diagonal4 + mean(abs(delta(i,j+1,t+tau(k,2)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_24_Diagonal5 = XG4_24_Diagonal5 + mean(abs(delta(i,j+1,t+tau(k,2)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_24_Diagonal6 = XG4_24_Diagonal6 + mean(abs(delta(i,j+1,t+tau(k,2)).*delta(i+1,j+1,t+tau(k,4))),3)/size(tau,1);
    XG4_24_Diagonal7 = XG4_24_Diagonal7 + mean(abs(delta(i+1,j+1,t+tau(k,2)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_24_Diagonal8 = XG4_24_Diagonal8 + mean(abs(delta(i+1,j+1,t+tau(k,2)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_24_Diagonal9 = XG4_24_Diagonal9 + mean(abs(delta(i+1,j+1,t+tau(k,2)).*delta(i+1,j+1,t+tau(k,4))),3)/size(tau,1);
end
XG4_24 = zeros(XG4xpixOUT, XG4ypixOUT);
XG4_24(1:4:end,1:4:end) = XG4_24_Central(1:end,1:end);
XG4_24(1:4:end,2:4:end) = XG4_24_Horizontal1(1:end,1:end-1);
XG4_24(1:4:end,3:4:end) = XG4_24_Horizontal2(1:end,1:end-1);
XG4_24(1:4:end,4:4:end) = XG4_24_Horizontal3(1:end,1:end-1);
XG4_24(2:4:end,1:4:end) = XG4_24_Vertical1(1:end-1,1:end);
XG4_24(3:4:end,1:4:end) = XG4_24_Vertical2(1:end-1,1:end);
XG4_24(4:4:end,1:4:end) = XG4_24_Vertical3(1:end-1,1:end);
XG4_24(2:4:end,2:4:end) = XG4_24_Diagonal1(1:end-1,1:end-1);
XG4_24(2:4:end,3:4:end) = XG4_24_Diagonal2(1:end-1,1:end-1);
XG4_24(2:4:end,4:4:end) = XG4_24_Diagonal3(1:end-1,1:end-1);
XG4_24(3:4:end,2:4:end) = XG4_24_Diagonal4(1:end-1,1:end-1);
XG4_24(3:4:end,3:4:end) = XG4_24_Diagonal5(1:end-1,1:end-1);
XG4_24(3:4:end,4:4:end) = XG4_24_Diagonal6(1:end-1,1:end-1);
XG4_24(4:4:end,2:4:end) = XG4_24_Diagonal7(1:end-1,1:end-1);
XG4_24(4:4:end,3:4:end) = XG4_24_Diagonal8(1:end-1,1:end-1);
XG4_24(4:4:end,4:4:end) = XG4_24_Diagonal9(1:end-1,1:end-1);
clear XG4_24_Central XG4_24_Horizontal1 XG4_24_Horizontal2 XG4_24_Horizontal3...
    XG4_24_Vertical1 XG4_24_Vertical2 XG4_24_Vertical3...
    XG4_24_Diagonal1 XG4_24_Diagonal2 XG4_24_Diagonal3...
    XG4_24_Diagonal4 XG4_24_Diagonal5 XG4_24_Diagonal6...
    XG4_24_Diagonal7 XG4_24_Diagonal8 XG4_24_Diagonal9;
%%%%%%%%%%% XG4_34
XG4_34_Central = zeros(xpix-2, ypix-2);
XG4_34_Horizontal1 = zeros(xpix-2, ypix-2);
XG4_34_Horizontal2 = zeros(xpix-2, ypix-2);
XG4_34_Horizontal3 = zeros(xpix-2, ypix-2);
XG4_34_Vertical1 = zeros(xpix-2, ypix-2);
XG4_34_Vertical2 = zeros(xpix-2, ypix-2);
XG4_34_Vertical3 = zeros(xpix-2, ypix-2);
XG4_34_Diagonal1 = zeros(xpix-2, ypix-2);
XG4_34_Diagonal2 = zeros(xpix-2, ypix-2);
XG4_34_Diagonal3 = zeros(xpix-2, ypix-2);
XG4_34_Diagonal4 = zeros(xpix-2, ypix-2);
XG4_34_Diagonal5 = zeros(xpix-2, ypix-2);
XG4_34_Diagonal6 = zeros(xpix-2, ypix-2);
XG4_34_Diagonal7 = zeros(xpix-2, ypix-2);
XG4_34_Diagonal8 = zeros(xpix-2, ypix-2);
XG4_34_Diagonal9 = zeros(xpix-2, ypix-2);
for k = 1:size(tau,1)
    XG4_34_Central = XG4_34_Central + mean(abs(delta(i,j,t+tau(k,3)).*delta(i,j+1,t+tau(k,4))),3)/size(tau,1);
    XG4_34_Horizontal1 = XG4_34_Horizontal1 + mean(abs(delta(i,j,t+tau(k,3)).*delta(i,j+1,t+tau(k,4))),3)/size(tau,1);
    XG4_34_Horizontal2 = XG4_34_Horizontal2 + mean(abs(delta(i,j+1,t+tau(k,3)).*delta(i,j+1,t+tau(k,4))),3)/size(tau,1);
    XG4_34_Horizontal3 = XG4_34_Horizontal3 + mean(abs(delta(i,j+1,t+tau(k,3)).*delta(i,j+1,t+tau(k,4))),3)/size(tau,1);
    XG4_34_Vertical1 = XG4_34_Vertical1 + mean(abs(delta(i,j,t+tau(k,3)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_34_Vertical2 = XG4_34_Vertical2 + mean(abs(delta(i+1,j,t+tau(k,3)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_34_Vertical3 = XG4_34_Vertical3 + mean(abs(delta(i+1,j,t+tau(k,3)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_34_Diagonal1 = XG4_34_Diagonal1 + mean(abs(delta(i,j+1,t+tau(k,3)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_34_Diagonal2 = XG4_34_Diagonal2 + mean(abs(delta(i,j+1,t+tau(k,3)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_34_Diagonal3 = XG4_34_Diagonal3 + mean(abs(delta(i,j+1,t+tau(k,3)).*delta(i+1,j+1,t+tau(k,4))),3)/size(tau,1);
    XG4_34_Diagonal4 = XG4_34_Diagonal4 + mean(abs(delta(i+1,j,t+tau(k,3)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_34_Diagonal5 = XG4_34_Diagonal5 + mean(abs(delta(i+1,j+1,t+tau(k,3)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_34_Diagonal6 = XG4_34_Diagonal6 + mean(abs(delta(i+1,j+1,t+tau(k,3)).*delta(i+1,j+1,t+tau(k,4))),3)/size(tau,1);
    XG4_34_Diagonal7 = XG4_34_Diagonal7 + mean(abs(delta(i+1,j,t+tau(k,3)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_34_Diagonal8 = XG4_34_Diagonal8 + mean(abs(delta(i+1,j+1,t+tau(k,3)).*delta(i+1,j,t+tau(k,4))),3)/size(tau,1);
    XG4_34_Diagonal9 = XG4_34_Diagonal9 + mean(abs(delta(i+1,j+1,t+tau(k,3)).*delta(i+1,j+1,t+tau(k,4))),3)/size(tau,1);
end
XG4_34 = zeros(XG4xpixOUT, XG4ypixOUT);
XG4_34(1:4:end,1:4:end) = XG4_34_Central(1:end,1:end);
XG4_34(1:4:end,2:4:end) = XG4_34_Horizontal1(1:end,1:end-1);
XG4_34(1:4:end,3:4:end) = XG4_34_Horizontal2(1:end,1:end-1);
XG4_34(1:4:end,4:4:end) = XG4_34_Horizontal3(1:end,1:end-1);
XG4_34(2:4:end,1:4:end) = XG4_34_Vertical1(1:end-1,1:end);
XG4_34(3:4:end,1:4:end) = XG4_34_Vertical2(1:end-1,1:end);
XG4_34(4:4:end,1:4:end) = XG4_34_Vertical3(1:end-1,1:end);
XG4_34(2:4:end,2:4:end) = XG4_34_Diagonal1(1:end-1,1:end-1);
XG4_34(2:4:end,3:4:end) = XG4_34_Diagonal2(1:end-1,1:end-1);
XG4_34(2:4:end,4:4:end) = XG4_34_Diagonal3(1:end-1,1:end-1);
XG4_34(3:4:end,2:4:end) = XG4_34_Diagonal4(1:end-1,1:end-1);
XG4_34(3:4:end,3:4:end) = XG4_34_Diagonal5(1:end-1,1:end-1);
XG4_34(3:4:end,4:4:end) = XG4_34_Diagonal6(1:end-1,1:end-1);
XG4_34(4:4:end,2:4:end) = XG4_34_Diagonal7(1:end-1,1:end-1);
XG4_34(4:4:end,3:4:end) = XG4_34_Diagonal8(1:end-1,1:end-1);
XG4_34(4:4:end,4:4:end) = XG4_34_Diagonal9(1:end-1,1:end-1);
clear XG4_34_Central XG4_34_Horizontal1 XG4_34_Horizontal2 XG4_34_Horizontal3...
    XG4_34_Vertical1 XG4_34_Vertical2 XG4_34_Vertical3...
    XG4_34_Diagonal1 XG4_34_Diagonal2 XG4_34_Diagonal3...
    XG4_34_Diagonal4 XG4_34_Diagonal5 XG4_34_Diagonal6...
    XG4_34_Diagonal7 XG4_34_Diagonal8 XG4_34_Diagonal9;

XC4 =  abs(XG4-XG4_12.*XG4_34-XG4_13.*XG4_24-XG4_14.*XG4_23);


fprintf('Finished. \n');
end

function [ XC2corrected, XC3corrected, XG4corrected, XC4corrected ]  = XC_pixelDistanceCorrection (XC2, XC3, XG4, XC4,  horCor, vertCor, diagCor1, diagCor2)
XC2corrected = XC2;
XC3corrected = XC3;
XG4corrected = XG4;
XC4corrected = XC4;
diagCor = (diagCor1 + diagCor2)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%
%%%%%%%%%%% Correct cross-correlation SOFI with pixel distance factors
%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf( 'Correcting cross-correlation SOFI with pixel distance factors...' );

%%%%%%%%%%% Correccion of 2nd-order cross-correlation
XC2corrected(1:2:end,2:2:end)=XC2corrected(1:2:end,2:2:end)*(horCor^(-1/4));
XC2corrected(2:2:end,1:2:end)=XC2corrected(2:2:end,1:2:end)*(vertCor^(-1/4));
XC2corrected(2:2:end,2:2:end)=XC2corrected(2:2:end,2:2:end)*(diagCor^(-1/4));
XC2corrected=imresize(imgaussfilt(imresize(XC2corrected,2),0.5),0.5);

%%%%%%%%%%% Correccion of 3rd-order cross-correlation
XC3corrected(1:3:end,2:3:end)=XC3corrected(1:3:end,2:3:end)*(horCor^(-2/6));
XC3corrected(1:3:end,3:3:end)=XC3corrected(1:3:end,3:3:end)*(horCor^(-2/6));
XC3corrected(2:3:end,1:3:end)=XC3corrected(2:3:end,1:3:end)*(vertCor^(-2/6));
XC3corrected(3:3:end,1:3:end)=XC3corrected(3:3:end,1:3:end)*(vertCor^(-2/6));
XC3corrected(2:3:end,2:3:end)=XC3corrected(2:3:end,2:3:end)*(diagCor1^(-2/6));
XC3corrected(3:3:end,3:3:end)=XC3corrected(3:3:end,3:3:end)*(diagCor1^(-2/6));
XC3corrected(2:3:end,3:3:end)=XC3corrected(2:3:end,3:3:end)*(diagCor2^(-2/6));
XC3corrected(3:3:end,2:3:end)=XC3corrected(3:3:end,2:3:end)*(diagCor2^(-2/6));
XC3corrected=imresize(imgaussfilt(imresize(XC3corrected,2),1.5),0.5);
XC3corrected(XC3corrected<0)=0;

%%%%%%%%%%% Correccion of 4th-order cross-correlation
XG4corrected(1:4:end,2:4:end)=XG4corrected(1:4:end,2:4:end)*(horCor^(-3/8));
XG4corrected(1:4:end,3:4:end)=XG4corrected(1:4:end,3:4:end)*(horCor^(-4/8));
XG4corrected(1:4:end,4:4:end)=XG4corrected(1:4:end,4:4:end)*(horCor^(-3/8));
XG4corrected(2:4:end,1:4:end)=XG4corrected(2:4:end,1:4:end)*(vertCor^(-3/8));
XG4corrected(3:4:end,1:4:end)=XG4corrected(3:4:end,1:4:end)*(vertCor^(-4/8));
XG4corrected(4:4:end,1:4:end)=XG4corrected(4:4:end,1:4:end)*(vertCor^(-3/8));
XG4corrected(2:4:end,2:4:end)=XG4corrected(2:4:end,2:4:end)*(horCor^(-2/8))*(vertCor^(-2/8))*(diagCor^(-1/8));
XG4corrected(2:4:end,3:4:end)=XG4corrected(2:4:end,3:4:end)*(horCor^(-2/8))*(vertCor^(-1/8))*(diagCor^(-2/8));
XG4corrected(2:4:end,4:4:end)=XG4corrected(2:4:end,4:4:end)*(horCor^(-2/8))*(vertCor^(-2/8))*(diagCor^(-1/8));
XG4corrected(3:4:end,2:4:end)=XG4corrected(3:4:end,2:4:end)*(horCor^(-1/8))*(vertCor^(-2/8))*(diagCor^(-2/8));
XG4corrected(3:4:end,3:4:end)=XG4corrected(3:4:end,3:4:end)*(horCor^(-2/8))*(vertCor^(-2/8))*(diagCor^(-2/8));
XG4corrected(3:4:end,4:4:end)=XG4corrected(3:4:end,4:4:end)*(horCor^(-1/8))*(vertCor^(-2/8))*(diagCor^(-2/8));
XG4corrected(4:4:end,2:4:end)=XG4corrected(4:4:end,2:4:end)*(horCor^(-2/8))*(vertCor^(-2/8))*(diagCor^(-1/8));
XG4corrected(4:4:end,3:4:end)=XG4corrected(4:4:end,3:4:end)*(horCor^(-2/8))*(vertCor^(-1/8))*(diagCor^(-2/8));
XG4corrected(4:4:end,4:4:end)=XG4corrected(4:4:end,4:4:end)*(diagCor^(-3/8));
XG4corrected=imresize(imgaussfilt(imresize(XG4corrected,2),4),0.5);
XG4corrected(XG4corrected<0)=0;

%%%%%%%%%%% Correccion of 4th-order cross-correlation cumulant
XC4corrected(1:4:end,2:4:end)=XC4corrected(1:4:end,2:4:end)*(horCor^(-3/8));
XC4corrected(1:4:end,3:4:end)=XC4corrected(1:4:end,3:4:end)*(horCor^(-4/8));
XC4corrected(1:4:end,4:4:end)=XC4corrected(1:4:end,4:4:end)*(horCor^(-3/8));
XC4corrected(2:4:end,1:4:end)=XC4corrected(2:4:end,1:4:end)*(vertCor^(-3/8));
XC4corrected(3:4:end,1:4:end)=XC4corrected(3:4:end,1:4:end)*(vertCor^(-4/8));
XC4corrected(4:4:end,1:4:end)=XC4corrected(4:4:end,1:4:end)*(vertCor^(-3/8));
XC4corrected(2:4:end,2:4:end)=XC4corrected(2:4:end,2:4:end)*(horCor^(-2/8))*(vertCor^(-2/8))*(diagCor^(-1/8));
XC4corrected(2:4:end,3:4:end)=XC4corrected(2:4:end,3:4:end)*(horCor^(-2/8))*(vertCor^(-1/8))*(diagCor^(-2/8));
XC4corrected(2:4:end,4:4:end)=XC4corrected(2:4:end,4:4:end)*(horCor^(-2/8))*(vertCor^(-2/8))*(diagCor^(-1/8));
XC4corrected(3:4:end,2:4:end)=XC4corrected(3:4:end,2:4:end)*(horCor^(-1/8))*(vertCor^(-2/8))*(diagCor^(-2/8));
XC4corrected(3:4:end,3:4:end)=XC4corrected(3:4:end,3:4:end)*(horCor^(-2/8))*(vertCor^(-2/8))*(diagCor^(-2/8));
XC4corrected(3:4:end,4:4:end)=XC4corrected(3:4:end,4:4:end)*(horCor^(-1/8))*(vertCor^(-2/8))*(diagCor^(-2/8));
XC4corrected(4:4:end,2:4:end)=XC4corrected(4:4:end,2:4:end)*(horCor^(-2/8))*(vertCor^(-2/8))*(diagCor^(-1/8));
XC4corrected(4:4:end,3:4:end)=XC4corrected(4:4:end,3:4:end)*(horCor^(-2/8))*(vertCor^(-1/8))*(diagCor^(-2/8));
XC4corrected(4:4:end,4:4:end)=XC4corrected(4:4:end,4:4:end)*(diagCor^(-3/8));
XC4corrected=imresize(imgaussfilt(imresize(XC4corrected,2),4),0.5);
XC4corrected(XC4corrected<0)=0;


fprintf('Finished. \n');
end

function [ XC2deconv, XC3deconv, XG4deconv, XC4deconv, averageDeconv ] = XC_deconvolution (XC2, XC3, XG4, XC4, average, measuredPSFstd)

fprintf( 'Calculating deconvolution of cross-correlation images...' );

psfXC2 = (fspecial('gaussian', [size(XC2,1), size(XC2,2)] , double(2*measuredPSFstd))).^2; %The factor 2 of the std deviation is due to the change in pixel size
psfXC3 = (fspecial('gaussian', [size(XC3,1), size(XC2,2)] , double(3*measuredPSFstd))).^3; %The factor 3 of the std deviation is due to the change in pixel size
psfXG4 = (fspecial('gaussian', [size(XG4,1), size(XG4,2)] , double(4*measuredPSFstd))).^4; %The factor 4 of the std deviation is due to the change in pixel size
psfXC4 = (fspecial('gaussian', [size(XC4,1), size(XC4,2)] , double(4*measuredPSFstd))).^4; %The factor 4 of the std deviation is due to the change in pixel size
psfAv  =  fspecial('gaussian', [size(average,1), size(average,2)] , double(measuredPSFstd));

XC2deconv = deconvblind(XC2, psfXC2);
XC3deconv = deconvblind(XC3, psfXC3);
XG4deconv = deconvblind(XG4, psfXG4);
XC4deconv = deconvblind(XC4, psfXC4);
averageDeconv = deconvblind(average, psfAv);

fprintf( 'Finished.\n' );

end


% Drift Correction
function [ dRegistered ] = DriftCorrection( d, roi )


% Drift analysis only in manually specified ROI
if roi==1
    dipshow(d(:,:,1),'percentile')
    p = ginput;
    close gcf
    dr = d(p(1,2):p(2,2),p(1,1):p(2,1),:);
else
    dr=d;
end



numstringlength=0;
s=[0;0];

sd = size(d);
nimages = sd(3);

d=imbin(d,0.5);
dr=imbin(dr,0.5);

dgr = gaussf(dr(:,:,1),5);
for i=2:nimages
    s = [s, findshift(dgr,gaussf(dr(:,:,i),5),'grs')];
    if mod(i,50)==0;
        numstring=num2str(i);
        for j=1:numstringlength;
            fprintf('\b')
        end
        fprintf(['Determining drift...', numstring, '/', num2str(nimages)])
        numstringlength=length(['Determining drift...', numstring, '/', num2str(nimages)]);
    end
end
clear dgr

fprintf('\n')

numstringlength=0;
for i=2:nimages
    d(:,:,i) = shift(d(:,:,i),s(:,i));
    if mod(i,50)==0;
        numstring=num2str(i);
        for j=1:numstringlength;
            fprintf('\b')
        end
        fprintf(['Correcting drift...', numstring, '/', num2str(nimages)])
        numstringlength=length(['Correcting drift...', numstring, '/', num2str(nimages)]);
    end
end

dRegistered = imbin(d,2);

fprintf('\n')

end

