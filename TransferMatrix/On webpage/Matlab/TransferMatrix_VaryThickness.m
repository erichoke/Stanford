% Copyright 2010 George F. Burkhard, Eric T. Hoke, Stanford University

%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.




% This program calculates the calculated short circuit current (assuming an
% internal quantum efficiency of 100%) for a device where the thickness 
% of one of the layers is varied. The program uses the transfer matrix
% method.

% The procedure was adapted from J. Appl. Phys Vol 86, No. 1 (1999) p.487
% and JAP 93 No. 7 p. 3693.

% George Burkhard and Eric Hoke February 5, 2010
% Modifications:
% 3/3/11 Parastic absorption (parasitic_abs) is calculated and made
% accessable outside of script. 
% 3/5/15 Improved precision of fundamental constants

function TransferMatrix_VaryThickness
%------------BEGIN USER INPUT PARAMETERS SPECITIFCATION---------------
%
lambda=350:5:800; % Wavelengths over which field patterns are calculated
stepsize = 1;   % The electric field is calculated at a latice of points (nm)
                % in the device cross section seperated by this distance

% Specify Layers in device (an arbitrary number of layers is permitted) and 
% thicknesses.
%
% Change these arrays to change the order or number of layers and/or
% thickness of layers. List the layers in the order that they appear in the
% device starting with the side the light is incident on.  THE NAMES OF THE
% LAYERS MUST CORRESPOND TO THE NAMES OF THE MATERIALS IN THE INDEX OF 
% REFRACTION LIBRARY FILE, 'Index_of_Refraction_library.xls'. The first
% layer must be the transparent substrate (glass) or 'Air' if the active
% layers are on the reflective electrode (rather than transparent electrode) side 
% of the device.  The layer thicknesses are in nanometers.

layers = {'SiO2' 'ITOsorizon' 'PEDOT' 'P3HTPCBMBlendDCB' 'Ca' 'Al'}; % Names of layers of materials starting from side light is incident from
thicknesses = [0 110 25 108 7 200];  % thickness of each corresponding layer in nm (thickness of the first layer, and layer which thickness is varied is irrelivant)

% Set plotGeneration to 'true' if you want to plot generation rate as a
% function of position in the device and output the calculated short circuit current
% under AM1.5G illumination (assuming 100% internal quantum efficiency)
% plotGeneration = true;
activeLayer = 4; % index of material layer where photocurrent is generated
LayerVaried = 4; % index of material where thickness will be varied
VaryThickness = 0:5:500; % thickness values for layer that will be varied in thickness in nm

%
%------------END USER INPUT PARAMETERS SPECIFICATION-------------------

% Load in 1sun AM 1.5 solar spectrum in mW/cm2nm
AM15_data=xlsread('AM15.xls');
AM15=interp1(AM15_data(:,1), AM15_data(:,2), lambda, 'linear', 'extrap');
    
% Load in index of refraction for each material
n = zeros(size(layers,2),size(lambda,2));
for index = 1:size(layers,2)
    n(index,:) = LoadRefrIndex(layers{index},lambda);
end

% Constants
h = 6.62606957e-34; 	% Js Planck's constant
c = 2.99792458e8;	% m/s speed of light
q = 1.60217657e-19;	% C electric charge

% Calculate Incoherent power transmission through substrate
% See Griffiths "Intro to Electrodynamics 3rd Ed. Eq. 9.86 & 9.87
T_glass=abs(4*1*n(1,:)./(1+n(1,:)).^2); 
R_glass=abs((1-n(1,:))./(1+n(1,:))).^2;

t = thicknesses;
t(1)=0;
Jsc=0*VaryThickness;

%% AT EACH LAYER THICKNESS THAT IS VARIED %%%%%%%%%%%%%%%%
for ThickInd= 1:length(VaryThickness)
    t(LayerVaried)=VaryThickness(ThickInd)
    % Calculate transfer matrices, and field at each wavelength and position

    t_cumsum=cumsum(t);
    x_pos=(stepsize/2):stepsize:sum(t); %positions to evaluate field
    %   x_mat specifies what layer number the corresponding point in x_pos is in:
    x_mat= sum(repmat(x_pos,length(t),1)>repmat(t_cumsum',1,length(x_pos)),1)+1; 
    R=lambda*0;
    E=zeros(length(x_pos),length(lambda));
    for l = 1:length(lambda)
      % Calculate transfer matrices for incoherent reflection and transmission at the first interface
        S=I_mat(n(1,l),n(2,l));
        for matindex=2:(length(t)-1)
            S=S*L_mat(n(matindex,l),t(matindex),lambda(l))*I_mat(n(matindex,l),n(matindex+1,l));
        end
        R(l)=abs(S(2,1)/S(1,1))^2; %JAP Vol 86 p.487 Eq 9 Power Reflection from layers other than substrate
        T(l)=abs(2/(1+n(1,l)))/sqrt(1-R_glass(l)*R(l)); %Transmission of field through glass substrate Griffiths 9.85 + multiple reflection geometric series

        % Calculate all other transfer matrices
        for material = 2:length(t) 
            xi=2*pi*n(material,l)/lambda(l);
            dj=t(material);
            x_indices=find(x_mat == material); %indices of points which are in the material layer considered
            x=x_pos(x_indices)-t_cumsum(material-1); %distance from interface with previous layer
            % Calculate S matrices (JAP Vol 86 p.487 Eq 12 and 13)
            S_prime=I_mat(n(1,l),n(2,l));
            for matindex=3:material
                S_prime=S_prime*L_mat(n(matindex-1,l),t(matindex-1),lambda(l))*I_mat(n(matindex-1,l),n(matindex,l));
            end
            S_doubleprime=eye(2);
            for matindex=material:(length(t)-1)
                S_doubleprime=S_doubleprime*I_mat(n(matindex,l),n(matindex+1,l))*L_mat(n(matindex+1,l),t(matindex+1),lambda(l));
            end
            % Normalized Field profile (JAP Vol 86 p.487 Eq 22)
            E(x_indices,l)=T(l)*(S_doubleprime(1,1)*exp(-1i*xi*(dj-x))+S_doubleprime(2,1)*exp(1i*xi*(dj-x))) ./(S_prime(1,1)*S_doubleprime(1,1)*exp(-1i*xi*dj)+S_prime(1,2)*S_doubleprime(2,1)*exp(1i*xi*dj));
        end 
    end

    % Absorption coefficient in cm^-1 (JAP Vol 86 p.487 Eq 23)
    a=zeros(length(t),length(lambda));
    for matindex=2:length(t)
        a(matindex,:)=4*pi*imag(n(matindex,:))./(lambda*1e-7);
    end

    % Calculate 100% IQE current achieved with device

    % Energy dissipation mW/cm3-nm at each position and wavelength (JAP Vol 86 p.487 Eq 22)
    ActivePos=find(x_mat == activeLayer);
    Q=repmat(a(activeLayer,:).*real(n(activeLayer,:)).*AM15,length(ActivePos),1).*(abs(E(ActivePos,:)).^2);

    % Exciton generation rate per second-cm3-nm at each position and wavelength
    Gxl=(Q*1e-3).*repmat(lambda*1e-9,length(ActivePos),1)/(h*c);
    if length(lambda)==1
        lambdastep= 1;
    else
        lambdastep=(max(lambda)-min(lambda))/(length(lambda)-1);
    end
    Gx=sum(Gxl,2)*lambdastep; % Exciton generation rate as a function of position/(sec-cm^3)

    % outputs predicted Jsc under AM1.5 illumination assuming 100% internal
    % quantum efficiency at all wavelengths
    Jsc(ThickInd)=sum(Gx)*stepsize*1e-7*q*1e3; %in mA/cm^2
end
plot(VaryThickness, Jsc)
title('Current Density obtained from 100% IQE')
xlabel('Layer thickness (nm)')
ylabel('Current Density (mA/cm^2)')
assignin('base','Thickness',VaryThickness);
assignin('base','Jsc',Jsc);

%------------------- Helper Functions ------------------------------------
% Function I_mat
% This function calculates the transfer matrix, I, for reflection and
% transmission at an interface between materials with complex dielectric 
% constant n1 and n2.
function I = I_mat(n1,n2)
r=(n1-n2)/(n1+n2);
t=2*n1/(n1+n2);
I=[1 r; r 1]/t;

% Function L_mat
% This function calculates the propagation matrix, L, through a material of
% complex dielectric constant n and thickness d for the wavelength lambda.
function L = L_mat(n,d,lambda)
xi=2*pi*n/lambda;
L=[exp(-1i*xi*d) 0; 0 exp(1i*xi*d)];

% Function LoadRefrIndex
% This function returns the complex index of refraction spectra, ntotal, for the
% material called 'name' for each wavelength value in the wavelength vector
% 'wavelengths'.  The material must be present in the index of refraction
% library 'Index_of_Refraction_library.xls'.  The program uses linear
% interpolation/extrapolation to determine the index of refraction for
% wavelengths not listed in the library.
function ntotal = LoadRefrIndex(name,wavelengths)

%Data in IndRefr, Column names in IndRefr_names
[IndRefr,IndRefr_names]=xlsread('Index_of_Refraction_library.xls');

% Load index of refraction data in spread sheet, will crash if misspelled
file_wavelengths=IndRefr(:,strmatch('Wavelength',IndRefr_names));
n=IndRefr(:,strmatch(strcat(name,'_n'),IndRefr_names));
k=IndRefr(:,strmatch(strcat(name,'_k'),IndRefr_names));

% Interpolate/Extrapolate data linearly to desired wavelengths
n_interp=interp1(file_wavelengths, n, wavelengths, 'linear', 'extrap');
k_interp=interp1(file_wavelengths, k, wavelengths, 'linear', 'extrap');

%Return interpolated complex index of refraction data
ntotal = n_interp+1i*k_interp; 
