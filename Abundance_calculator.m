%% Calculate the concentration in the "layer" 
%{

Code to calculate abundances in layer of model (two layer system) assuming
some bulk Mantle composition (described in Sramek et al. 2013). 

Written by Scott Wipperfurth
University of Maryland

Last updated: July 25, 2018






References in Code: 
Arevalo, R., McDonough, W.F., Luong, M., 2009. The K/U ratio of the silicate
    Earth: Insights into mantle composition, structure and thermal evolution. 
    Earth and Planetary Science Letters 278, 361–369.
    https://doi.org/10.1016/j.epsl.2008.12.023

Huang, Y., Chubakov, V., Mantovani, F., Rudnick, R.L., McDonough, W.F., 
    2013. A reference Earth model for the heat-producing elements and 
    associated geoneutrino flux. Geochemistry, Geophysics, Geosystems 14,
    2003–2029. https://doi.org/10.1002/ggge.20129

Javoy, M., Kaminski, E., Guyot, F., Andrault, D., Sanloup, C., Moreira,
    M., Labrosse, S., Jambon, A., Agrinier, P., Davaille, A., Jaupart, C.,
    2010. The chemical composition of the Earth: Enstatite chondrite models.
    Earth and Planetary Science Letters 293, 259–268.
    https://doi.org/10.1016/j.epsl.2010.02.033

McDonough, W.F., Sun, S. -s., 1995. The composition of the Earth. Chem. 
    Geol. 120, 223–253. https://doi.org/10.1016/0009-2541(94)00140-4

Šrámek, O., McDonough, W.F., Kite, E.S., Leki?, V., Dye, S.T., Zhong, S., 
    2013. Geophysical and geochemical constraints on geoneutrino fluxes 
    from Earth’s mantle. Earth and Planetary Science Letters 361, 356–366. 
    https://doi.org/10.1016/j.epsl.2012.11.001

Wipperfurth, S.A., Guo, M., Šrámek, O., McDonough, W.F., 2018. Earth’s 
    chondritic Th/U: Negligible fractionation during accretion, core 
    formation, and crust–mantle differentiation. Earth and Planetary Science
    Letters 498, 196–202. https://doi.org/10.1016/j.epsl.2018.06.029


%}
%% ---- Define parameters ----
% -- Move to folder with file if not already in it --
cd(fileparts(matlab.desktop.editor.getActiveFilename));

clear; clc; close all;

% -- Load Preliminary Reference Earth Model (PREM) -- Dziewonski and Anderson (1981)
% from http://ds.iris.edu/spud/earthmodel/9991844
% Format: radius depth density 
PREM = csvread('PREM_1s.csv'); 

% -- Redefine PREM in 1km bins --
interval = 1; %1 should be good but you can change 
PREM_2 = max(PREM(:,1)):-interval:1; PREM_2 = PREM_2'; %radius (km)
index = knnsearch(PREM(:,1),PREM_2); 
PREM_2(:,2) = flipud(PREM_2); % depth (km)
PREM_2(:,3:10) = PREM(index,3:10); 


% -- Define indexes for each part of earth in PREM_2 -- (large density jumps)
jump.moho = 3.38068;  %density jump from crust to mantle
jump.cmb = 10.02942; %density jump from mantle to core
idx.crust = 1:find(PREM_2(:,3)==jump.moho)-1; 
idx.mantle = find(PREM_2(:,3)==jump.moho):find(PREM_2(:,3)==jump.cmb)-1; 
idx.core = find(PREM_2(:,3)==jump.cmb):length(PREM_2); 


% -- Mass of Layers (kg) -- (calculated from PREM, see below)
mass.mantle = 4.0347*10^24; 
mass.crust = 3.2372*10^22; 
mass.bse = 4.0671*10^24; 


% -- Define thickness of model layer (km) -- 
layer.thick = (1:interval:(max(PREM_2(idx.mantle,1)-min(PREM_2(idx.mantle,1)))))'; 


% -- Define Uranium and Potassium mass ratio -- (from NIST values)
mass_ratio.U = 1/139.567; % U235/U mass ratio (not molar!)
mass_ratio.K = 0.0001196; % K40/K mass ratio (not molar!)

% -- Define K/U and Th/U -- 
abund_ratio.KU = 13800; % Arevalo et al. (2009)
abund_ratio.ThU = 3.9; % Wipperfurth et al. (2018)

% -- Define Heat production --
hp.U235 = 5.68402e-4; %W/kg-isotope
hp.U238 = 9.4946e-5; 
hp.Th232 = 2.6368e-5;
hp.K40 = 2.8761e-5; 

clear jump
%% -- Define Bulk Mantle compositions --
%  ------------ Ambient Mantle -------------s
% Ambient mantle = Arevalo et al. (2013) (kg/kg)
ambient.U = 8.3*10^-9; 
ambient.U238 = ambient.U * (1 - mass_ratio.U); 
ambient.U235 = ambient.U * mass_ratio.U; 
ambient.Th232 = 24*10^-9; 
ambient.K = 110*10^-6; 
ambient.K40 = ambient.K * mass_ratio.K; 


%  ------------ Bulk Mantle -------------
% low endmember, middle value, and high endmember
% Mantle values for high+low from Sramek et al. (2013)

% High Endmember = McDonough and Sun (1995) (kg/kg)
bm.high.U = 12*10^-9; %kg/kg
bm.high.U238 = bm.high.U*(1-mass_ratio.U); 
bm.high.U235 = bm.high.U*(mass_ratio.U); 
bm.high.Th232 = 46*10^-9; 
bm.high.K = 192*10^-6; 
bm.high.K40 = bm.high.K *mass_ratio.K; 


% Low Endmember = Javoy (2005) (kg/kg)
bm.low.U = 4.1*10^-9; 
bm.low.U238 = bm.low.U*(1-mass_ratio.U); 
bm.low.U235 = bm.low.U*(mass_ratio.U);
bm.low.Th232 = 8.4*10^-9; 
bm.low.K = 57*10^-6; 
bm.low.K40 = bm.low.K * mass_ratio.K; 


% Medium value = 200 ppm K in Bulk Earth 
% --Arbitrary ppm K between high and low models (bulk mantle K should be ~
% 125 ppm). Assuming bulk Earth K we solve for bulk Mantle K. 
mass_K.crust = 1.52*10^-2 * mass.crust; %Huang et al. 2013
mass_K.bse = 200*10^-6 * mass.bse; 
mass_K.mantle = mass_K.bse - mass_K.crust; 

bm.med.K = mass_K.mantle/mass.mantle; %kg/kg 
bm.med.K40 = bm.med.K * mass_ratio.K; 
bm.med.U = bm.med.K / abund_ratio.KU; 
bm.med.U238 = bm.med.U * (1-mass_ratio.U); 
bm.med.U235 = bm.med.U * mass_ratio.U; 
bm.med.Th232 = bm.med.U * abund_ratio.ThU; 


%% -- Mass balance with set layer thickness --

% -- Calculate thickness and mass of each PREM_2 layer --
mass_balance.volume = 4/3*pi*(PREM_2(:,1)*1000).^3; %m^3
%mass_balance.propVol = mass_balance.volume(idx.mantle)./sum(mass_balance.volume(idx.mantle));
mass_balance.mass.mass = abs(diff(mass_balance.volume)) .* PREM_2(2:end,3)*1000; %kg


% -- Calculate mass for each possible layer thickness --
x = length(layer.thick); %thickness of mantle 
for i = 1:length(layer.thick)-1
mass_balance.mass.layer(i,1) = sum(mass_balance.mass.mass((x:-1:(x-layer.thick(i))))); 
mass_balance.propVol(i,1) = sum(mass_balance.volume((x:-1:(x-layer.thick(i)))));
end


% -- Calculate mass proportion of mantle for each possible layer thickness --
layer.prop = mass_balance.mass.layer./sum(mass_balance.mass.mass(idx.mantle)); 
layer.inv_prop = 1-layer.prop; 
layer.propVol = mass_balance.volume(idx.mantle)./sum(mass_balance.volume(idx.mantle)); 


% -- Calculate chemistry of our layer for each possible layer thickness --% (ppm)

n = {'U','U238','U235','Th232','K','K40'}; %loop through isotopes
for i = 1:length(n)
    for j = 1:length(layer.thick)-1
    chem(j,1).(n{i}) =  ((bm.high(:).(n{i}) - (ambient(:).(n{i}) * layer.inv_prop(j)))/layer.prop(j))*10^9; 
    chem(j,2).(n{i}) =  ((bm.med(:).(n{i}) - (ambient(:).(n{i}) * layer.inv_prop(j)))/layer.prop(j))*10^9; 
    chem(j,3).(n{i}) =  ((bm.low(:).(n{i}) - (ambient(:).(n{i}) * layer.inv_prop(j)))/layer.prop(j))*10^9;    
    
    end
end

% -- Calculate heat production (W/kg-rock) --
x = 10^-9; 
for n = 1:3
    for j = 1:length(layer.thick)-1
    chem(j,n).hp = chem(j,n).U238*x*hp.U238 + chem(j,n).U235*x*hp.U235...
            + chem(j,n).Th232*x*hp.Th232 + chem(j,n).K40*x*hp.K40;  
    
    chem(j,n).massProp = layer.prop(j); 
    chem(j,n).volProp = layer.propVol(j); 
    end
end

clear n j i x index 

%% -- Mass Balance with set mass % of mantle --
% this calculation assumes a spherical shell completely surrounding the
% core

layer.percent.values= 5:10; % percent of mantle that is layer

for i = 1:length(layer.percent.values)
    x = knnsearch(layer.prop, layer.percent.values(i)/100); 
    layer.percent.thick(i) = layer.thick(x); %thickness of layer with some percent mass
end


return

%% -- Calculate chemical values back through time --

yr = 1:1:4543; % Ma

n = {'U','U238','U235','Th232','K','K40'}; %loop through isotopes

for i = 1:length(n)
    for z = 1:length(yr)
        for j = 1:length(layer.thick)-1
    
            chemYr(j,1,z) = 5;% high Q model
              
        
        end 
    end
end















return
%% Plots

for i = 1:length(chem)
    x(i) = chem(i,3).K40; 
end




plot(layer.thick(1:end-1),layer.prop)
xlabel('Layer thickness (km)')
ylabel('Proportion of mantle')

% - 
figure
plot(layer.thick(1:end-1),x); hold on
ylabel('Conc. of K40 (ppm)')
axis([-inf inf 0 10])
xlabel('Layer thickness (km)')
title('Thickness vs aK40 (low BSE)')


% - 
y = x>0; 
figure
plot(layer.thick(y),layer.prop(y))
title('Layer thickness/proportion with aK40 >0 (low BSE)')
ylabel('Mass proportion of layer')
xlabel('Layer thickness (km)')














