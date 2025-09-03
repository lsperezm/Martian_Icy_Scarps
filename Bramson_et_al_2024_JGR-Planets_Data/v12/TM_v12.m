%% This code calculates temperatures, equilibrium depths and sublimation of
% ice given the parameters listed in a parameter file required to load in

% The base of this code is version 10zA6. To see the history of this code
% before this point, see that code. Enough changes occuring here that
% starting fresh.

% As of 04/04/2018 updates to use flat terrain values

% Using January 2021 for protrusion and trough wall tests...small
% modifications in directory structure but thats it

%% 1. Read in parameter files and set up a bunch of variables and constants
clear;
% Read in parameter file:
LayerTestName = 'Test1';
RunFolderName = 'Scarps/South1';

%%
FoldersPath = strcat(RunFolderName, '/', LayerTestName, '/ ');
paramFileName = strcat('paramFile_v12_', LayerTestName, '.txt');

paramFile = strcat(FoldersPath, paramFileName);
s = readStruct(paramFile); %creates a structure array by reading data of the file

%%

% Grab some important values from the structure
f = s.impexp;  % Set variable for how explicit vs implicit to numerically integrate
dt = s.dt; 
Tref = s.initSurfTemp; % First orbital timestep will use this reference temp, after that will use mean surface temperature from orbital run immediately before

% Load Laskar 2004 solutions (backwards so it is marching forward in time)
load Laskar2004Solutions_backwards2.mat;
L04Lsperi = L04Lsp - pi; % Laskar calculates at aphelion
% Number of loops it will do in total
numLaskarRuns = ceil((s.Laskar_end+1-s.Laskar_begin)/s.Laskar_step);
eqDepths = zeros(1, numLaskarRuns);
dz_excessIceTimestep = zeros(1, numLaskarRuns);

flatVisSaved = 0; % Will get written over, just to have something to pass into calcOrbitalParams in the initial or flat cases where we haven't calculated a flatVisSaved or flatIRSaved
flatIRSaved = 0; % Will get written over, just to have something to pass into calcOrbitalParams in the initial or flat cases where we haven't calculated a flatVisSaved or flatIRSaved

%% 2. Setting variables
sigmasb = 5.6705119e-8;          % Stefan Boltzmann constant W m^-2 K^-4
avogadro = 6.022140857e23; % molecules/mole
molarMassH2O = 0.01801528; % kg/mole
boltzmann = 1.38064852e-23; % Boltzmann Constant, m2 kg s-2 K-1
densityIce = s.density_ice; % 917 is non-porous, from Dundas et al. 2015, or just use value in param file
PtoRhoFactor = (molarMassH2O / avogadro) / boltzmann; % to make go faster
SF12_a1 = -1.27; % Parameterization from Schorghofer and Forget (2012)
SF12_b1 = 0.139; % Parameterization  from Schorghofer and Forget (2012)
SF12_c1 = -0.00389; % Parameterization from Schorghofer and Forget (2012)
porosity_param = (s.iceDust)*(1/(1-s.icePorosity-s.iceDust))/(1-s.regolithPorosity);

%% 3. Set up file and folder for storing output
% Will store outputs in folders organized by the date it was ran within an
% Output folder. Each file that is saved will contain the same string in
% its filename with date and time it was run
cd(FoldersPath)
if exist('Output', 'dir') == 0
    mkdir('Output');
end
cd Output
if exist(datestr(now, 'yyyy-mm-dd'), 'dir') == 0
    mkdir(datestr(now, 'yyyy-mm-dd'));
end
cd ../../../..

strFileLocationParam = strcat(FoldersPath, 'Output/', datestr(now, 'yyyy-mm-dd'), '/', datestr(now, 'yyyy-mm-dd_HHMM'));

if s.saveParamFile == 'y'
    copyfile(paramFile, strcat(strFileLocationParam, '_ParamFile.txt'));
end

fileName = strcat(strFileLocationParam, s.outputSaveFile);   % if want to input file name manually, use something like 'Folder/file.txt'
fileID = fopen(fileName,'w');
% Write a header to the output file with each value that will be saved
% later in a tab delineated text file
%fprintf(fileID, ['LaskarTimestep\t' 'time_years\t' 'ecc\t' 'obl\t' 'Lsp\t' 'Tsurf\t' 'Ticesheet\t' 'atm_v\t' 'icetable_v\t' 'iter counter\t' 'equilD\t' 'zeqConverge\t' 'PFice\t' 'track1\t' 'track2\t' 'track3\t' 'changeLag\t' 'lagThickness\t' 'dzExcessIce\t' 'totalExcessIceLoss\n']);
%fmt = '%6.0f\t%12.0f\t%12.9f\t%12.9f\t%12.9f\t%12.9f\t%12.9f\t%12.9f\t%12.9f\t%2.0f\t%12.9f\t%12.9f\t%1.0f\t%12.9f\t%12.9f\t%1.0f\t%12.9f\t%12.9f\t%12.9f\t%20.19f\n';
fprintf(fileID, ['LaskarTimestep\t' 'time_years\t' 'ecc\t' 'obl\t' 'Lsp\t' 'Tsurf\t' 'SF12\t' 'ForcedConv\t' 'FreeConv\t' 'TotalSubl\n']);
fmt = '%6.0f\t%12.0f\t%12.9f\t%12.9f\t%12.9f\t%12.9f\t%12.9f\t%20.9f\t%20.9f\t%20.9f\n';

fclose(fileID);

%% 4. Get the numerical grid of layers set up
% Calculate thermophysical properties for dry lag, pf, and excess ice
% This is run just once at the beginning to get vectors of k_input,
% rhoc_input, diurnalSkinDepths, annualSkinDepths and allKappas
layerProperties_setup; % This just run once at beginning

% Prepare vectors for flat cases (generally excess ice or whatever
% stratigraphy is input
k_flat = [s.flat_TI_input(1)*s.flat_TI_input(1)/s.rhoc_flat_input(1); s.flat_TI_input(2)*s.flat_TI_input(2)/s.rhoc_flat_input(2); s.flat_TI_input(3)*s.flat_TI_input(3)/s.rhoc_flat_input(3)];
rhoc_flat = [s.rhoc_flat_input(1); s.rhoc_flat_input(2); s.rhoc_flat_input(3)];
allKappas_flat = k_flat ./ rhoc_flat;
diurnalSkinDepths_flat = sqrt(allKappas.*lengthOfDay/pi); % meters
annualSkinDepths_flat = sqrt(allKappas.*lengthOfYear/pi); % meters

% Get initial lag thickness - assumes no pore filling ice at start
% Either use a user-specific value or performs initial run to determine the
% equilibrium depth of ice for the first timestep and start with a dry lag
% of that thickness on a flat terrain with main albedo
if s.useInputDepth == 'n'
    z_ei = annualSkinDepths(1)*annualLayers; % look within 6 (or whatever annualLayers is) skin depths of dry lag and have no pore filling ice
    z_pf = annualSkinDepths(1)*annualLayers; % no pore-filling ice- set pf value to be same as z_ei
    [k, rhoc, Kappa, dz, z_ei_index, z_pf_index, depthsAtLayerBoundaries, depthsAtMiddleOfLayers, nLayers] = setLayerProperties13(z_ei, z_ei, k_input, rhoc_input, diurnalSkinDepths, annualSkinDepths, s.layerGrowth, s.dailyLayers, s.annualLayers);

    ecc = L04ecc(s.Laskar_begin);
    Lsp = L04Lsperi(s.Laskar_begin);
    obl = L04obl(s.Laskar_begin);
    
    oblDeg = rad2deg(obl);
    if oblDeg > 10 && oblDeg < 28
        atmosPressSF12 = exp(SF12_a1 + SF12_b1*(oblDeg-28));
    elseif oblDeg > 28 && oblDeg < 50
        atmosPressSF12 = exp(SF12_a1 + SF12_b1*(oblDeg-28) + SF12_c1*(oblDeg-28)^2);
    end
    atmosPressSF12 = atmosPressSF12 * s.atmosFactor;
    
    slope = 0;
    slope_aspect = 0;
    albedo = s.albedo;
    
    findEquilibriumDepth;
    
    if PFICE == 0
        % case where no equilibrium depth found within several annual skin
        % depths
        fprintf('No equilibrium depth found. Ice not stable at all. Setting excess ice interface at 6 meters depth.');
        z_ei = 6;
        z_pf = 6;
    else 
        z_ei = z_eq; % Will set initial thickness of dry lag to be equilibrium depth
        z_pf = z_eq; % No pore-filling ice to start, excess ice starts at equilibrium depth, z_pf >= z_ei is that condition
    end
    
elseif s.useInputDepth == 'y'
    z_ei = s.initialLagThickness;
    z_pf = s.initialLagThickness;
else
    fprintf('Error in parameter file specifying which depth to use');
end

% This sets up the grids, will be run each time needs to change layer thicknesses
[k, rhoc, Kappa, dz, z_ei_index, z_pf_index, depthsAtLayerBoundaries, depthsAtMiddleOfLayers, nLayers] = setLayerProperties13(z_ei, z_pf, k_input, rhoc_input, diurnalSkinDepths, annualSkinDepths, s.layerGrowth, s.dailyLayers, s.annualLayers);

%% 5. Run model throughout orbital timesteps in Laskar's solutions
counterOrbital = 1;
for L04Step = s.Laskar_begin:s.Laskar_step:s.Laskar_end
    
    tic
    ecc = L04ecc(L04Step);
    Lsp = L04Lsperi(L04Step);
    obl = L04obl(L04Step);
    timebeforepresent = L04time(L04Step);
    fprintf('\nTesting %12.0f years ago when obliquity was %12.9f, ecc was %12.9f and Lsp was %12.9f.\n', timebeforepresent, rad2deg(obl), ecc, rad2deg(Lsp));
    
    % Calculate atmospheric water vapor density based on Schorghofer and
    % Forget 2012 scheme, which just uses obliquity. Get the atmospheric
    % partial pressure here and then convert to density once mean surface
    % temperature is calculated
    oblDeg = rad2deg(obl);
    if oblDeg > 10 && oblDeg < 28
        atmosPressSF12 = exp(SF12_a1 + SF12_b1*(oblDeg-28));
    elseif oblDeg > 28 && oblDeg < 50
        atmosPressSF12 = exp(SF12_a1 + SF12_b1*(oblDeg-28) + SF12_c1*(oblDeg-28)^2);
    end
    atmosPressSF12 = atmosPressSF12 * s.atmosFactor;
    
    %% 6a. Calculate flat terrain first if a sloped case
    % It will output values needed to flatIRSaved and flatVisSaved
    if s.slope ~= 0 && s.slope_rerad == 1
        
        % For v12, since I know I'm using this for the NPLD, will make a
        % purely icy stratigraphy for calculating surroundings
        z_ei_old = z_ei;
        z_pf_old = z_pf;
        z_ei = 0;
        z_pf = 0;
        
        % Makes sure input is just excess ice values
        %{
        k_ei = [k_input(3); k_input(3); k_input(3)];
        rhoc_ei = [rhoc_input(3); rhoc_input(3); rhoc_input(3)];
        diurnalSkinDepths_ei = [diurnalSkinDepths(3); diurnalSkinDepths(3); diurnalSkinDepths(3)];
        annualSkinDepths_ei = [annualSkinDepths(3); annualSkinDepths(3); annualSkinDepths(3)]; 
        %}
        [k, rhoc, Kappa, dz, z_ei_index, z_pf_index, depthsAtLayerBoundaries, depthsAtMiddleOfLayers, nLayers] = setLayerProperties13(z_ei, z_pf, k_flat, rhoc_flat, diurnalSkinDepths_flat, annualSkinDepths_flat, s.layerGrowth, s.dailyLayers, s.annualLayers);
        
        slope = 0;
        slope_aspect = 0;
        albedo = s.flatAlbedo;
        
        % Calculate orbital parameters and insolation values
        [sf, IRdown, visScattered, flatVis, flatIR, sky, CO2FrostPoints, nStepsInYear, lsWrapped, hr, ltst] = calcOrbitalParams(slope, slope_aspect, ecc, obl, Lsp, s, dt, flatVisSaved, flatIRSaved);

        % Run tempCalculator
        [Temps, Tsurf, frostMasses, flatIRSaved, flatVisSaved] = tempCalculator(nLayers, nStepsInYear, s, dt, k, dz, rhoc, f, Kappa, sf, visScattered, flatVis, albedo, IRdown,flatIR, sky, CO2FrostPoints, depthsAtMiddleOfLayers, sigmasb, slope, Tref);
        
        reorgDiurnal;
        fprintf('step 6a');
        z_ei = z_ei_old;
        z_pf = z_pf_old;
    end
        
    %% 6b. Calculate equilibrium depth, set stratigraphy and calculate temperatures
    if z_ei > 0
        % There is a lag deposit, so calculate these things while solving
        % for the equilibrium depth of pore filling ice in the
        % findEquilibriumDepth script but first set up layers
        [k, rhoc, Kappa, dz, z_ei_index, z_pf_index, depthsAtLayerBoundaries, depthsAtMiddleOfLayers, nLayers] = ...
            setLayerProperties13(z_ei, z_pf, k_input, rhoc_input, diurnalSkinDepths, annualSkinDepths, s.layerGrowth, s.dailyLayers, s.annualLayers);
     
        slope = s.slope;
        slope_aspect = s.slope_aspect;
        albedo = s.albedo;
        
        findEquilibriumDepth
        
        reorgDiurnal;
    else
        % No lag - only have excess ice all the way to the surface
        z_ei = 0;
        z_pf = 0;
        % Makes sure input is just excess ice values
        %{
        k_ei = [k_input(3); k_input(3); k_input(3)];
        rhoc_ei = [rhoc_input(3); rhoc_input(3); rhoc_input(3)];
        diurnalSkinDepths_ei = [diurnalSkinDepths(3); diurnalSkinDepths(3); diurnalSkinDepths(3)];
        annualSkinDepths_ei = [annualSkinDepths(3); annualSkinDepths(3); annualSkinDepths(3)]; 
        %}
        [k, rhoc, Kappa, dz, z_ei_index, z_pf_index, depthsAtLayerBoundaries, depthsAtMiddleOfLayers, nLayers] = setLayerProperties13(z_ei, z_pf, k_input, rhoc_input, diurnalSkinDepths, annualSkinDepths, s.layerGrowth, s.dailyLayers, s.annualLayers);
        
        % Calculate temperatures!
        slope = s.slope;
        slope_aspect = s.slope_aspect;
        albedo = s.albedo;
        
        % Calculate orbital parameters and insolation values
        [sf, IRdown, visScattered, flatVis, flatIR, sky, CO2FrostPoints, nStepsInYear, lsWrapped, hr, ltst] = calcOrbitalParams(slope, slope_aspect, ecc, obl, Lsp, s, dt, flatVisSaved, flatIRSaved);

        % Run tempCalculator
        [Temps, Tsurf, frostMasses, flatIRSaved, flatVisSaved] = tempCalculator(nLayers, nStepsInYear, s, dt, k, dz, rhoc, f, Kappa, sf, visScattered, flatVisSaved, albedo, IRdown, flatIRSaved, sky, CO2FrostPoints, depthsAtMiddleOfLayers, sigmasb, slope, Tref);

        reorgDiurnal;
    end
    

    %% 7. Calculate ice retreat in Mars Years!
    [TotalFree, TotalForced, TotalSublimation] = ...
        sublCalculator(L04Step, s, Ls0, L04obl, atmosPressSF12, SLOPEDaverageDiurnalSurfTemps, REGminDiurnalSurfTemps, REGdiurnalTsurfCurves, dt, densityIce, avogadro, s.windSpeed);
    
    % Convert to annual ice loss in Earth Years! (what Laskar timesteps are
    % in and diffusive ice loss calculations)
    au  = 1.4959787061e11;           % Size of one astronomical unit in meters
    gc  = 6.67259e-11;               % Gravitational constant
    sm  = 1.9891e30;                 % Solar mass
    MarsYearLength = 2*pi/sqrt(gc*sm/(s.a*au)^3);   % length of Mars year in seconds using Kepler's 3rd law
    EarthYearLength = 2*pi*sqrt(au^3/(gc*sm));   % Length of one Earth year in seconds
    
    TotalFree = TotalFree*EarthYearLength/MarsYearLength;
    TotalForced = TotalForced*EarthYearLength/MarsYearLength;
    TotalSublimation = TotalSublimation*EarthYearLength/MarsYearLength;
    
    %% 8 Set up any other variables that we'll save
    meanTsurf = mean(Tsurf); % mean surface temperature
    
    %% 9. Save outputs
    %{
    fileID = fopen(fileName,'a');
    fprintf(fileID, fmt, [L04Step timebeforepresent ecc oblDeg rad2deg(Lsp) Tsurf_mean iceTableTemps_mean atmosDensitySF12 iceTable_rhov_mean groundicecounter-1 eqDepths(counterOrbital) z_eqConverge PFICE tempTracker1 tempTracker2 tempTracker3 delta_lagThickness old_lagThickness dz_excessIceTimestep(counterOrbital) dz_excessIceCumulate(counterOrbital)]);
    fclose(fileID);
    %}
    fileID = fopen(fileName,'a');
    fprintf(fileID, fmt, [L04Step timebeforepresent ecc oblDeg rad2deg(Lsp) meanTsurf atmosPressSF12 TotalForced TotalFree TotalSublimation]);
    fclose(fileID);
    
    fprintf('\nYoure %5.0f / %5.0f of the way there for orbital solutions!\n\n', counterOrbital, numLaskarRuns);
    counterOrbital = counterOrbital + 1;
    
    toc
end