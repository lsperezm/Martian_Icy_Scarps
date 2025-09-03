function [TotalFree, TotalForced, TotalSublimation] = sublCalculator(L04Step, s, Ls0, L04obl, atmosPressSF12, SLOPEDaverageDiurnalSurfTemps, REGminDiurnalSurfTemps, REGdiurnalTsurfCurves, dt, densityIce, avogadro, u_wind)
% FREE AND FORCED CONVECTION

%% Calculate Partial Pressure H2O scheme

    PPH2O_x=[1:60,119.5:.1:120.5,150:360];
    PPH2O_y=0*PPH2O_x+5;
    PPH2O_y(61:71)=75;
    PPH2O_x=PPH2O_x';
    PPH2O_y=PPH2O_y';
    PPH2O_f = fit(PPH2O_x,PPH2O_y,'smoothingspline');
    
    PPH2Otest_x = Ls0;
    PPH2Otest_y=PPH2O_f(PPH2Otest_x);
    PPH2Otest_y(find(PPH2Otest_y<1e-2)) = 0;
    mean(PPH2Otest_y);
    
    % g is the acceleration due to gravity (m/s^2)
    g = 3.71;
    HcH = 1;
    
    % Partial Pressure of water at present day
    PPH2O_present = g * PPH2Otest_y * 10^-6 * 1000 * (44/18) * (1/(1-exp(-HcH)));   % 1000 is the density of water
    
    obl_present = L04obl(end);
    
    SF12_a1 = -1.27; % Parameterization from Schorghofer and Forget (2012)
    SF12_b1 = 0.139; % Parameterization  from Schorghofer and Forget (2012)
    SF12_c1 = -0.00389; % Parameterization from Schorghofer and Forget (2012)
    oblDeg_present = rad2deg(obl_present);
    if oblDeg_present > 10 && oblDeg_present < 28
        atmosPressSF12_present = exp(SF12_a1 + SF12_b1*(oblDeg_present-28));
    elseif oblDeg_present > 28 && oblDeg_present < 50
        atmosPressSF12_present = exp(SF12_a1 + SF12_b1*(oblDeg_present-28) + SF12_c1*(oblDeg_present-28)^2);
    end
    atmosPressSF12_present = atmosPressSF12_present * s.atmosFactor;
    

% Vary relative to present day and Schorghofer and Forget 2012
if s.atmosFactor == 0
    PPH2O = PPH2O_present * 0; % 0 at all times
else
    PPH2O = PPH2O_present * (atmosPressSF12/atmosPressSF12_present);
end

%% Begin free and forced convection calculations
numDays_subl = size(SLOPEDaverageDiurnalSurfTemps,2);

% Switched to cells to be able to inspect diurnal curves if we want
%m_free_thickness = zeros(1,numDays);
%m_forced_thickness = zeros(1,numDays);
%total = zeros(1,numDays);
clear m_forced_thickness diurnal_m_forced_thickness m_free_thickness diurnal_m_free_thickness total summedForcedDayTotals summedFreeDayTotals summedDayTotals
summedDayTotals = zeros(1,numDays_subl);
summedForcedDayTotals = zeros(1,numDays_subl);
summedFreeDayTotals = zeros(1,numDays_subl);

counter = 1;
%% Loop over every day
for day=1:1:numDays_subl
    
    T_surf_subl = SLOPEDaverageDiurnalSurfTemps(day);
    
    reg_Tmin = REGminDiurnalSurfTemps(day);
    reg_all_Tsurfs = REGdiurnalTsurfCurves{day};
    
    numTimesteps_subl = size(reg_all_Tsurfs,1);
    
    diurnal_m_free_thickness = zeros(numTimesteps_subl,1);
    diurnal_m_forced_thickness = zeros(numTimesteps_subl,1);
    diurnal_total = zeros(numTimesteps_subl,1);
    
    for timestep = 1:1:numTimesteps_subl
        reg_Tsurf = reg_all_Tsurfs(timestep);
        
        %% Compute temperatures used in equations
        % T_surf is the ice surface temperature (K) according to Dundas 2010, BUT
        % in our case surface temp is not necessarily ice temp since our ice is
        % buried by a little bit of dust, so will just use the ice table temp
        
        % T_atm is the near-surface atmospheric temperature (K) and is composed of
        % Tmin, Treg and b
        % Tmin is the most recent diurnal minimum in the regional surface temp
        % Treg is the current regional surface temperature
        % b = 0.2 was found to roughly match the ~20-25K peak temp diff observed
        % by Viking sites in summertime
        T_min = reg_Tmin;
        T_reg = reg_Tsurf;
        b_exp = 0.2;
        T_atm_subl = T_min^b_exp * T_reg^(1-b_exp);
        
        % T_bl is the boundary layer temperature (K) and Dundas 2010 says according
        % to Hecht 2002 it can be set to the average of the ice (surface)
        % and atmospheric temperatures
        T_bl = (T_surf_subl + T_atm_subl)/2;
        
        %% Sensible heat forced
        % Equation 5 from Dundas 2010: m_forced = (Mw/k*T) * A * u_wind * (e_sat-e)
        % Mw is the molecular mass of water (kg)
        % k is the Boltzmann's constant (J/K)
        % T is the temperature of the boundary layer
        % A is the drag coefficient of ice (dimensionless)- is an equation but can
        % just use 0.002
        % u_wind is the wind speed (m/s)
        % e_sat is the saturated water vapor partial pressure (Pa)
        % I will use saturation vapor pressure using temperature Tsurf
        % e is the atmospheric water vapor partial pressure (Pa)
        % I will use the value outputed by Schorghofer and Forget 2012
        % but for now, Bapst was saying Colin said it was ~1 Pa...?
        molecMassWater = 2.99e-26;
        boltzmann = 1.38064852e-23; % Boltzmann Constant, m2 kg s-2 K-1
        T_subl = T_bl;
        A_dragcoef = 0.002;
        % u_wind = 2.5; % Dundas 2010 varied wind speed from 0-7.5 m/s This
        % value now input in param file
        e_sat = 611 .* exp( (-51058/8.31).*(1./T_surf_subl - 1/273.16) ); % equation to compute saturation vapor pressure for given temperature
        
        %e = 1; % Pa
        e_watvap = PPH2O(counter);
        %e = atmosVapPressSF12;
        
        m_forced = (molecMassWater/(boltzmann*T_subl)) * A_dragcoef * u_wind * (e_sat-e_watvap); % (kg/m^2??)
        
        if m_forced > 0
            diurnal_m_forced_thickness(timestep) = (m_forced * dt / densityIce); % units of meters
        else
            diurnal_m_forced_thickness(timestep) = 0;
        end
        
        %% Sensible heat free convection
        % Equation 6 from Dundas 2010: m_free = 0.14 * deltaNu * eta_ave * D *((delta_rho/rho)*(g/nu^2) *(nu/D))^(1/3)
        % P_atm is the atmospheric pressure (Pa) - converted 7 mbar to Pa but Mars'
        % atmospheric surface pressure ranges from ~4-8.7 mbar depending on
        % season and location according to https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
        % moved to above for calculated saturated surface gas density
        % rho_surf is the density of saturated surface gas (kg/m^3)
        % rho_atm is the atmospheric density (kg/m^3)
        % didn't use 0.02 kg/m^3 from https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
        % will calculate ourselves adding rho_H2O + rho_CO2 but basically get 0.02 anyways
        % rho_ave is the average of rho_atm and rho_surf (kg/m^3)
        P_atm = 700;
        
        molarMassH2O = 0.01801528; % kg/mole
        molarMassCO2 = 0.04401; %kg/mole
        % avogadro = 6.022140857e23; %molecules/mole already defined in TM12
        
        rho_surf = (((molarMassCO2/avogadro)*P_atm) + (((molarMassH2O-molarMassCO2)/avogadro)*e_sat)) / (boltzmann * T_surf_subl);
        rho_atm = (((molarMassCO2/avogadro)*P_atm) + (((molarMassH2O-molarMassCO2)/avogadro)*e_watvap)) / (boltzmann * T_atm_subl);
        rho_ave = (rho_atm + rho_surf)/2;
        
        % deltaEta is the difference between atmospheric and surface gas water mass
        % fraction (dimensionless)
        % I really don't understand what to do here so, if this is the same as
        % delta C from Hecht 2002 then can use 0.026 kg/m^3 divided by total
        % density of Mars atmsphere to make dimensionless??
        %deltaEta = 0.026/rho_atm;
        % or from Chittenden 2008: deltaEta = rho_sat - rho_atm / rho_ice
        %deltaEta = (rho_surf_H2O - rho_atm) / densityIce;
        % or from Chittenden 2008 but different values;
        rho_bl_H2O = e_sat * (molarMassH2O / avogadro) / (boltzmann * T_bl);
        rho_atm_H2O = e_watvap * (molarMassH2O / avogadro) / (boltzmann * T_atm_subl);
        rho_surf_H2O = e_sat * (molarMassH2O / avogadro) / (boltzmann * T_surf_subl);
        deltaEta = (rho_surf_H2O - rho_atm_H2O) / rho_atm;
        
        % D is the diffusion coefficient for H2O in CO2 (m^2/s)
        DiffCoef_H2OinCO2 = (1.387e-5)*((T_bl / 273.15)^(3/2)) * (10^5/P_atm);
        
        % delta rho/rho is the atmospheric and surface gas density difference,
        % divided by a reference density (dimensionless)
        % molarMassCO2 (m_c in the eq) is the molar mass of CO2 (kg)
        % molarMassH2O(m_w in the eq) is the molar mass of water (kg)
        % T_surf, T_atm, P_atm, e_sat and e are given above
        delta_rho_over_rho_num = (molarMassCO2*P_atm*((T_surf_subl/T_atm_subl)-1)) + ((molarMassCO2-molarMassH2O)*(e_sat - ((T_surf_subl/T_atm_subl)*e_watvap)));
        delta_rho_over_rho_denom = 0.5*((molarMassCO2*P_atm*((T_surf_subl/T_atm_subl)+1)) - ((molarMassCO2-molarMassH2O)*((T_surf_subl/T_atm_subl)*e_watvap + e_sat)));
        delta_rho_over_rho = delta_rho_over_rho_num/delta_rho_over_rho_denom;
        if delta_rho_over_rho < 0
            delta_rho_over_rho = 0;
        end
        
        % g is the acceleration due to gravity (m/s^2)
        g = 3.71;
        
        % nu is the kinematic viscosity (m^2/s)
        % R is the universal gas constant (J/(K mol))
        % m_c, T_bl and P_atm are given above
        R_univGasConst = 8.314;
        nu = (1.48e-5)*(R_univGasConst*T_bl/(molarMassCO2*P_atm))*((240+293.15)/(240+T_bl))*(T_bl/293.15)^(3/2);
        
        x = delta_rho_over_rho*(g/nu^2)*(nu/DiffCoef_H2OinCO2);
        xcuberoot = sign(x).*abs(x.^(1/3));
        m_free = 0.14 * deltaEta * rho_ave * DiffCoef_H2OinCO2 * xcuberoot; % (kg/m^2??)
        
        if m_free > 0
            diurnal_m_free_thickness(timestep) = (m_free * dt / densityIce); % units of meters
        else
            diurnal_m_free_thickness(timestep) = 0;
        end
        
        TSURFS(counter) = T_surf_subl;
        TATMS(counter) = T_atm_subl;
        TBLS(counter) = T_bl;
        ESATS(counter) = e_sat;
        RHOAVES(counter) = rho_ave;
        DELTAETAS(counter) = deltaEta;
        DS(counter) = DiffCoef_H2OinCO2;
        DELTARHOOVERRHOS(counter)= delta_rho_over_rho;
        NUS(counter) = nu;
        XCUBEROOTS(counter) = xcuberoot;
        MFORCED(counter) = m_forced;
        MFREE(counter) = m_free;
        counter = counter + 1;
    end
    
    m_forced_thickness{day} = diurnal_m_forced_thickness;
    m_free_thickness{day} = diurnal_m_free_thickness;
    total{day} = m_forced_thickness{day} + m_free_thickness{day};
    summedForcedDayTotals(day) = sum(diurnal_m_forced_thickness);
    summedFreeDayTotals(day) = sum(diurnal_m_free_thickness);
    summedDayTotals(day) = sum(total{day});
end

TotalFree = sum(summedFreeDayTotals);
TotalForced = sum(summedForcedDayTotals);
TotalSublimation = sum(summedDayTotals);

fprintf('\nFinal forced convection sublimation is %10.9f mm\n', TotalForced*1000);
fprintf('\nFinal free convection sublimation is %10.9f mm\n', TotalFree*1000);
fprintf('\nFinal total sublimation is %10.9f mm\n', TotalSublimation*1000);
end

