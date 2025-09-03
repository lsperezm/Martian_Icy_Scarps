    %% Create diurnal information but move stuff around first so that the day for Ls = 0 starts first
    % but at the first timestep past midnight for which the day occurs
    % where Ls becomes 0
    
    Ls0_startOfDay_temp = find(floor(lsWrapped)==0,1); % First index after Ls wrapped from 360 back to 0
    currenthr_temp = hr(Ls0_startOfDay_temp);
    if Ls0_startOfDay_temp == 1
        beforehr_temp = hr(nStepsInYear);
    else
        beforehr_temp = hr(Ls0_startOfDay_temp-1);
    end
    while ~(currenthr_temp < 0 && beforehr_temp > 0)
        if Ls0_startOfDay_temp == 1
            Ls0_startOfDay_temp = nStepsInYear;
        else
            Ls0_startOfDay_temp = Ls0_startOfDay_temp - 1;
        end
        currenthr_temp = hr(Ls0_startOfDay_temp);
        if Ls0_startOfDay_temp == 1
            beforehr_temp = hr(nStepsInYear);
        else
            beforehr_temp = hr(Ls0_startOfDay_temp-1);
        end
    end
    Ls0_startOfDay = Ls0_startOfDay_temp;
    
    %% Wrap main variables around to start with the start of the day that strattles Ls = 0
    % Ones not used are commented out but could use those for certain
    % purposes
    clear TsurfLs0 Ls0 hrLs0
    TsurfLs0(1:nStepsInYear-Ls0_startOfDay+1) = Tsurf(Ls0_startOfDay:end);
    TsurfLs0(nStepsInYear-Ls0_startOfDay+2:nStepsInYear) = Tsurf(1:Ls0_startOfDay-1);
    Ls0(1:nStepsInYear-Ls0_startOfDay+1) = lsWrapped(Ls0_startOfDay:end);
    Ls0(nStepsInYear-Ls0_startOfDay+2:nStepsInYear) = lsWrapped(1:Ls0_startOfDay-1);
    %iceTempsLs0(1:nStepsInYear-Ls0_startOfDay+1) = iceTableTemps(Ls0_startOfDay:end);
    %iceTempsLs0(nStepsInYear-Ls0_startOfDay+2:nStepsInYear) = iceTableTemps(1:Ls0_startOfDay-1);
    %FrostMassLs0(1:nStepsInYear-Ls0_startOfDay+1) = frostMasses(Ls0_startOfDay:end);
    %FrostMassLs0(nStepsInYear-Ls0_startOfDay+2:nStepsInYear) = frostMasses(1:Ls0_startOfDay-1);
    %TempsLs0(:, 1:nStepsInYear-Ls0_startOfDay+1) = Temps(:, Ls0_startOfDay:end);
    %TempsLs0(:, nStepsInYear-Ls0_startOfDay+2:nStepsInYear) = Temps(:, 1:Ls0_startOfDay-1);
    hrLs0(1:nStepsInYear-Ls0_startOfDay+1) = hr(Ls0_startOfDay:end); % hr is from - pi (midnight) to pi (back to midnight)
    hrLs0(nStepsInYear-Ls0_startOfDay+2:nStepsInYear) = hr(1:Ls0_startOfDay-1);
    %ltstLs0(1:nStepsInYear-Ls0_startOfDay+1) = ltst(Ls0_startOfDay:end); % ltst is from 0 (midnight) to 24 (with 12 = noon)
    %ltstLs0(nStepsInYear-Ls0_startOfDay+2:nStepsInYear) = ltst(1:Ls0_startOfDay-1);
    % Don't use these because the temperature calculator doesn't use the
    % wrapped/reorganized values
    %flatIRSavedLs0(1:nStepsInYear-Ls0_startOfDay+1) = flatIRSaved(Ls0_startOfDay:end);
    %flatIRSavedLs0(nStepsInYear-Ls0_startOfDay+2:nStepsInYear) = flatIRSaved(1:Ls0_startOfDay-1);
    %flatVisSavedLs0(1:nStepsInYear-Ls0_startOfDay+1) = flatVisSaved(Ls0_startOfDay:end);
    %flatVisSavedLs0(nStepsInYear-Ls0_startOfDay+2:nStepsInYear) = flatVisSaved(1:Ls0_startOfDay-1);
    
    clear beginDayIndex;
    beginDayIndex(1) = 1;
    dayIndex = 2;
    for n = 2:size(hrLs0, 2)
        if hrLs0(n) < 0 && hrLs0(n-1) > 0
            beginDayIndex(dayIndex) = n;
            dayIndex = dayIndex + 1;
        end
    end
    
    numDays = max(size(beginDayIndex));
    %averageDiurnalTemps = zeros(nLayers, numDays);
    averageDiurnalSurfTemps = zeros(1, numDays);
    
    % Also now computes minimum and maximum diurnal temperatures
    %minDiurnalTemps = zeros(nLayers, numDays);
    minDiurnalSurfTemps = zeros(1, numDays);
    %maxDiurnalTemps = zeros(nLayers, numDays);
    %maxDiurnalSurfTemps = zeros(1, numDays);
    diurnalLs = zeros(1, numDays);
    
    % Put each day's surface diurnal temperature curve into columns to be
    % easy to grab along with each minimum value (so index of the max and
    % min and diurnal Ls value corresponds to the column for that same
    % day's temps)
    numStepsInDay = beginDayIndex(2)-beginDayIndex(1);
    %diurnalTsurfCurves = zeros(numStepsInDay, numDays);
    
    for n = 1:numDays
    
        if n == numDays
            %averageDiurnalTemps(:,n) = mean(TempsLs0(:, beginDayIndex(n):size(Temps,2)), 2);
            averageDiurnalSurfTemps(n) = mean(TsurfLs0(beginDayIndex(n):size(Temps,2)));
            %minDiurnalTemps(:,n) = min(TempsLs0(:, beginDayIndex(n):size(Temps,2)),[], 2);
            minDiurnalSurfTemps(n) = min(TsurfLs0(beginDayIndex(n):size(Temps,2)));
            %maxDiurnalTemps(:,n) = max(TempsLs0(:, beginDayIndex(n):size(Temps,2)),[], 2);
            %maxDiurnalSurfTemps(n) = max(TsurfLs0(beginDayIndex(n):size(Temps,2)));
            diurnalLs(n) = median(Ls0(beginDayIndex(n):size(Temps,2)));
            %diurnalTsurfCurves(:,n) = TsurfLs0(beginDayIndex(n):size(Tsurf,2))';
            diurnalTsurfCurves{n} = TsurfLs0(beginDayIndex(n):size(Tsurf,2))';
        else
            %averageDiurnalTemps(:,n) = mean(TempsLs0(:, beginDayIndex(n):beginDayIndex(n+1)-1), 2);
            averageDiurnalSurfTemps(n) = mean(TsurfLs0(beginDayIndex(n):beginDayIndex(n+1)-1));
            %minDiurnalTemps(:,n) = min(TempsLs0(:, beginDayIndex(n):beginDayIndex(n+1)-1), [], 2);
            minDiurnalSurfTemps(n) = min(TsurfLs0(beginDayIndex(n):beginDayIndex(n+1)-1));
            %maxDiurnalTemps(:,n) = max(TempsLs0(:, beginDayIndex(n):beginDayIndex(n+1)-1), [], 2);
            %maxDiurnalSurfTemps(n) = max(TsurfLs0(beginDayIndex(n):beginDayIndex(n+1)-1));
            diurnalLs(n) = median(Ls0(beginDayIndex(n):beginDayIndex(n+1)-1));
            %diurnalTsurfCurves(:,n) = TsurfLs0(beginDayIndex(n):beginDayIndex(n+1)-1)';
            diurnalTsurfCurves{n} = TsurfLs0(beginDayIndex(n):beginDayIndex(n+1)-1)';
        end
    
    end

    %averageDiurnalAllTemps = [averageDiurnalSurfTemps; averageDiurnalTemps];
    
   
% Things needed for free/forced convection calculation
if slope == 0
    REGminDiurnalSurfTemps = minDiurnalSurfTemps;
    REGdiurnalTsurfCurves = diurnalTsurfCurves;
    REGdiurnalLs = diurnalLs;
    
    %REGaverageDiurnalSurfTemps = averageDiurnalSurfTemps;
else
    SLOPEDaverageDiurnalSurfTemps = averageDiurnalSurfTemps;
    SLOPEDdiurnalTsurfCurves = diurnalTsurfCurves;
    SLOPEDdiurnalLs = diurnalLs;
    
    %SLOPEDminDiurnalSurfTemps = minDiurnalSurfTemps;
end