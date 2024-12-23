%% Rate of Split-belt Treadmill Adaptation Calculator
% Originating Author: Andrew Hagen
% Last Revised: 8/8/2024

% Uses multiple established methods to calculate rate of adaptation for step length asymmetry, double support asymmetry, and peak propulsion asymmetry. Each method has their advantages and disadvantages.
% For each variable, this script will calculate rate of adaptation using a single exponential model, a double exponential model, and a plateau threshold method
% Additionally, each method will be tested using no cropping of initial steps, and cropping the first 6 initial steps (to allow belts to speed up) - advantages and disadvantages to both in terms of model accuracy and rate calculations depending on your data
% Each method of rate calculations will also have a "raw data" output, and an output that utilized a rolling mean to pre-smooth the data.

% Toolboxes required: Optimization Toolbox, Global Optimization Toolbox, Particle Swarm Optimization Toolbox - Brian Birge
tic
clc 
clearvars -except SBNT_Data
close all

%% Data import and initialization (this example is from a matlab structure)

% load('SBNT_Data.mat')

% Initialize empty matrices to store the data
stepLengthAsym = [];
smoothedStepLengthAsym = [];
peakPropulsionAsym = [];
smoothedPeakPropulsionAsym = [];
dblSupAsym = [];
smoothedDblSupAsym = [];

% Initiate empty results tables with variable names
    % Define variable names for 6 data sets, 2 conditions (raw and cropped), and metrics (magnitude, rate, AIC or plateau and steps)
    dataSetNames = {'StepLengthAsym', 'SmoothedStepLengthAsym','PeakPropulsionAsym', 'SmoothedPeakPropulsionAsym','DblSupAsym', 'SmoothedDblSupAsym'};
    conditions = {'Raw', 'Crop'};
    metrics = {'Magnitude', 'Rate', 'AIC'};
    metricsPlateau = {'Plateau', 'Steps'};
    
    % Create full variable names
    doubleExpVarNames = {};
    singleExpVarNames = {};
    plateauVarNames = {};
    
    for i = 1:length(dataSetNames)
        for j = 1:length(conditions)
            for k = 1:length(metrics)
                % Construct variable name for exponential methods
                varName = strcat(dataSetNames{i}, '_', conditions{j}, '_', metrics{k});
                doubleExpVarNames{end+1} = varName;  % Append to list for double exponential model
                singleExpVarNames{end+1} = varName;  % Append to list for single exponential model
            end
        end
        for l = 1:length(metricsPlateau)
        % Construct variable name for plateau method
        varName = strcat(dataSetNames{i}, '_', metricsPlateau{l});
        plateauVarNames{end+1} = varName;  % Append to list for plateau method
        end
    end
    
% Initialize tables with explicit variable names including text elements
% Add columns for participant and visit names
doubleExpVarNames = [{'Participant', 'Visit', 'Phase'}, doubleExpVarNames];
singleExpVarNames = [{'Participant', 'Visit', 'Phase'}, singleExpVarNames];
plateauVarNames = [{'Participant', 'Visit', 'Phase'}, plateauVarNames];

% Initialize tables with the additional text columns
numRows = 196; % 49 participants each with 2 visits and 2 phases per visit
doubleExpTable = table('Size', [numRows, length(doubleExpVarNames)], 'VariableTypes', [repmat("string", 1, 3), repmat("double", 1, length(doubleExpVarNames) - 3)],'VariableNames', doubleExpVarNames);

singleExpTable = table('Size', [numRows, length(singleExpVarNames)], 'VariableTypes', [repmat("string", 1, 3), repmat("double", 1, length(singleExpVarNames) - 3)], 'VariableNames', singleExpVarNames);

plateauTable = table('Size', [numRows, length(plateauVarNames)], 'VariableTypes', [repmat("string", 1, 3), repmat("double", 1, length(plateauVarNames) - 3)], 'VariableNames', plateauVarNames);

%% Loop through each Participant, Visit, and Phase

% Loop through each Participant
participants = fieldnames(SBNT_Data);
for p = 1:length(participants)
participant = participants{p};
fittedTimeseries.(participant) = struct(); % Initialize layer for nested structure

    % Loop through each Visit for the current Participant
    visits = fieldnames(SBNT_Data.(participant));
    for v = 1:length(visits)
    visit = visits{v};
    fittedTimeseries.(participant).(visit) = struct(); % Initialize layer for nested structure

        % Loop through the Phases 'Adaptation' and 'Deadaptation'
        phases = {'Adaptation', 'Deadaptation'};
        for ph = 1:length(phases)
        phase = phases{ph};
        fittedTimeseries.(participant).(visit).(phase) = struct(); % Initialize layer for nested structure

             % Extract the data for StepLengthAsymNorm, PeakPropulsionAsym, and DblSupAsym and also create rolling mean presmoothed copies
            if isfield(SBNT_Data.(participant).(visit).(phase), 'StepLengthAsym_Norm')
                stepLengthAsym = SBNT_Data.(participant).(visit).(phase).StepLengthAsym_Norm;
                smoothedStepLengthAsym = movmean(stepLengthAsym, 5, 'Endpoints', 'shrink');
            end
            
            if isfield(SBNT_Data.(participant).(visit).(phase), 'PeakPropulsionAsym')
                peakPropulsionAsym = SBNT_Data.(participant).(visit).(phase).PeakPropulsionAsym;
                smoothedPeakPropulsionAsym = movmean(peakPropulsionAsym, 5, 'Endpoints', 'shrink');
            end
            
            if isfield(SBNT_Data.(participant).(visit).(phase), 'DoubleSupportAsym')
                dblSupAsym = SBNT_Data.(participant).(visit).(phase).DoubleSupportAsym;
                smoothedDblSupAsym = movmean(dblSupAsym, 5, 'Endpoints', 'shrink');
            end

             % Initialize row vectors to store the results for the current participant and visit
             resultsRowDoubleExp = nan(1, 36);  % 39 columns for 6 data sets * 2 conditions (raw and cropped) * 3 metrics ('Magnitude', 'Rate', 'AIC') 
             resultsRowSingleExp = nan(1, 36);  % 39 columns for 6 data sets * 2 conditions (raw and cropped) * 3 metrics ('Magnitude', 'Rate', 'AIC') 
             resultsRowPlateau = nan(1, 12);    % 27 columns for 6 data sets * 2 metrics ('Plateau', 'Steps') 

            % Data sets to process
            dataSets = {stepLengthAsym, smoothedStepLengthAsym, peakPropulsionAsym, smoothedPeakPropulsionAsym, dblSupAsym, smoothedDblSupAsym};
            dataSetNames = {'stepLengthAsym', 'smoothedStepLengthAsym', 'peakPropulsionAsym', 'smoothedPeakPropulsionAsym', 'dblSupAsym', 'smoothedDblSupAsym'};
            % Loop over each of the 6 data sets for different rate calculation methods
            for i = 1:length(dataSets)
                dataSet = dataSets{i};
                dataSetName = dataSetNames{i};
                
                %% Exponential modeling to calculate rate and magnitude of adaptation  (SLA Dblsup and PropAsym): Exponential fitting 
                % Rate of adaptation is also a measure of sensorimotor integration: The rate that SLA is restored represents the rate at which the nervous system integrates sensory feedback and updates motor output - Kuhman et al. 2022 IBRO Neuroscience Reports
                % This fits exponential models data timeseries using the particle swarm algorithm - Rashid et al. 2020. Brain Sciences.
                % fitSingleExpPSO and fitDoubleExpPSO functions are from Rashid et al. and only slightly adapted by myself - https://github.com/GallVp/knkTools/tree/master/exponentialModels 
               
                % Parameters of the fitted model(FitParams): [a_s, b_s, a_f, b_f, c] 
                % In this usage the parameters of the exponential models represent:
                % The initial asymmetry: (a_s + a_f + c)
                % Change in symmetry (magnitude) from beginning to the end of the training phase: (a_s + a_f)
                % The number of steps taken to complete 50% of the change: (ln(2)/abs(b_s))
                % The asymmetry at the end of the training phase: (c)
    
                % % Initialize structure nested layers to save fitted timseries data for exponential models
                fittedTimeseries.(participant).(visit).(phase).(dataSetName) = struct();
                fittedTimeseries.(participant).(visit).(phase).(dataSetName).DoubleExp = struct();
                fittedTimeseries.(participant).(visit).(phase).(dataSetName).SingleExp = struct();

                % Process raw data (1:end) for double-exponential fit
                if ~isempty(dataSet)
                    [Fitted, ~, FitParams, ~, Aic] = fitDoubleExpPSO(dataSet(1:end));
                    AdaptMagnitude = FitParams(1) + FitParams(3);
                    AdaptRate = log(2) / abs(FitParams(2));
                    % Store results in the appropriate columns for the current data set
                    resultsRowDoubleExp((i-1)*6 + 1) = AdaptMagnitude;
                    resultsRowDoubleExp((i-1)*6 + 2) = AdaptRate;
                    resultsRowDoubleExp((i-1)*6 + 3) = Aic;
                    fittedTimeseries.(participant).(visit).(phase).(dataSetName).DoubleExp.Raw = Fitted;
    
                    % Process cropped data (7:end) for double-exponential fit
                    [Fitted, ~, FitParams, ~, Aic] = fitDoubleExpPSO(dataSet(7:end));
                    AdaptMagnitude = FitParams(1) + FitParams(3);
                    AdaptRate = log(2) / abs(FitParams(2));
                    % Store results in the appropriate columns for the current data set
                    resultsRowDoubleExp((i-1)*6 + 4) = AdaptMagnitude;
                    resultsRowDoubleExp((i-1)*6 + 5) = AdaptRate;
                    resultsRowDoubleExp((i-1)*6 + 6) = Aic;
                    fittedTimeseries.(participant).(visit).(phase).(dataSetName).DoubleExp.Crop = Fitted;
                end
                
                % Process raw data (1:end) for single-exponential fit
                if ~isempty(dataSet)
                    [Fitted, ~, FitParams, ~, Aic] = fitSingleExpPSO(dataSet(1:end));
                    AdaptMagnitude = FitParams(1);
                    AdaptRate = log(2) / abs(FitParams(2));
                    % Store results in the appropriate columns for the current data set
                    resultsRowSingleExp((i-1)*6 + 1) = AdaptMagnitude;
                    resultsRowSingleExp((i-1)*6 + 2) = AdaptRate;
                    resultsRowSingleExp((i-1)*6 + 3) = Aic;
                    fittedTimeseries.(participant).(visit).(phase).(dataSetName).SingleExp.Raw = Fitted;
                    
                    % Process cropped data (7:end) for single-exponential fit
                    [Fitted, ~, FitParams, ~, Aic] = fitSingleExpPSO(dataSet(7:end));
                    AdaptMagnitude = FitParams(1);
                    AdaptRate = log(2) / abs(FitParams(2));
                    % Store results in the appropriate columns for the current data set
                    resultsRowSingleExp((i-1)*6 + 4) = AdaptMagnitude;
                    resultsRowSingleExp((i-1)*6 + 5) = AdaptRate;
                    resultsRowSingleExp((i-1)*6 + 6) = Aic;
                    fittedTimeseries.(participant).(visit).(phase).(dataSetName).SingleExp.Crop = Fitted;
                end
    
                %% Plateau Threshold Rate calculation
                % Defined the plateau range for each participant as the mean ± 2 SDs of the last 30 strides (60 steps) of adaptation
                % Analyzes the adaptation data to find the number of steps until 10 consecutive steps fell within the participants plateau range. This step number is our rate.

                if ~isempty(dataSet)
                     % Define plateau range using mean ± 2 SD of the last 30 strides (60 steps)
                     plateauRangeMean = mean(dataSet(end-59:end));
                     plateauRangeSD = std(dataSet(end-59:end));
                     plateauRange = [plateauRangeMean - 2*(plateauRangeSD), plateauRangeMean + 2*(plateauRangeSD)];
        
                    % Check if the next 10 consecutive steps fall within the plateau range for raw data.
                    for step = 1:length(dataSet)
                        if all(dataSet(step:min(step+9, length(dataSet))) >= plateauRange(1) & dataSet(step:min(step+9, length(dataSet))) <= plateauRange(2))
                            % If condition is met, store the current step number as the steps to plateau and then end the loop
                            stepsToPlateau = step; 
                            break;
                        end
                    end
        
                    % Store raw plateau data in the plateauRow
                    resultsRowPlateau((i-1)*2 + 1) = plateauRangeMean;
                    resultsRowPlateau((i-1)*2 + 2) = stepsToPlateau;

                end
            end

        %% Add Result Rows to table
        % Add participant, visit, and phase to the begining of the row
        resultsRowDoubleExp = [{participant}, {visit}, {phase}, num2cell(resultsRowDoubleExp)];
        resultsRowSingleExp = [{participant}, {visit}, {phase}, num2cell(resultsRowSingleExp)];
        resultsRowPlateau = [{participant}, {visit}, {phase}, num2cell(resultsRowPlateau)];

        % Add row to the tables
        rowIndex = ((p - 1) * 2 * 2) + ((v - 1) * 2) + ph;
        doubleExpTable(rowIndex,:) = resultsRowDoubleExp;
        singleExpTable(rowIndex,:) = resultsRowSingleExp;
        plateauTable(rowIndex,:) = resultsRowPlateau;
    
% close participant, visit, phase, loop
        end
    end
end

%% Data Export

exportData = inputdlg('Save and Export Data to Group xlsx? Y/N: ');
if strcmp(exportData,'Y')

    % Export Rate Calculation Tables
      writetable(doubleExpTable, 'DoubleExp_RateTable.xlsx')
      writetable(singleExpTable, 'SingleExp_RateTable.xlsx')
      writetable(plateauTable, 'Plateau_RateTable.xlsx')

    % Export Fitted Timeseries Data
    % Define the length for NaN padding and target length for interpolation
    paddingLength = 1300;
    interpolationLength = 1000;
    
    % Export Fitted Timeseries Data
    participantNames = fieldnames(fittedTimeseries); % Get the participant names
    
    for i = 1:length(dataSetNames)
        for j = 1:length(conditions)
            % Prepare data for export
            dataToExportDoubleExp = [];
            dataToExportSingleExp = [];
            dataToExportDoubleExpInterp = [];
            dataToExportSingleExpInterp = [];
            
            for p = 1:length(participantNames)
                participantName = participantNames{p};
                visitNames = fieldnames(fittedTimeseries.(participantName)); % Get the visit names for the current participant
                
                for v = 1:length(visitNames)
                    visitName = visitNames{v};
                    
                    % Padding with NaN to 1300 and combining phases for DoubleExp
                    dataAdaptationDouble = fittedTimeseries.(participantName).(visitName).Adaptation.(dataSetNames{i}).DoubleExp.(conditions{j});
                    dataDeadaptationDouble = fittedTimeseries.(participantName).(visitName).Deadaptation.(dataSetNames{i}).DoubleExp.(conditions{j});
                    dataAdaptationDoublePadded = [dataAdaptationDouble, nan(1, paddingLength - length(dataAdaptationDouble))];
                    dataDeadaptationDoublePadded = [dataDeadaptationDouble, nan(1, paddingLength - length(dataDeadaptationDouble))];
                    dataCombinedDoublePadded = [dataAdaptationDoublePadded, dataDeadaptationDoublePadded];
                    
                    % Padding with NaN to 1300  and combining phases for SingleExp
                    dataAdaptationSingle = fittedTimeseries.(participantName).(visitName).Adaptation.(dataSetNames{i}).SingleExp.(conditions{j});
                    dataDeadaptationSingle = fittedTimeseries.(participantName).(visitName).Deadaptation.(dataSetNames{i}).SingleExp.(conditions{j});
                    dataAdaptationSinglePadded = [dataAdaptationSingle, nan(1, paddingLength - length(dataAdaptationSingle))];
                    dataDeadaptationSinglePadded = [dataDeadaptationSingle, nan(1, paddingLength - length(dataDeadaptationSingle))];
                    dataCombinedSinglePadded = [dataAdaptationSinglePadded, dataDeadaptationSinglePadded];
                    
                    % Interpolation (or extrapolation) to 1000 points and 
                    dataAdaptationDoubleInterp = interp1(1:length(dataAdaptationDouble), dataAdaptationDouble, linspace(1, length(dataAdaptationDouble), interpolationLength), 'linear', 'extrap');
                    dataDeadaptationDoubleInterp = interp1(1:length(dataDeadaptationDouble), dataDeadaptationDouble, linspace(1, length(dataDeadaptationDouble), interpolationLength), 'linear', 'extrap');
                    dataCombinedDoubleInterp = [dataAdaptationDoubleInterp, dataDeadaptationDoubleInterp];
                    
                    dataAdaptationSingleInterp = interp1(1:length(dataAdaptationSingle), dataAdaptationSingle, linspace(1, length(dataAdaptationSingle), interpolationLength), 'linear', 'extrap');
                    dataDeadaptationSingleInterp = interp1(1:length(dataDeadaptationSingle), dataDeadaptationSingle, linspace(1, length(dataDeadaptationSingle), interpolationLength), 'linear', 'extrap');
                    dataCombinedSingleInterp = [dataAdaptationSingleInterp, dataDeadaptationSingleInterp];
    
                    % Append to export matrix
                    dataToExportDoubleExp = [dataToExportDoubleExp; [{participantName}, {visitName}, num2cell(dataCombinedDoublePadded)]];
                    dataToExportSingleExp = [dataToExportSingleExp; [{participantName}, {visitName}, num2cell(dataCombinedSinglePadded)]];
                    dataToExportDoubleExpInterp = [dataToExportDoubleExpInterp; [{participantName}, {visitName}, num2cell(dataCombinedDoubleInterp)]];
                    dataToExportSingleExpInterp = [dataToExportSingleExpInterp; [{participantName}, {visitName}, num2cell(dataCombinedSingleInterp)]];
                end
            end
            
            % Export to Excel with a sheet for each dataset and condition
            sheetNameDoublePadded = sprintf('%s_%s_DExp_Pad', dataSetNames{i}, conditions{j});
            sheetNameSinglePadded = sprintf('%s_%s_SExp_Pad', dataSetNames{i}, conditions{j});
            sheetNameDoubleInterp = sprintf('%s_%s_DExp_Intp', dataSetNames{i}, conditions{j});
            sheetNameSingleInterp = sprintf('%s_%s_SExp_Intp', dataSetNames{i}, conditions{j});
            
            writecell(dataToExportDoubleExp, 'DoubleExp_FittedData_Padded.xlsx', 'Sheet', sheetNameDoublePadded);
            writecell(dataToExportSingleExp, 'SingleExp_FittedData_Padded.xlsx', 'Sheet', sheetNameSinglePadded);
            writecell(dataToExportDoubleExpInterp, 'DoubleExp_FittedData_Interp.xlsx', 'Sheet', sheetNameDoubleInterp);
            writecell(dataToExportSingleExpInterp, 'SingleExp_FittedData_Interp.xlsx', 'Sheet', sheetNameSingleInterp);
        end
    end
    disp('Exported Successfully :)')
end

toc

%% Functions

    function [ fitted, resids, fitParams, linearCIs, aic] = fitDoubleExpPSO(dataSeries, fixJthTheta, atThetaJ)%% fitDoubleExpPSO Fit the double exponential model using the particle swarm optimization algorithm
    %% fitDoubleExpPSO Fit the double exponential model using the particle swarm optimization algorithm
    
        %   Inputs:
        %   a. dataSeries: Step length symmetry data as a vector.
        %   b. fixThetaJ (optional): Index of the jth theta. This is used when
        %   confidence intervals are to be obtained using the method which does not
        %   make the assumption of linearisation.
        %   c. atThetaJStar (optional): Vale to be used for the jth theta.
        %
        %   Output:
        %   a. fitted: Fitted values of symmetry.
        %   b. resids: Model residuals.
        %   c. fitParams: Parameters of the fitted model. [a_s, b_s, a_f, b_f, c]
    
        %   d. linearCIs: CIs based on the assumption of linearisation.
        %   e. aic: AIC value for the fitted model.
        %
        %   Toolboxes required: Optimization Toolbox, Global Optimization Toolbox, Particle Swarm Optimization Toolbox - Brian Birge
        %
        %   Copyright (c) <2018> <Usman Rashid>
        %   Licensed under the MIT License. See LICENSE in the project root for
        %   license information.
        %
        %   https://github.com/GallVp/knkTools/tree/master/exponentialModels
        %   Article on Model:  https://doi.org/10.3390/brainsci10100737
        %
        %   Adapted by Andrew Hagen
        
        %% Constants
        NUM_DCSN_STRDS       = 50; % Number of strides used to detect overall trend
        LAMBDA               = 1e3;
        
        
        %% Detect direction of symmetry series and set upper and lower bounds
        
        % These bounds correspond to No-overshoot.
        if(mean(dataSeries(1:NUM_DCSN_STRDS)) < mean(dataSeries(length(dataSeries) - NUM_DCSN_STRDS:end)))
            lb              = [-1,   -log(2),   -1,     -log(2),    -1];
            ub              = [0,    0,         0,      0,          1];
        else
            lb              = [0,   -log(2),    0,      -log(2),    -1];
            ub              = [1,    0,         1,      0,          1];
        end
        
        % Fix jth theta at the provided value if additional arguments are passed
        if nargin > 1
            lb(fixJthTheta) = atThetaJ;
            ub(fixJthTheta) = atThetaJ;
        end
        
        %% Set up particle swarm algorithm for optimisation
        rng shuffle;
        xData               = 1:length(dataSeries);
        fun                 = @(x)expCostFunc(xData, dataSeries, x);
        numParam            = 5;
        hybridopts          = optimoptions('fmincon', 'Display', 'none',...
            'Algorithm', 'interior-point',...
            'FunctionTolerance', 1e-6);
        options             = optimoptions('particleswarm',...
            'SwarmSize', numParam * 3,...
            'UseParallel', false,...
            'Display', 'none',...
            'HybridFcn', {@fmincon, hybridopts},...
            'FunctionTolerance', 1e-3);
        
        %% Run optimisation twice and chose the solution with the smaller cost
        fitParam1           = particleswarm(fun, numParam, lb, ub, options);
        fitParam2           = particleswarm(fun, numParam, lb, ub, options);
        c1                  = expCostFunc(xData, dataSeries, fitParam1);
        c2                  = expCostFunc(xData, dataSeries, fitParam2);
        
        if c1 < c2
            fitParamA       = fitParam1;
            cA              = c1;
        else
            fitParamA       = fitParam2;
            cA              = c2;
        end
        
        %% Detect direction of symmetry series and set upper and lower bounds
        
        % These bounds correspond to Yes-overshoot.
        if(mean(dataSeries(1:NUM_DCSN_STRDS)) < mean(dataSeries(length(dataSeries) - NUM_DCSN_STRDS:end)))
            lb              = [0,   -log(2),    -1,     -log(2),    -1];
            ub              = [1,    0,         0,      0,          1];
        else
            lb              = [-1,   -log(2),   0,      -log(2),    -1];
            ub              = [0,    0,         1,      0,          1];
        end
        
        % Fix jth theta at the provided value if additional arguments are passed
        if nargin > 1
            lb(fixJthTheta) = atThetaJ;
            ub(fixJthTheta) = atThetaJ;
        end
        
        %% Run optimisation twice and chose the solution with the smaller cost
        fitParam3           = particleswarm(fun, numParam, lb, ub, options);
        fitParam4           = particleswarm(fun, numParam, lb, ub, options);
        c3                  = expCostFunc(xData, dataSeries, fitParam3);
        c4                  = expCostFunc(xData, dataSeries, fitParam4);
        
        if c3 < c4
            fitParamB       = fitParam3;
            cB              = c3;
        else
            fitParamB       = fitParam4;
            cB              = c4;
        end
        
        % Chose the best solution from No-overshoot and Yes-overshoot solutions
        if cA < cB
            fitParams       = fitParamA;
        else
            fitParams       = fitParamB;
        end
        
        %% Assign outputs
        fitted              = fitParams(1) .* exp(xData.*fitParams(2)) + fitParams(3) .* exp(xData.*fitParams(4)) + fitParams(5);
        xData               = xData';
        J                   = [exp(fitParams(2).*xData) xData.*fitParams(1).*exp(fitParams(2).*xData)...
            exp(fitParams(4).*xData) xData.*fitParams(3).*exp(fitParams(4).*xData)...
            ones(length(xData), 1)];
        resids              = dataSeries - fitted;
        linearCIs           = nlparci(fitParams, resids,'jacobian', J);
        aic                 = 2*(numParam+1) + length(resids)*log(sum(resids.^2));
        
        % Cost function
            function c      = expCostFunc(X, Y, P)
                fitY        = P(1) .* exp(X.*P(2)) + P(3) .* exp(X.*P(4)) + P(5);
                ratePen     = max(0, P(4) - P(2) + 1e-3).^2;
                amountPen   = max(0, abs(P(1) + P(3)) - 2).^2;
                c           = sum((Y - fitY).^2) + LAMBDA*ratePen + LAMBDA*amountPen;
            end
    
        %% Plot Results
%         subplot(2, 2, 1)
%         plot(1:length(dataSeries), dataSeries, 'k.');
%         hold
%         plot(1:length(fitted), fitted, 'r-', 'LineWidth', 1.5);
%         title('Double exp. model')
%         box off;
%         ylabel('Symmetry')
%         ax = axis;
%         axis([1 length(dataSeries) ax(3) ax(4)]);
%         xlabel('Stride No.')
%         
%         subplot(2, 2, 2)
%         plot(1:length(dataSeries), resids, 'k.');
%         box off;
%         ax = axis;
%         axis([1 length(dataSeries) ax(3) ax(4)]);
%         xlabel('Stride No.')
%         
%         subplot(2, 2, 3)
%         histogram(resids, 10, 'FaceColor', [1 1 1], 'LineWidth', 1.5);
%         box off;
%         xlabel('Symmetry');
%         ylabel('Frequency');
%         
%         subplot(2, 2, 4)
%         qh = qqplot(resids);
%         title('');
%         ylabel('Quantiles of residuals');
%         xlabel('Standard normal quantiles');
%         qh(1).MarkerEdgeColor = [0 0 0];
    
    end
    
    function [ fitted, resids, fitParams, linearCIs, aic] = fitSingleExpPSO(dataSeries, fixJthTheta, atThetaJ)
    %%fitSingleExpPSO Fit the single exponential model using the particle swarm algorithm
    
        %   Inputs:
        %   a. dataSeries: Step length symmetry data as a vector.
        %   b. fixThetaJ (optional): Index of the jth theta. This is used when
        %   confidence intervals are to be obtained using the method which does not
        %   make the assumption of linearisation.
        %   c. atThetaJStar (optional): Vale to be used for the jth theta.
        %
        %   Output:
        %   a. fitted: Fitted values of symmetry.
        %   b. resids: Model residuals.
        %   c. fitParams: Parameters of the fitted model. [a_s, b_s, a_f, b_f, c]
        %   d. linearCIs: CIs based on the assumption of linearisation.
        %   e. aic: AIC value for the fitted model.
        %
        %   Toolboxes required: Optimization Toolbox, Global Optimization Toolbox, Particle Swarm Optimization Toolbox - Brian Birge
        %
        %   Copyright (c) <2018> <Usman Rashid>
        %   Licensed under the MIT License. See LICENSE in the project root for
        %   license information.
        %
        %   https://github.com/GallVp/knkTools/tree/master/exponentialModels
        %   Article on Model:  https://doi.org/10.3390/brainsci10100737
        %
        %   Adapted by Andrew Hagen
        
     %% Constants
    NUM_DCSN_STRDS       = 50; % Number of strides used to detect overall trend
    
    %% Detect direction of symmetry series and set upper and lower bounds
    if(mean(dataSeries(1:NUM_DCSN_STRDS)) < mean(dataSeries(length(dataSeries) - NUM_DCSN_STRDS:end)))
        lb              = [-2,  -log(2),    -1];
        ub              = [0,   0,          1];
    else
        lb              = [0,   -log(2),    -1];
        ub              = [2,   0,          1];
    end
    
    % Fix jth theta at the provided value if additional arguments are passed
    if nargin > 1
        lb(fixJthTheta) = atThetaJ;
        ub(fixJthTheta) = atThetaJ;
    end
    
    %% Set up particle swarm algorithm for optimisation
    rng shuffle;
    xData               = 1:length(dataSeries);
    fun                 = @(x)expCostFunc(xData, dataSeries, x);
    numParam            = 3;
    hybridopts          = optimoptions('fmincon', 'Display', 'none',...
        'Algorithm', 'interior-point',...
        'FunctionTolerance', 1e-6);
    options             = optimoptions('particleswarm',...
        'SwarmSize', numParam * 3,...
        'UseParallel', false,...
        'Display', 'none',...
        'HybridFcn', {@fmincon, hybridopts},...
        'FunctionTolerance', 1e-3);
    
    %% Run optimisation twice and chose the solution with the smaller cost
    fitParam1           = particleswarm(fun, numParam, lb, ub, options);
    fitParam2           = particleswarm(fun, numParam, lb, ub, options);
    
    c1                  = expCostFunc(xData, dataSeries, fitParam1);
    c2                  = expCostFunc(xData, dataSeries, fitParam2);
    
    if c1 < c2
        fitParams       = fitParam1;
    else
        fitParams       = fitParam2;
    end
    
    %% Assign outputs
    fitted              = fitParams(1) .* exp(xData.*fitParams(2)) + fitParams(3);
    xData               = xData';
    J                   = [exp(fitParams(2).*xData) xData.*fitParams(1).*exp(fitParams(2).*xData) ones(length(xData), 1)];
    resids              = dataSeries - fitted;
    linearCIs           = nlparci(fitParams, resids,'jacobian', J);
    aic                 = 2*(numParam+1) + length(resids)*log(sum(resids.^2));
    
    %% Define the cost function
        function c      = expCostFunc(X, Y, P)
            fitY        = P(1) .* exp(X.*P(2)) + P(3);
            c           = sum((Y - fitY).^2);
        end
    
        %% Plot Results
%         subplot(2, 2, 1)
%         plot(1:length(dataSeries), dataSeries, 'k.');
%         hold
%         plot(1:length(fitted), fitted, 'r-', 'LineWidth', 1.5);
%         title('Single exp. model')
%         box off;
%         ylabel('Symmetry')
%         ax = axis;
%         axis([1 length(dataSeries) ax(3) ax(4)]);
%         xlabel('Stride No.')
%         
%         subplot(2, 2, 2)
%         plot(1:length(dataSeries), resids, 'k.');
%         box off;
%         ax = axis;
%         axis([1 length(dataSeries) ax(3) ax(4)]);
%         xlabel('Stride No.')
%         
%         subplot(2, 2, 3)
%         histogram(resids, 10, 'FaceColor', [1 1 1], 'LineWidth', 1.5);
%         box off;
%         xlabel('Symmetry');
%         ylabel('Frequency');
%         
%         subplot(2, 2, 4)
%         qh = qqplot(resids);
%         title('');
%         ylabel('Quantiles of residuals');
%         xlabel('Standard normal quantiles');
%         qh(1).MarkerEdgeColor = [0 0 0];
    
    end
