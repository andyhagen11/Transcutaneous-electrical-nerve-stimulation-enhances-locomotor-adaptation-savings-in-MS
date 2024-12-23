%% Vicon Nexus Split-Belt Data Analysis %%
% Originating Author: Andrew Hagen
% Last Revised: 8/20/2024

% This script uses data from Vicon Nexus during a split-belt treadmill adaptation paradigm to calculate a variety of gait parameters and their asymmetries including Phase Coordination Index (PCI)...
% double support percent, step length asymmetry (SLA), ground reaction forces, joint angles, joint moments, and joint powers

% Calculates rate and magnitude of adaptation or deadaptation using exponential fitting for step length asymmetry, double support percent, and propulsion
    % Further outputs xlsx files with the timeseries of adaptation for these variables

% This script was used in a study in people with multiple sclerosis, and labeled each limb as more affected or less unaffected

% Prompts: Condition - Visit1_NoTENS, Visit1_TENS, Visit1_NoTENS, Visit2_TENS
%          Phase - Baseline, Adaptation, Deadaptation
%          Affected Limb: L or R
%          Body mass (kg): To calculate forces relative to body weight

% At the end of the script, creates multi-layer nested data structure that is additive for all participants 
    % Additionally pumps out another data structure with each gait cycle interpolated to be 100 data points long

    %Layer 1 = Participant
    %Layer 2 = Visit (Visit1_NoTENS... etc.)
    %Layer 3 = Phase (Baseline (2min), Adaptation (10min) , Deadaptation (10 min) - will have to run seperately for each phase since they are separate Vicon trials
    %Layer 4 = More affected or Less affected
    %Layer 5 = Joint/force/variable
    %Layer 6 = Stance or Swing
    %Layer 7 = X, Y or Z plane (only for force and joints)
    %Layer 8 = Data with each row being one gait cycle

    %Additionally calculated adaptation series data and TimepointData is added to the strucutre but does not have multiple layers

%% Script initalization

% IF YOU WANT TO ADD TO PREVIOUS STRUCTURE MAKE SURE YOU HAVE THE DATA STRUCTURES SBNT_Data and SBNT_Data_Interp LOADED IN THE WORKSPACE
    % It takes way too long to reload it for each run - so keep it loaded and make sure it doesn't clear 

% Must have a Nexus trail open on this machine
    % This trail must have gait events (heel strike and toe off) labeled and the dynamic gait model processed
    % Using the MATLAB integration in the Nexus pipeline makes this easy
   
% Table variable abbreviations
    % LA = Less Affected Limb
    % MA =  More Affected limb
    % B = Baseline - last 30 strides 
    % EAd = Early Adapt- strides 6–30 
    % MAd = Mid Adapt - middle 30 strides
    % LAd = Late Adapt - last 30 strides
    % EDe = Early Deadapt - strides 11–31
    % MDe = Mid Deadapt - middle 30 strides
    % LDe = Late Deadapt - last 30 strides

clc 
clearvars -except SBNT_Data SBNT_Data_Interp 
close all

 % User Initializes Visit and Phases
     VisitChoices = {'Visit1_NoTENS', 'Visit1_TENS', 'Visit2_NoTENS', 'Visit2_TENS'};
     selection = listdlg('PromptString', 'Select a visit:', 'SelectionMode', 'single', 'ListString', VisitChoices);
     Visit = VisitChoices{selection};
    Phase = questdlg('Select a phase:', 'User Input', 'Baseline', 'Adaptation', 'Deadaptation','Baseline');
    if strcmp(Phase,'Baseline') 
        Timepoints = {'B'};
    elseif strcmp(Phase,'Adaptation') 
        Timepoints = {'EAd','MAd','LAd'}; % Early mid and late adapt
    elseif strcmp(Phase,'Deadaptation') 
        Timepoints = {'EDe','MDe','LDe'}; % Early mid and late deadapt
    end

%% Set up connection with Vicon Nexus and set path
 %Make sure Nexus is open with the desired trial 
    vicon = ViconNexus;

    SubjectName = vicon.GetSubjectNames;
    SubNum = char(SubjectName);
    Rate = vicon.GetFrameRate;
    Time = 1/Rate;
    Alimb = inputdlg('Affected Limb? (L or R): '); 
      Mass = inputdlg('Body Mass in Kg: ');
      Mass = str2double(Mass{1});
      BodyWeight = (Mass * 9.81);

     RootPath = 'R:\Split-belt fNIRS TENS (SBNT) Study\Group Data\';
     cd(RootPath)

%% Import Heel Strike and Toe Off times - these come in frames (named according to affected limb)

if strcmp(Alimb,'L')
    MA_HeelStrikes = vicon.GetEvents(SubNum, 'Left', 'Foot Strike');
    LA_HeelStrikes = vicon.GetEvents(SubNum, 'Right', 'Foot Strike');
    MA_ToeOffs = vicon.GetEvents(SubNum, 'Left', 'Foot Off');
    LA_ToeOffs = vicon.GetEvents(SubNum, 'Right', 'Foot Off');

elseif strcmp(Alimb,'R')
    LA_HeelStrikes = vicon.GetEvents(SubNum, 'Left', 'Foot Strike');
    MA_HeelStrikes = vicon.GetEvents(SubNum, 'Right', 'Foot Strike');
    LA_ToeOffs = vicon.GetEvents(SubNum, 'Left', 'Foot Off');
    MA_ToeOffs = vicon.GetEvents(SubNum, 'Right', 'Foot Off');

end

%% Import Model Data (named according to affected limb)

if strcmp(Alimb,'L')
    MA_AnkleAng = vicon.GetModelOutput(SubNum, 'LAnkleAngles'); 
    LA_AnkleAng = vicon.GetModelOutput(SubNum, 'RAnkleAngles');
    MA_KneeAng = vicon.GetModelOutput(SubNum, 'LKneeAngles');
    LA_KneeAng = vicon.GetModelOutput(SubNum, 'RKneeAngles');
    MA_HipAng = vicon.GetModelOutput(SubNum, 'LHipAngles');
    LA_HipAng = vicon.GetModelOutput(SubNum, 'RHipAngles');
    MA_AnkleMom = vicon.GetModelOutput(SubNum, 'LAnkleMoment'); 
    LA_AnkleMom = vicon.GetModelOutput(SubNum, 'RAnkleMoment');
    MA_KneeMom = vicon.GetModelOutput(SubNum, 'LKneeMoment');
    LA_KneeMom = vicon.GetModelOutput(SubNum, 'RKneeMoment');
    MA_HipMom = vicon.GetModelOutput(SubNum, 'LHipMoment');
    LA_HipMom = vicon.GetModelOutput(SubNum, 'RHipMoment');
    MA_AnklePow = vicon.GetModelOutput(SubNum, 'LAnklePower'); 
    LA_AnklePow = vicon.GetModelOutput(SubNum, 'RAnklePower');
    MA_KneePow = vicon.GetModelOutput(SubNum, 'LKneePower');
    LA_KneePow = vicon.GetModelOutput(SubNum, 'RKneePower');
    MA_HipPow = vicon.GetModelOutput(SubNum, 'LHipPower');
    LA_HipPow = vicon.GetModelOutput(SubNum, 'RHipPower');
    LA_Fx = vicon.GetDeviceChannel(1,1,1);
    LA_Fy = vicon.GetDeviceChannel(1,1,2);
    LA_Fz = vicon.GetDeviceChannel(1,1,3);
    MA_Fx = vicon.GetDeviceChannel(2,1,1);
    MA_Fy = vicon.GetDeviceChannel(2,1,2);
    MA_Fz = vicon.GetDeviceChannel(2,1,3);

elseif strcmp(Alimb,'R')
    LA_AnkleAng = vicon.GetModelOutput(SubNum, 'LAnkleAngles'); 
    MA_AnkleAng = vicon.GetModelOutput(SubNum, 'RAnkleAngles');
    LA_KneeAng = vicon.GetModelOutput(SubNum, 'LKneeAngles');
    MA_KneeAng = vicon.GetModelOutput(SubNum, 'RKneeAngles');
    LA_HipAng = vicon.GetModelOutput(SubNum, 'LHipAngles');
    MA_HipAng = vicon.GetModelOutput(SubNum, 'RHipAngles');
    LA_AnkleMom = vicon.GetModelOutput(SubNum, 'LAnkleMoment'); 
    MA_AnkleMom = vicon.GetModelOutput(SubNum, 'RAnkleMoment');
    LA_KneeMom = vicon.GetModelOutput(SubNum, 'LKneeMoment');
    MA_KneeMom = vicon.GetModelOutput(SubNum, 'RKneeMoment');
    LA_HipMom = vicon.GetModelOutput(SubNum, 'LHipMoment');
    MA_HipMom = vicon.GetModelOutput(SubNum, 'RHipMoment');
    LA_AnklePow = vicon.GetModelOutput(SubNum, 'LAnklePower'); 
    MA_AnklePow = vicon.GetModelOutput(SubNum, 'RAnklePower');
    LA_KneePow = vicon.GetModelOutput(SubNum, 'LKneePower');
    MA_KneePow = vicon.GetModelOutput(SubNum, 'RKneePower');
    LA_HipPow = vicon.GetModelOutput(SubNum, 'LHipPower');
    MA_HipPow = vicon.GetModelOutput(SubNum, 'RHipPower');
    MA_Fx = vicon.GetDeviceChannel(1,1,1);
    MA_Fy = vicon.GetDeviceChannel(1,1,2);
    MA_Fz = vicon.GetDeviceChannel(1,1,3);
    LA_Fx = vicon.GetDeviceChannel(2,1,1);
    LA_Fy = vicon.GetDeviceChannel(2,1,2);
    LA_Fz = vicon.GetDeviceChannel(2,1,3);

end

    % Make force output in all planes into one matrix (like for joints)
    % Flipped X and Y to match our specific coordinate system of joints (e.g. Row 1 is flexion/extension for joint and propulsion/braking for force)
    MA_Force = [MA_Fy; MA_Fx; MA_Fz];
    LA_Force = [LA_Fy; LA_Fx; LA_Fz];

  disp('Data imported successfully.');

%% Gait Cycle Duration [s] : - Heel strike to ipsilateral Heel strike 

    MA_StrideT = (diff([MA_HeelStrikes(:,end) MA_HeelStrikes],[],2));
    LA_StrideT = (diff([LA_HeelStrikes(:,end) LA_HeelStrikes],[],2));
    
    MA_StrideT(:,1) = []; % Crop first columns (always 0)
    LA_StrideT(:,1) = []; 
    
    % Confirm stride count is the same so we can work with matrix
        if length(MA_StrideT) ~= length(LA_StrideT)
            if length(MA_StrideT) > length(LA_StrideT)
                MA_StrideT(length(MA_StrideT)) = [];
            else LA_StrideT(length(LA_StrideT)) = [];  
            end
        end
    
    % Convert from frames to seconds 
    MA_StrideTFrames = cast(MA_StrideT, "double");
    MA_StrideT = MA_StrideTFrames/Rate;
    LA_StrideTFrames = cast(LA_StrideT, "double");
    LA_StrideT = LA_StrideTFrames/Rate;
    

%% Step Duration [s] - Heel strike to contralateral heel strike

    %START W/ MA STEP
    if LA_HeelStrikes(1) > MA_HeelStrikes(1) % Starting w/ a MA step (i.e. MA heel strike)
        MA_StepDur= zeros(1,length(MA_HeelStrikes)-1);
        LA_StepDur= zeros(1,length(MA_HeelStrikes)-1);
                for i = 1:(length(MA_HeelStrikes)-1) % The # of heelstrikes will always be equal or +1 for MA_Heelstrikes depending on last step
                    MA_StepDur(i) = LA_HeelStrikes(i) - MA_HeelStrikes(i);
                    LA_StepDur(i) = MA_HeelStrikes(i+1) - LA_HeelStrikes(i);
                end
    
    %START w/ LA STEP
    elseif MA_HeelStrikes(1) > LA_HeelStrikes(1) % Starting w/ a LA step (i.e. LA heel strike)
            MA_StepDur= zeros(1,length(LA_HeelStrikes)-1);
            LA_StepDur= zeros(1,length(LA_HeelStrikes)-1);
                for i = 1:(length(LA_HeelStrikes)-1) % The # of heelstrikes will always be equal or +1 for LA_Heelstrikes depending on last step
                    LA_StepDur(i) = MA_HeelStrikes(i) - LA_HeelStrikes(i);
                    MA_StepDur(i) = LA_HeelStrikes(i+1) - MA_HeelStrikes(i);
                end
    end
    
    % Confirm step count is the same so we can work with matrix
        if length(MA_StepDur) ~= length(LA_StepDur)
            if length(MA_StepDur) > length(LA_StepDur)
                MA_StepDur(length(MA_StepDur)) = [];
            else LA_StepDur(length(LA_StepDur)) = [];    
            end
        end
    
    % Convert from frames to seconds 
    MA_StepDurFrames = cast(MA_StepDur, "double");
    MA_StepDur = MA_StepDurFrames/Rate;
    LA_StepDurFrames = cast(LA_StepDur, "double");
    LA_StepDur = LA_StepDurFrames/Rate;
 
%% Swing Time [s] - Toe off to ipsilateral heel strike

    if length(MA_HeelStrikes) == length(MA_ToeOffs) % Always will start with a toe off so heel strikes <= toe offs. Therefore do n or n-1 of MA_HeelStrikes iterations.
        MA_SwingT= zeros(1,length(MA_HeelStrikes)-1);
        for i = 1:(length(MA_HeelStrikes)-1) 
        MA_SwingT(i) = MA_HeelStrikes(i) - MA_ToeOffs(i); % Since we always start with a toe off we can subtract from same index.
        end
    else
        MA_SwingT= zeros(1,length(MA_HeelStrikes));
        for i = 1:(length(MA_HeelStrikes))
        MA_SwingT(i) = MA_HeelStrikes(i) - MA_ToeOffs(i); % Since we always start with a toe off we can subtract from same index.
        end
    end
    
    if length(LA_HeelStrikes) == length(LA_ToeOffs)
        LA_SwingT= zeros(1,length(LA_HeelStrikes)-1);% Always will start with a toe off so heel strikes > or = toe offs. Therefore do n or n-1 of LA_HeelStrikes iterations.
        for i = 1:(length(LA_HeelStrikes)-1) 
        LA_SwingT(i) = LA_HeelStrikes(i) - LA_ToeOffs(i); % Since we always start with a toe off we can subtract from same index.
        end
    else
        LA_SwingT= zeros(1,length(LA_HeelStrikes));
        for i = 1:length(LA_HeelStrikes) 
        LA_SwingT(i) = LA_HeelStrikes(i) - LA_ToeOffs(i); % Since we always start with a toe off we can subtract from same index.
        end
    end
    
    % Confirm swing count is the same so we can work with matrix
        if length(MA_SwingT) ~= length(LA_SwingT)
            if length(MA_SwingT) > length(LA_SwingT)
                MA_SwingT(length(MA_SwingT)) = [];
    
            else LA_SwingT(length(LA_SwingT)) = []; 
            end
        end
    
    % Convert from frames to seconds
    MA_SwingTFrames = cast(MA_SwingT, "double");
    MA_SwingT = MA_SwingTFrames/Rate;
    LA_SwingTFrames = cast(LA_SwingT, "double");
    LA_SwingT = LA_SwingTFrames/Rate;
    
    % Determining Short and Long Swing Times 
    if mean(MA_SwingT) < mean(LA_SwingT)
            ShortTime = MA_StepDur; LongTime = LA_StrideT;
        else
            ShortTime = LA_StepDur; LongTime = MA_StrideT;
    end

%% Calculating phi (degrees) and Phase Coordination Index (PCI) for each timepoint
 
     Phi= zeros(1,(length(LongTime)-1));
    for i = 1:(length(LongTime)-1)
        Phi(2*i-1) = 360*ShortTime(i)/LongTime(i);  
        Phi(2*i) = 360*ShortTime(i)/LongTime(i+1);
    end 

    % Loop through PCI calc for each timepoint
    TimepointData.Phi = TimepointDataGrab(Phi,Phase,Timepoints,'Step');
    for i = 1:length(Timepoints)
        CurrentTimepoint = Timepoints{i};

        Phi_diff = abs(TimepointData.Phi.(CurrentTimepoint) - 180); % Absolute difference of phi values
        Phi_abs = mean(Phi_diff); % Mean value of absolute differences (degrees)
        Pphi_abs = 100*Phi_abs/180; % Percentage-converted phi_ABS (%)
        Phi_cv = (std(Phi)/mean(Phi))*100; % Coefficient of variation of mean of phi (%)
        TimepointData.PCI.(CurrentTimepoint) = Phi_cv + Pphi_abs; % Calculate PCI
    end


 %% Calculating Terminal Double Support Time - (Contralateral HS to Ipsilateral TO) as a percentage of stride time
 % There are two periods of double support per stride cycle and we define slow double support as occurring at the end of the slow limb's stance (i.e., the time from fast leg HS to slow leg TO), and fast double support at the end of the fast limb's stance (i.e., the time from slow leg HS to fast leg TO).
   % Initialize variables
   MA_DblSup= zeros(1,length(MA_HeelStrikes)-1);
   LA_DblSup= zeros(1,length(LA_HeelStrikes)-1);
   DblSupAsym = zeros(1,(2*(length(MA_DblSup))-2));

    % START W/ MA STEP
    if LA_HeelStrikes(1) > MA_HeelStrikes(1) % Starting w/ a MA step (i.e. MA heel strike)
        for i = 1:(length(LA_HeelStrikes)-1) 
            MA_DblSup(i) = double((MA_ToeOffs(i+1)-LA_HeelStrikes(i)))/MA_StrideTFrames(i); % Since LA is second step we need to subtract from MA_ToeOffs i+1 
            LA_DblSup(i) = double((LA_ToeOffs(i)-MA_HeelStrikes(i)))/LA_StrideTFrames(i); % Since we always start with a toe off we can subtract from same index.
        end
            % Calculate Asymmetry (MA - LA)
        for i=1:length(LA_DblSup)-1
             DblSupAsym(2*i-1)= (MA_DblSup(i)-LA_DblSup(i));
             DblSupAsym(2*i)=(MA_DblSup(i)-LA_DblSup(i+1));
        end
    
    % START w/ LA STEP
    elseif MA_HeelStrikes(1) > LA_HeelStrikes(1) % Starting w/ a right step (i.e. right heel strike)
        for i = 1:(length(MA_HeelStrikes)-1) 
            MA_DblSup(i) = double((MA_ToeOffs(i)-LA_HeelStrikes(i)))/MA_StrideTFrames(i); % Since we always start with a toe off we can subtract from same index.
            LA_DblSup(i) = double((LA_ToeOffs(i+1)-MA_HeelStrikes(i)))/LA_StrideTFrames(i); %  Since MA is second step we need to subtract from LA_ToeOffs i+1 
        end
            % Calculate Asymmetry (MA - LA)
        for i=1:length(MA_DblSup)-1
             DblSupAsym(2*i-1)= (MA_DblSup(i)-LA_DblSup(i));
             DblSupAsym(2*i)=(MA_DblSup(i)-LA_DblSup(i+1));
        end
    end
    

    % Loop through each timepoint
    TimepointData.DoubleSupportAsym = TimepointDataGrab(DblSupAsym,Phase,Timepoints,'Step');

%% Step Length Asymetry (total, positional contribution, timing contribution, and velocity contribution)
    %Contribution analyis based on temporal and spatial decomposition of SLA detailed in:
    %(Finley et al. 2015) Spatial and Temporal Control Contribute to Step Length Asymmetry During Split-Belt Adaptation and Hemiparetic Gait

    % If you only set one variable = vicon.GetTrajectory it will output only the x direction. If you want all variables use "[X, Y, Z, Exists]" as the output
    % Because of our calibration, our x axis goes the width of the treadmill and y axis goes the length of the treadmill (opposite of what is typical) so we grab the y values from our data and label them as x values for this code. 

    [L_AnkY,L_Ankx] = vicon.GetTrajectory(SubNum, 'LANK'); % Switched the x and y axis because our calibration was opposite
    [R_AnkY,R_Ankx] = vicon.GetTrajectory(SubNum, 'RANK');
    [L_PsisY,L_Psisx] = vicon.GetTrajectory(SubNum, 'LPSI');
    [R_PsisY,R_Psisx] = vicon.GetTrajectory(SubNum, 'RPSI');

    %Rather than using ASIS markers (like in Finley et al) we will use an average of PSIS markers (proxy for COM) to avoid hip motion being a confound
    COM = mean([R_Psisx; L_Psisx], 1); 

    %Created a copy of MA & LA heel strike matrices as lhs and rhs so it is not affected for other analyses
    %Convert back to L and R instead of MA and LA to complete SLA calculations
    if strcmp(Alimb,'L')
        Rhs = LA_HeelStrikes;
        Lhs = MA_HeelStrikes;
        Rto = LA_ToeOffs;
        Lto = MA_ToeOffs;
    elseif strcmp(Alimb,'R')
        Lhs = LA_HeelStrikes;
        Rhs = MA_HeelStrikes;
        Lto = LA_ToeOffs;
        Rto = MA_ToeOffs;
    end

    %Determine Start Leg
    if Rhs(1)<Lhs(1)
        StartLeg = 1;
    else
        StartLeg = 2;
    end
    
    % Find trajectory locations at each heel strike relative to COM (ankle - COM)
    %There should always be one more slow leg heel strike than fast leg heel strike and this is made sure of by below if statement.
    %Check if the start leg matches the slow leg; if not then skip first event. Then check if slow leg has one more event than the fast leg: if not then remove final fast leg event because it will not be used in calculations.

    if StartLeg==1 && strcmpi(Alimb,'L')==1 %Right is first and slow
        if length(Rhs)<length(Lhs)
            Lhs = Lhs(1:length(Rhs));
            Lto = Lto(1:length(Rto));
        end
        if length(Rhs)==length(Lhs)
            Lhs = Lhs(1:end-1);
            Lto = Lto(1:end-1);
        end
        X_f = L_Ankx-COM; 
        X_s = R_Ankx-COM;
    
    elseif StartLeg==1 && strcmpi(Alimb,'R')==1 %Right is first but fast
        Rhs = Rhs(2:end);
        if length(Lhs)<length(Rhs)
            Rhs = Rhs(1:length(Lhs));
            Rto = Rto(1:length(Lto));
        end
        if length(Rhs)==length(Lhs)
            Rhs = Rhs(1:end-1);
            Rto = Rto(1:end-1);
        end
        X_s = L_Ankx-COM;
        X_f= R_Ankx-COM;
    
    elseif StartLeg==2 && strcmpi(Alimb,'R')==1 %Left is first and slow
        if length(Lhs)<length(Rhs)
            Rhs = Rhs(1:length(Lhs));
            Rto = Rto(1:length(Lto));
        end
        if length(Rhs)==length(Lhs)
            Rhs = Rhs(1:end-1);
            Rto = Rto(1:end-1);
        end
        X_s = L_Ankx-COM;
        X_f = R_Ankx-COM;
    
    elseif StartLeg==2 && strcmpi(Alimb,'L')==1 %Left is first but fast
        Lhs = Lhs(2:end);
        if length(Rhs)<length(Lhs)
            Lhs = Lhs(1:length(Rhs));
            Lto = Lto(1:length(Rto));
        end
        if length(Rhs)==length(Lhs)
            Lhs = Lhs(1:end-1);
            Lto = Lto(1:end-1);
        end
        X_f = L_Ankx-COM;
        X_s = R_Ankx-COM;
    end

    % Combine MA and LA variables 
    NumStrikes = [length(Rhs), length(Lhs)];
    CombinedHS = {Rhs,Lhs};
    CombinedTO = {Rto,Lto};
    CombinedAnk = {R_Ankx,L_Ankx};
    FastLeg = strcmpi(Alimb,{'R','L'});
    
    % Initialize variables to be used in loops
    Xs_SHS1 = zeros(1,NumStrikes(~FastLeg)-1); % Slow leg ankle-COM at slow leg heelstrike 1
    Xf_SHS1 = Xs_SHS1; % Fast leg ankle-COM at slow leg heelstrike 1
    Xs_SHS2 = Xs_SHS1; % Slow leg ankle-COM at slow leg heelstrike 2 (SHS1 shifted by one position)
    Xf_SHS2 = Xs_SHS1; % Fast leg ankle-COM at slow leg heelstrike 2 (SHS1 shifted by one position)
    Xf_FTO1 = Xs_SHS1; % Fast leg ankle-COM at fast leg toe off 1
    Xs_STO1 = Xs_SHS1; % Slow leg ankle-COM at slow leg toe off 1
    Xs_FHS = zeros(1,NumStrikes(FastLeg));
    Xf_FHS = Xs_FHS;
    T_s = zeros(1,NumStrikes(~FastLeg)-1);
    T_f = T_s;
    
    %Ankle location relative to COM location at slow leg heel strikes 
    for i = 1:NumStrikes(~FastLeg)-1
        Xs_SHS1(i) = CombinedAnk{~FastLeg}(CombinedHS{~FastLeg}(i))-COM(CombinedHS{~FastLeg}(i));
        Xf_SHS1(i) = CombinedAnk{FastLeg}(CombinedHS{~FastLeg}(i))-COM(CombinedHS{~FastLeg}(i));
        Xs_SHS2(i) = CombinedAnk{~FastLeg}(CombinedHS{~FastLeg}(i+1))-COM(CombinedHS{~FastLeg}(i+1));
        Xf_SHS2(i) = CombinedAnk{FastLeg}(CombinedHS{~FastLeg}(i+1))-COM(CombinedHS{~FastLeg}(i+1));
    end
    
    %Ankle location relative to COM location at slow leg toe off
    for i = 1:NumStrikes(~FastLeg)-1
        Xs_STO1(i) = CombinedAnk{~FastLeg}(CombinedTO{~FastLeg}(i))-COM(CombinedTO{~FastLeg}(i));
    end
    
    %Ankle location relative to COM location at fast leg heel strike
    for i = 1:NumStrikes(FastLeg)
        Xs_FHS(i) = CombinedAnk{~FastLeg}(CombinedHS{FastLeg}(i))-COM(CombinedHS{FastLeg}(i));
        Xf_FHS(i) = CombinedAnk{FastLeg}(CombinedHS{FastLeg}(i))-COM(CombinedHS{FastLeg}(i));
    end
    
    %Ankle location relative to COM location at fast leg toe off
    for i = 1:NumStrikes(FastLeg)
        Xf_FTO1(i) = CombinedAnk{FastLeg}(CombinedTO{FastLeg}(i))-COM(CombinedTO{FastLeg}(i));
    end
    
    %Step Length
    SL_s = Xs_SHS2-Xf_SHS2; %Step length slow leg (ankle markers position at slow leg heel strike)
    SL_f = Xf_FHS-Xs_FHS; %Step length fast leg (ankle markers position at fast leg heel strike)
    
    %Step Length Asym Normalized (percentage relative to stride length): SL_f - SL_s / SL_f + SL_s
    StepLengthAsymNorm = zeros(1,(2*(length(SL_s)))-2);
    for i = 1:length(SL_s)-1
         StepLengthAsymNorm(2*i-1) = (SL_f(i)-SL_s(i))/(SL_f(i)+SL_s(i));
         StepLengthAsymNorm(2*i) = (SL_f(i)-SL_s(i+1))/(SL_f(i)+SL_s(i+1));
    end
    
    %Step Length Asym Raw: SL_f - SL_s 
    StepLengthAsymRaw = zeros(1,(2*(length(SL_s)))-2);
    for i = 1:length(SL_s)-1
         StepLengthAsymRaw(2*i-1) = (SL_f(i)-SL_s(i));
         StepLengthAsymRaw(2*i) = (SL_f(i)-SL_s(i+1));
    end

    %Step Times
    for i = 1:NumStrikes(~FastLeg)-1
        T_s(i) = (CombinedHS{~FastLeg}(i+1)-CombinedHS{FastLeg}(i)); % time of slow HS2 - time of fast HS
        T_s(i) = cast(T_s(i),"double")/Rate; % change variable type for precision
        T_f(i) = (CombinedHS{FastLeg}(i)-CombinedHS{~FastLeg}(i)); % time of fast HS - time of slow HS1
        T_f(i) = cast(T_f(i),"double")/Rate; % change variable type for precision
    end
    
    % Step Timing Asymmetry (percentage relative to stride time)
    StepTimeAsymNorm = zeros(1,(2*(length(T_s)))-2);
    for i = 1:length(T_s)-1
        StepTimeAsymNorm(2*i-1) = (T_f(i)-T_s(i))/(T_f(i)+T_s(i));
        StepTimeAsymNorm(2*i) = (T_f(i)-T_s(i+1))/(T_f(i)+T_s(i+1));
    end

    % Average foot velocity with respect to COM while during stance phase
    V_s = (Xs_SHS1-Xs_FHS)./T_s;
    V_f = (Xf_FHS-Xf_SHS2)./T_f;
    
    % Step Velocity Asymmetry
    StepVelocityAsym = zeros(1,(2*(length(T_s)))-2);
    for i = 1:length(T_s)-1
        StepVelocityAsym(2*i-1) = (V_f(i)-V_s(i))/(V_f(i)+V_s(i));
        StepVelocityAsym(2*i) = (V_f(i)-V_s(i+1))/(V_f(i)+V_s(i+1));
    end

    % Step Position Asymmetry (also step position contribution to SlA, (percentage relative to stride) 
    % Positive step position contributions indicate that the unaffected limb progresses further in front of the trunk than the affected limb.
    Alpha_s = Xs_SHS2-Xf_FHS;
    Alpha_f = Xf_FHS-Xs_SHS1;
    StepPositionContrib = zeros(1,(2*(length(Xf_FHS)))-2);
    for i = 1:length(Xf_FHS)-1
        StepPositionContrib(2*i-1) = (Alpha_f(i)-Alpha_s(i));
        StepPositionContrib(2*i) = (Alpha_f(i)-Alpha_s(i+1));
    end

    % Step Timing Contribution to SLA 
    % Negative step time contributions represent shorter stance-phase durations in the more affected versus less affected limb.
    StepTimeContrib = zeros(1,(2*(length(T_s)))-2);
    for i = 1:length(T_s)-1
        StepTimeContrib(2*i-1) = ((V_s(i)+V_f(i))/2)*(T_s(i)-T_f(i));
        StepTimeContrib(2*i) = ((V_s(i+1)+V_f(i))/2)*(T_s(i+1)-T_f(i));
    end
    
    % Step Velocity Contribution to SLA 
    % % Negative step velocity contributions are indicative of greater angular velocity of the less affected limb, and this strategy would be expected from individuals with marked weakness of the paretic plantar flexor
    StepVelocityContrib = zeros(1,(2*(length(T_s)))-2);
    for i = 1:length(T_s)-1
        StepVelocityContrib(2*i-1) = ((T_s(i)+T_f(i))/2)*(V_s(i)-V_f(i));
        StepVelocityContrib(2*i) = ((T_s(i+1)+T_f(i))/2)*(V_s(i+1)-V_f(i));
    end
      
    StepNum = (1:length(StepLengthAsymRaw));
    
    % Loop Through Each Timepoint for SLA Contributtions
    TimepointData.MA_StepLength = TimepointDataGrab(SL_f,Phase,Timepoints,'Step');
    TimepointData.LA_StepLength = TimepointDataGrab(SL_s,Phase,Timepoints,'Step');
    TimepointData.StepLengthAsym_Normalized = TimepointDataGrab(StepLengthAsymNorm,Phase,Timepoints,'Step');
    TimepointData.StepLengthAsym_Raw = TimepointDataGrab(StepLengthAsymRaw,Phase,Timepoints,'Step');
    TimepointData.StepTimeAsym = TimepointDataGrab(StepTimeAsymNorm,Phase,Timepoints,'Step');
    TimepointData.StepPositionContrib = TimepointDataGrab(StepPositionContrib,Phase,Timepoints,'Step');
    TimepointData.StepTimeContrib = TimepointDataGrab(StepTimeContrib,Phase,Timepoints,'Step');
    TimepointData.StepVelocityContrib = TimepointDataGrab(StepVelocityContrib,Phase,Timepoints,'Step');

    %SLA Table 
    % SLA_Contrib_Results= table(StepNum', StepLengthAsymNorm', StepLengthAsymRaw',StepPositionContrib', StepTimeContrib', StepVelocityContrib');
    % SLA_Contrib_Results.Properties.VariableNames = {'StepNum', 'StepLengthAsymNorm','StepLengthAsymRaw', 'StepPositionContrib', 'StepTimeContrib', 'StepVelocityContrib'};

%% Propulsion and Braking Asymmetry - relative to body weight 
    % X and Y force data are flipped 
     % Initialize with Nans, Set WindowLength to the length of the first gait cycle - for later interpolation
        MA_MaxWindowLength = length(MA_Force(1:3,(MA_HeelStrikes(1)*10):(MA_ToeOffs(2)*10))); % Need to multiple by 10 so index is correct (force collects at 1000hz)
        MA_XForce_BW = nan((length(MA_HeelStrikes)-1), MA_MaxWindowLength);
        MA_YForce_BW = nan((length(MA_HeelStrikes)-1), MA_MaxWindowLength);
        MA_ZForce_BW = nan((length(MA_HeelStrikes)-1), MA_MaxWindowLength);
        LA_MaxWindowLength = length(LA_Force(1:3,(LA_HeelStrikes(1)*10):(LA_ToeOffs(2)*10))); % Need to multiple by 10 so index is correct (force collects at 1000hz)
        LA_XForce_BW = nan((length(LA_HeelStrikes)-1), LA_MaxWindowLength);
        LA_YForce_BW = nan((length(LA_HeelStrikes)-1), LA_MaxWindowLength);
        LA_ZForce_BW = nan((length(LA_HeelStrikes)-1), LA_MaxWindowLength);

   for i = 1:(length(MA_HeelStrikes)-1)
       TempWindow = MA_Force(1:3,(10*(MA_HeelStrikes(i))):(10*(MA_ToeOffs(i+1)))); % Find window of forces from each stance phase (heel strike to toe off) - 10x due to collection at 1000hz instead of 100hz (markers)
       MA_XForce_BW(i,1:length(TempWindow)) = (TempWindow(1,:))/BodyWeight; % Propoulsion/braking forces
       MA_YForce_BW(i,1:length(TempWindow)) = (TempWindow(2,:))/BodyWeight; % Mediolateral force
       MA_ZForce_BW(i,1:length(TempWindow)) = (TempWindow(3,:))/BodyWeight; % Vertical forces
       MA_PeakPropulsion_BW(i) = min(TempWindow(1,:))/BodyWeight;
       MA_PeakBraking_BW(i) = max(TempWindow(1,:))/BodyWeight;
       MA_PeakVForce_BW(i) = max(TempWindow(3,:))/BodyWeight;
   end
   for i = 1:(length(LA_HeelStrikes)-1)
       TempWindow = LA_Force(1:3,(10*(LA_HeelStrikes(i))):(10*(LA_ToeOffs(i+1)))); % Find window of forces from each stance phase (heel strike to toe off) - 10x due to collection at 1000hz instead of 100hz (markers)
       LA_XForce_BW(i,1:length(TempWindow)) = (TempWindow(1,:))/BodyWeight; % Propulsion/braking forces
       LA_YForce_BW(i,1:length(TempWindow)) = (TempWindow(2,:))/BodyWeight; % Mediolateral force
       LA_ZForce_BW(i,1:length(TempWindow)) = (TempWindow(3,:))/BodyWeight; % Vertical forces
       LA_PeakPropulsion_BW(i) = min(TempWindow(1,:))/BodyWeight;
       LA_PeakBraking_BW(i) = max(TempWindow(1,:))/BodyWeight;
       LA_PeakVForce_BW(i) = max(TempWindow(3,:))/BodyWeight;
   end

   PeakPropulsionAsym = zeros(1,(2*(length(MA_PeakPropulsion_BW)))-2);
   PeakBrakingAsym = zeros(1,(2*(length(MA_PeakBraking_BW)))-2);
   PeakVForceAsym = zeros(1,(2*(length(MA_PeakVForce_BW)))-2);
    for i = 1:length(LA_PeakPropulsion_BW)-1
        PeakPropulsionAsym(2*i-1) = (MA_PeakPropulsion_BW(i) - LA_PeakPropulsion_BW(i));
        PeakPropulsionAsym(2*i) = (MA_PeakPropulsion_BW(i) - LA_PeakPropulsion_BW(i+1));
        PeakBrakingAsym(2*i-1) = MA_PeakBraking_BW(i) - LA_PeakBraking_BW(i);
        PeakBrakingAsym(2*i) = MA_PeakBraking_BW(i) - LA_PeakBraking_BW(i+1);
        PeakVForceAsym(2*i-1) = MA_PeakVForce_BW(i) - LA_PeakVForce_BW(i);
        PeakVForceAsym(2*i) = MA_PeakVForce_BW(i) - LA_PeakVForce_BW(i+1);
    end

   % Loop through timepoints
   TimepointData.PeakPropulsionAsym = TimepointDataGrab(PeakPropulsionAsym,Phase,Timepoints,'Step');
   TimepointData.PeakBrakingAsym = TimepointDataGrab(PeakBrakingAsym,Phase,Timepoints,'Step');
   TimepointData.PeakVForceAsym = TimepointDataGrab(PeakVForceAsym,Phase,Timepoints,'Step');

%% Magnitude of Adaptation (SLA, DblSup, PropAsym): late (de)adapt - early (de)adapt
    if strcmp(Phase,'Adaptation')
        SLAAdaptMagnitude = TimepointData.StepLengthAsym_Normalized.LAd - TimepointData.StepLengthAsym_Normalized.EAd;
        PropAdaptMagnitude = TimepointData.PeakPropulsionAsym.LAd - TimepointData.PeakPropulsionAsym.EAd;
        DblSupAdaptMagnitude = TimepointData.DoubleSupportAsym.LAd - TimepointData.DoubleSupportAsym.EAd;
        TimepointData.AdaptationMagnitudes = struct('StepLengthAsym', SLAAdaptMagnitude, 'PeakPropulsionAsym', PropAdaptMagnitude, 'DoubleSupportAsym', DblSupAdaptMagnitude);
    elseif strcmp(Phase,'Deadaptation')
        SLAAdaptMagnitude = TimepointData.StepLengthAsym_Normalized.LDe - TimepointData.StepLengthAsym_Normalized.EDe;
        PropAdaptMagnitude = TimepointData.PeakPropulsionAsym.LDe - TimepointData.PeakPropulsionAsym.EDe;   
        DblSupAdaptMagnitude = TimepointData.DoubleSupportAsym.LDe - TimepointData.DoubleSupportAsym.EDe;
        TimepointData.AdaptationMagnitudes = struct('StepLengthAsym', SLAAdaptMagnitude, 'PeakPropulsionAsym', PropAdaptMagnitude, 'DoubleSupportAsym', DblSupAdaptMagnitude);
    end

    
%% Rate of Adaptation and Magnitude of Adaptation (SLA Dblsup and PropAsym): Exponential fitting 
    %  Rate of adaptation is also a measure of sensorimotor integration: The rate that SLA is restored represents the rate at which the nervous system integrates sensory feedback and updates motor output - Kuhman et al. 2022 IBRO Neuroscience Reports
    %  This fits a single exponential model to SLA timeseries using the particle swarm algorithm - Rashid et al. 2020. Brain Sciences.

    % Parameters of the fitted model(FitParams): [a_s, b_s, a_f, b_f, c] 

    % In this usage the parameters of the double exponential model represent:
    % The initial asymmetry: (a_s + a_f + c)
    % Change in symmetry (magnitude) from beginning to the end of the training phase: (a_s + a_f)
    % The number of steps taken to complete 50% of the change: (ln(2)/abs(b_s))
    % The asymmetry at the end of the training phase: (c)

     [DblSup_Fitted, DblSup_Resids, DblSup_FitParams, DblSup_LinearCIs, DblSup_Aic] = fitSingleExpPSO(DblSupAsym(1:end));  % Depending on data, cropping the first few steps may improve fit
     [Prop_Fitted, Prop_Resids, Prop_FitParams, Prop_LinearCIs, Prop_Aic] = fitSingleExpPSO(PeakPropulsionAsym(1:end));  
     [SLA_Fitted, SLA_Resids, SLA_FitParams, SLA_LinearCIs, SLA_Aic] = fitSingleExpPSO(StepLengthAsymNorm(1:end));  

     if strcmp(Phase,'Adaptation') || strcmp(Phase,'Deadaptation')
        % Calculate Magnitude of Adaptation 
        SLAAdaptMagnitude = SLA_FitParams(1)+ SLA_FitParams(3);
        PropAdaptMagnitude = Prop_FitParams(1)+ Prop_FitParams(3);
        DblSupAdaptMagnitude = DblSup_FitParams(1)+ DblSup_FitParams(3);
        TimepointData.ModelFit_AdaptationMagnitudes = struct('StepLengthAsym', SLAAdaptMagnitude, 'PeakPropulsionAsym', PropAdaptMagnitude, 'DoubleSupportAsym', DblSupAdaptMagnitude);
        
        % Calculate Rate of Adaptation (how many steps it takes to get to 50% of total adaptation)
        SLAAdaptRate = log(2)/abs(SLA_FitParams(2));
        PropAdaptRate = log(2)/abs(Prop_FitParams(2));
        DblSupAdaptRate = log(2)/abs(DblSup_FitParams(2));
        TimepointData.ModelFit_AdaptationRates = struct('StepLengthAsym', SLAAdaptRate, 'PeakPropulsionAsym', PropAdaptRate, 'DoubleSupportAsym', DblSupAdaptRate);
     end
% 
%% Create Data Structures
ExportStruct = inputdlg('Save and Export Data to Structures? Y/N: ');
 if strcmp(ExportStruct,'Y')
   % Make sure SBNT_Data and SBNT_Data_Interp are in the workspace if you want to add on another participant, otherwise it will create a new structure
    SBNT_Data.(SubNum).(Visit).(Phase).MA_AnkleAngles = CreateDataStruct(MA_AnkleAng,MA_HeelStrikes,MA_ToeOffs);
    SBNT_Data.(SubNum).(Visit).(Phase).LA_AnkleAngles = CreateDataStruct(LA_AnkleAng,LA_HeelStrikes,LA_ToeOffs); 
    SBNT_Data.(SubNum).(Visit).(Phase).MA_KneeAngles = CreateDataStruct(MA_KneeAng,MA_HeelStrikes, MA_ToeOffs); 
    SBNT_Data.(SubNum).(Visit).(Phase).LA_KneeAngles = CreateDataStruct(LA_KneeAng,LA_HeelStrikes,LA_ToeOffs); 
    SBNT_Data.(SubNum).(Visit).(Phase).MA_HipAngles = CreateDataStruct(MA_HipAng,MA_HeelStrikes,MA_ToeOffs); 
    SBNT_Data.(SubNum).(Visit).(Phase).LA_HipAngles = CreateDataStruct(LA_HipAng,LA_HeelStrikes,LA_ToeOffs); 
    SBNT_Data.(SubNum).(Visit).(Phase).MA_AnkleMoments = CreateDataStruct(MA_AnkleMom,MA_HeelStrikes,MA_ToeOffs); 
    SBNT_Data.(SubNum).(Visit).(Phase).LA_AnkleMoments = CreateDataStruct(LA_AnkleMom,LA_HeelStrikes,LA_ToeOffs);
    SBNT_Data.(SubNum).(Visit).(Phase).MA_KneeMoments = CreateDataStruct(MA_KneeMom,MA_HeelStrikes,MA_ToeOffs);
    SBNT_Data.(SubNum).(Visit).(Phase).LA_KneeMoments = CreateDataStruct(LA_KneeMom,LA_HeelStrikes,LA_ToeOffs);
    SBNT_Data.(SubNum).(Visit).(Phase).MA_HipMoments = CreateDataStruct(MA_HipMom,MA_HeelStrikes,MA_ToeOffs);
    SBNT_Data.(SubNum).(Visit).(Phase).LA_HipMoments = CreateDataStruct(LA_HipMom,LA_HeelStrikes,LA_ToeOffs);
    SBNT_Data.(SubNum).(Visit).(Phase).MA_AnklePowers = CreateDataStruct(MA_AnklePow,MA_HeelStrikes,MA_ToeOffs);
    SBNT_Data.(SubNum).(Visit).(Phase).LA_AnklePowers = CreateDataStruct(LA_AnklePow,LA_HeelStrikes,LA_ToeOffs);
    SBNT_Data.(SubNum).(Visit).(Phase).MA_KneePowers = CreateDataStruct(MA_KneePow,MA_HeelStrikes,MA_ToeOffs);
    SBNT_Data.(SubNum).(Visit).(Phase).LA_KneePowers = CreateDataStruct(LA_KneePow,LA_HeelStrikes,LA_ToeOffs);
    SBNT_Data.(SubNum).(Visit).(Phase).MA_HipPowers = CreateDataStruct(MA_HipPow,MA_HeelStrikes,MA_ToeOffs);
    SBNT_Data.(SubNum).(Visit).(Phase).LA_HipPowers = CreateDataStruct(LA_HipPow,LA_HeelStrikes,LA_ToeOffs);
    SBNT_Data.(SubNum).(Visit).(Phase).MA_StanceGRF = struct('X',MA_XForce_BW,'Y',MA_YForce_BW,'Z',MA_ZForce_BW);
    SBNT_Data.(SubNum).(Visit).(Phase).LA_StanceGRF = struct('X',LA_XForce_BW,'Y',LA_YForce_BW,'Z',LA_ZForce_BW);
    SBNT_Data.(SubNum).(Visit).(Phase).Phi = Phi;
    SBNT_Data.(SubNum).(Visit).(Phase).MA_StepLength = SL_f;
    SBNT_Data.(SubNum).(Visit).(Phase).LA_StepLength = SL_s;
    SBNT_Data.(SubNum).(Visit).(Phase).StepLengthAsym_Norm = StepLengthAsymNorm;
    SBNT_Data.(SubNum).(Visit).(Phase).StepLengthAsym_Raw = StepLengthAsymRaw;
    SBNT_Data.(SubNum).(Visit).(Phase).StepLengthAsym_Fitted = SLA_Fitted;
    SBNT_Data.(SubNum).(Visit).(Phase).StepTimeAsym_Norm = StepTimeAsymNorm;
    SBNT_Data.(SubNum).(Visit).(Phase).StepPositionContribution = StepPositionContrib;
    SBNT_Data.(SubNum).(Visit).(Phase).StepTimeContribution = StepTimeContrib;
    SBNT_Data.(SubNum).(Visit).(Phase).StepVelocityContribution = StepVelocityContrib;
    SBNT_Data.(SubNum).(Visit).(Phase).DoubleSupportAsym = DblSupAsym;
    SBNT_Data.(SubNum).(Visit).(Phase).DoubleSupportAsym_Fitted = DblSup_Fitted;
    SBNT_Data.(SubNum).(Visit).(Phase).PeakPropulsionAsym = PeakPropulsionAsym;
    SBNT_Data.(SubNum).(Visit).(Phase).PeakPropulsionAsym_Fitted = Prop_Fitted;
    SBNT_Data.(SubNum).(Visit).(Phase).PeakBrakingAsym = PeakBrakingAsym;
    SBNT_Data.(SubNum).(Visit).(Phase).PeakVerticalForceAsym = PeakVForceAsym;

    % Create Interpolated Structure 
        SBNT_Data_Interp.(SubNum).(Visit).(Phase) = InterpolateNestedStruct((SBNT_Data.(SubNum).(Visit).(Phase)), (Phase), 100);

    % Add TimepointData to Structures
        SBNT_Data.(SubNum).(Visit).(Phase).TimepointData = TimepointData;
        SBNT_Data_Interp.(SubNum).(Visit).(Phase).TimepointData = TimepointData;
 end
%% Adaptation Data Table
% Create variable names for each data point
VariableNames = cell(1, (length(SBNT_Data_Interp.(SubNum).(Visit).(Phase).StepLengthAsym_Norm(:))));
for i = 1:length(SBNT_Data_Interp.(SubNum).(Visit).(Phase).StepLengthAsym_Norm(:))
    VariableNames{i} = ['GC_' num2str(i)];
end

SLAAdaptationTable = array2table(SBNT_Data_Interp.(SubNum).(Visit).(Phase).StepLengthAsym_Norm(:)', 'VariableNames', VariableNames);
DblSupAdaptationTable = array2table(SBNT_Data_Interp.(SubNum).(Visit).(Phase).DoubleSupportAsym(:)','VariableNames', VariableNames);
PropAdaptationTable = array2table(SBNT_Data_Interp.(SubNum).(Visit).(Phase).PeakPropulsionAsym(:)','VariableNames' ,VariableNames);
PhiAdaptationTable = array2table(SBNT_Data_Interp.(SubNum).(Visit).(Phase).Phi(:)','VariableNames',VariableNames);
SLAFittedAdaptationTable = array2table(SBNT_Data_Interp.(SubNum).(Visit).(Phase).StepLengthAsym_Fitted(:)','VariableNames',VariableNames);

% Add 'Participant' and 'Visit' columns
SLAAdaptationTable = addvars(SLAAdaptationTable, repmat({SubNum}, height(SLAAdaptationTable), 1), repmat({Visit}, height(SLAAdaptationTable), 1), 'Before', 1, 'NewVariableNames', {'Participant', 'Visit'});
DblSupAdaptationTable = addvars(DblSupAdaptationTable, repmat({SubNum}, height(DblSupAdaptationTable), 1), repmat({Visit}, height(DblSupAdaptationTable), 1), 'Before', 1, 'NewVariableNames', {'Participant', 'Visit'});
PropAdaptationTable = addvars(PropAdaptationTable, repmat({SubNum}, height(PropAdaptationTable), 1), repmat({Visit}, height(PropAdaptationTable), 1), 'Before', 1, 'NewVariableNames', {'Participant', 'Visit'});
PhiAdaptationTable = addvars(PhiAdaptationTable, repmat({SubNum}, height(PhiAdaptationTable), 1), repmat({Visit}, height(PhiAdaptationTable), 1), 'Before', 1, 'NewVariableNames', {'Participant', 'Visit'});
SLAFittedAdaptationTable = addvars(SLAFittedAdaptationTable, repmat({SubNum}, height(SLAFittedAdaptationTable), 1), repmat({Visit}, height(SLAFittedAdaptationTable), 1), 'Before', 1, 'NewVariableNames', {'Participant', 'Visit'});

%% Make Timepoint Data Table 
% Initialize an empty table
TimepointDataTable = {(SubNum),(Visit),(Alimb)};
TimepointDataTable = cell2table(TimepointDataTable,"VariableNames",{'Participant','Visit','More_Affected_Limb'});
% Loop through each field of the original structure
Fields = fieldnames(TimepointData);
for i = 1:length(Fields)
    % Extract the nested structure
    SubStruct = TimepointData.(Fields{i});
    SubFields = fieldnames(SubStruct);
    
    % Initialize cell arrays to store values and variable names
    Values = cell(1, length(SubFields));
    VarNames = cell(1, length(SubFields));

    % Loop through the SubFields and extract values and create variable names
    for j = 1:length(SubFields)
        Values{j} = SubStruct.(SubFields{j});
        VarNames{j} = [Fields{i}, '_', SubFields{j}];
    end
    
    NewVariable = table(Values{:}, 'VariableNames', VarNames);
    TimepointDataTable = [TimepointDataTable, NewVariable];
end

%% Export Trial to Group Data
% Data exported in this section is added onto existing data tables

ExportData = inputdlg('Save and Export Data to Group xlsx? Y/N: ');
 if strcmp(ExportData,'Y')
    % TimepointData
      if strcmp(Phase,'Baseline')
        TimepointExcelData = readtable([RootPath '\Timepoint_Data.xlsx'], 'Sheet', 'Baseline');
        % Append new data to the existing table
        TimepointExcelData = [TimepointExcelData; TimepointDataTable];
        % Save the updated data to the same sheet
        writetable(TimepointExcelData, [RootPath '\Timepoint_Data.xlsx'], 'Sheet', 'Baseline');

      elseif strcmp(Phase,'Adaptation')
        TimepointExcelData = readtable([RootPath '\Timepoint_Data.xlsx'], 'Sheet', 'Adaptation');
        % Append new data to the existing table
        TimepointExcelData = [TimepointExcelData; TimepointDataTable];
        % Save the updated data to the same sheet
        writetable(TimepointExcelData, [RootPath '\Timepoint_Data.xlsx'], 'Sheet', 'Adaptation');

      elseif strcmp(Phase,'Deadaptation')
        TimepointExcelData = readtable([RootPath '\Timepoint_Data.xlsx'], 'Sheet', 'Deadaptation');
        % Append new data to the existing table
        TimepointExcelData = [TimepointExcelData; TimepointDataTable];
        % Save the updated data to the same sheet
        writetable(TimepointExcelData, [RootPath '\Timepoint_Data.xlsx'], 'Sheet', 'Deadaptation');
       end

    %Adaptation Profiles
    if strcmp(Phase,'Baseline')
         SLAExcelData = readtable([RootPath '\StepLengthAsym_Adaptation_Data.xlsx'], 'Sheet', 'Baseline');
         DblSupExcelData = readtable([RootPath '\DoubleSupportAsym_Adaptation_Data.xlsx'], 'Sheet', 'Baseline');
         PropExcelData = readtable([RootPath '\PeakPropulsionAsym_Adaptation_Data.xlsx'], 'Sheet', 'Baseline');
         PhiExcelData = readtable([RootPath '\Phi_Adaptation_Data.xlsx'], 'Sheet', 'Baseline');
         SLAFittedExcelData = readtable([RootPath '\StepLengthAsym_Fitted_Adaptation_Data.xlsx'], 'Sheet', 'Baseline');
         % Append new data to the existing table
           SLAExcelData = [SLAExcelData; SLAAdaptationTable];
           DblSupExcelData = [DblSupExcelData; DblSupAdaptationTable];
           PropExcelData = [PropExcelData; PropAdaptationTable];
           PhiExcelData = [PhiExcelData; PhiAdaptationTable];
           SLAFittedExcelData = [SLAFittedExcelData; SLAFittedAdaptationTable];
         % Save the updated data to the same sheet
           writetable(SLAExcelData, [RootPath '\StepLengthAsym_Adaptation_Data.xlsx'], 'Sheet', 'Baseline', 'WriteMode', 'overwrite');
           writetable(DblSupExcelData, [RootPath '\DoubleSupportAsym_Adaptation_Data.xlsx'], 'Sheet', 'Baseline', 'WriteMode', 'overwrite');
           writetable(PropExcelData, [RootPath '\PeakPropulsionAsym_Adaptation_Data.xlsx'], 'Sheet', 'Baseline', 'WriteMode', 'overwrite');
           writetable(PhiExcelData, [RootPath '\Phi_Adaptation_Data.xlsx'], 'Sheet', 'Baseline', 'WriteMode', 'overwrite')
           writetable(SLAFittedExcelData, [RootPath '\StepLengthAsym_Fitted_Adaptation_Data.xlsx'], 'Sheet', 'Baseline', 'WriteMode', 'overwrite')

    elseif strcmp(Phase,'Adaptation')
         SLAExcelData = readtable([RootPath '\StepLengthAsym_Adaptation_Data.xlsx'], 'Sheet', 'Adaptation');
         DblSupExcelData = readtable([RootPath '\DoubleSupportAsym_Adaptation_Data.xlsx'], 'Sheet', 'Adaptation');
         PropExcelData = readtable([RootPath '\PeakPropulsionAsym_Adaptation_Data.xlsx'], 'Sheet', 'Adaptation');
         PhiExcelData = readtable([RootPath '\Phi_Adaptation_Data.xlsx'], 'Sheet', 'Adaptation');
         SLAFittedExcelData = readtable([RootPath '\StepLengthAsym_Fitted_Adaptation_Data.xlsx'], 'Sheet', 'Adaptation');
         % Append new data to the existing table
           SLAExcelData = [SLAExcelData; SLAAdaptationTable];
           DblSupExcelData = [DblSupExcelData; DblSupAdaptationTable];
           PropExcelData = [PropExcelData; PropAdaptationTable];
           PhiExcelData = [PhiExcelData; PhiAdaptationTable];
           SLAFittedExcelData = [SLAFittedExcelData; SLAFittedAdaptationTable];
         % Save the updated data to the same sheet
           writetable(SLAExcelData, [RootPath '\StepLengthAsym_Adaptation_Data.xlsx'], 'Sheet', 'Adaptation', 'WriteMode', 'overwrite');
           writetable(DblSupExcelData, [RootPath '\DoubleSupportAsym_Adaptation_Data.xlsx'], 'Sheet', 'Adaptation', 'WriteMode', 'overwrite');
           writetable(PropExcelData, [RootPath '\PeakPropulsionAsym_Adaptation_Data.xlsx'], 'Sheet', 'Adaptation', 'WriteMode', 'overwrite');
           writetable(PhiExcelData, [RootPath '\Phi_Adaptation_Data.xlsx'], 'Sheet', 'Adaptation', 'WriteMode', 'overwrite')
           writetable(SLAFittedExcelData, [RootPath '\StepLengthAsym_Fitted_Adaptation_Data.xlsx'], 'Sheet', 'Adaptation', 'WriteMode', 'overwrite')


    elseif strcmp(Phase,'Deadaptation')
        SLAExcelData = readtable([RootPath '\StepLengthAsym_Adaptation_Data.xlsx'], 'Sheet', 'Deadaptation');
         DblSupExcelData = readtable([RootPath '\DoubleSupportAsym_Adaptation_Data.xlsx'], 'Sheet', 'Deadaptation');
         PropExcelData = readtable([RootPath '\PeakPropulsionAsym_Adaptation_Data.xlsx'], 'Sheet', 'Deadaptation');
         PhiExcelData = readtable([RootPath '\Phi_Adaptation_Data.xlsx'], 'Sheet', 'Deadaptation');
         SLAFittedExcelData = readtable([RootPath '\StepLengthAsym_Fitted_Adaptation_Data.xlsx'], 'Sheet', 'Deadaptation');
         % Append new data to the existing table
           SLAExcelData = [SLAExcelData; SLAAdaptationTable];
           DblSupExcelData = [DblSupExcelData; DblSupAdaptationTable];
           PropExcelData = [PropExcelData; PropAdaptationTable];
           PhiExcelData = [PhiExcelData; PhiAdaptationTable];
           SLAFittedExcelData = [SLAFittedExcelData; SLAFittedAdaptationTable];
         % Save the updated data to the same sheet
           writetable(SLAExcelData, [RootPath '\StepLengthAsym_Adaptation_Data.xlsx'], 'Sheet', 'Deadaptation', 'WriteMode', 'overwrite');
           writetable(DblSupExcelData, [RootPath '\DoubleSupportAsym_Adaptation_Data.xlsx'], 'Sheet', 'Deadaptation', 'WriteMode', 'overwrite');
           writetable(PropExcelData, [RootPath '\PeakPropulsionAsym_Adaptation_Data.xlsx'], 'Sheet', 'Deadaptation', 'WriteMode', 'overwrite');
           writetable(PhiExcelData, [RootPath '\Phi_Adaptation_Data.xlsx'], 'Sheet', 'Deadaptation', 'WriteMode', 'overwrite')
           writetable(SLAFittedExcelData, [RootPath '\StepLengthAsym_Fitted_Adaptation_Data.xlsx'], 'Sheet', 'Deadaptation', 'WriteMode', 'overwrite')

    end
    disp('Exported Successfully :)')
 end

% %% Plot Raw SLA, Position Contribution, Time Contribution and Velocity Contribution
% SLA_Plot = figure(2);
% width=1500;
% height=1000;
% set(gcf,'position',[0,0,width,height])
% hold on
% set(gca,'FontSize',14);
% ylabel('Asymmetry (Absolute Difference in mm)','FontName','Arial')
% xlabel('Step Number','FontName','Arial')
% title('Step Length Asymmetry Contributions','FontName','Arial')
% 
% plot(SLA_Contrib_Results{:,1},SLA_Contrib_Results{:,3},'b','LineWidth',2.0)
% plot(SLA_Contrib_Results{:,1},SLA_Contrib_Results{:,4},'r','LineWidth',0.5)
% plot(SLA_Contrib_Results{:,1},SLA_Contrib_Results{:,5},'m','LineWidth',0.5)
% plot(SLA_Contrib_Results{:,1},SLA_Contrib_Results{:,6},'g','LineWidth',0.5)
% legend('Step Length Asymmetry','Spatial Component','Step Time Component','Step Velocity Component')
% 
% 
%     if strcmp(ExportFig,'Y')
%       PlotSLA = 'SLA Figure' ;
%         savefig(SLA_Plot,PlotSLA);
%         print([RootPath, '\', TrialName, ' SLA Figure'],'-dpng','-r0')
%     else
%     end


%% Functions

% function to grab timepoint data of a variable 
function TimepointData = TimepointDataGrab(Variable,Phase,Timepoints,DataRate) %Data rate is used to confirm the frequency of data (1 point per step or 1 point for stride) E.g. SLA 1 point for every step while PCI has one point for every stride
TimepointData = struct();
    % Loop through each timepoint
    for i = 1:length(Timepoints)
        CurrentTimepoint = Timepoints{i};
        TimepointData.(CurrentTimepoint) = struct();

        if strcmp(Phase, 'Baseline') 
            if strcmp(DataRate, 'Stride')
            Variable_B = Variable(end-29:end);
            elseif strcmp(DataRate, 'Step')
            Variable_B = Variable(end-59:end);
            end
            TimepointData.(CurrentTimepoint) = mean(Variable_B);
        
        elseif strcmp(Phase, 'Adaptation') 
            if strcmp(CurrentTimepoint, 'EAd')
                if strcmp(DataRate, 'Stride')
                Variable_EAd = Variable(6:30);
                elseif strcmp(DataRate, 'Step')
                Variable_EAd = Variable(12:60);
                end
                TimepointData.(CurrentTimepoint) = mean(Variable_EAd);

            elseif strcmp(CurrentTimepoint, 'MAd')
               if strcmp(DataRate, 'Stride')
               Variable_MAd = Variable(fix(length(Variable)/2) - 14 : fix(length(Variable)/2) + 15);
               elseif strcmp(DataRate, 'Step')
               Variable_MAd = Variable(fix(length(Variable)/2) - 29 : fix(length(Variable)/2) + 30); 
               end
               TimepointData.(CurrentTimepoint) = mean(Variable_MAd);
                
            elseif strcmp(CurrentTimepoint, 'LAd')
                if strcmp(DataRate, 'Stride')
                Variable_LAd = Variable(end-29:end);
                elseif strcmp(DataRate, 'Step')
                Variable_LAd = Variable(end-59:end); 
                end
                TimepointData.(CurrentTimepoint) = mean(Variable_LAd);
            end  
            
        elseif strcmp(Phase, 'Deadaptation') 
            if strcmp(CurrentTimepoint, 'EDe')
                if strcmp(DataRate, 'Stride')
                Variable_EDe = Variable(6:30);
                elseif strcmp(DataRate, 'Step')
                Variable_EDe = Variable(12:60);
                end
                TimepointData.(CurrentTimepoint) = mean(Variable_EDe);

            elseif strcmp(CurrentTimepoint, 'MDe')
               if strcmp(DataRate, 'Stride')
               Variable_MDe = Variable(fix(length(Variable)/2) - 14 : fix(length(Variable)/2) + 15);
               elseif strcmp(DataRate, 'Step')
               Variable_MDe = Variable(fix(length(Variable)/2) - 29 : fix(length(Variable)/2) + 30); 
               end
               TimepointData.(CurrentTimepoint) = mean(Variable_MDe);
                
            elseif strcmp(CurrentTimepoint, 'LDe')
                if strcmp(DataRate, 'Stride')
                Variable_LDe = Variable(end-29:end);
                elseif strcmp(DataRate, 'Step')
                Variable_LDe = Variable(end-59:end); 
                end
                TimepointData.(CurrentTimepoint) = mean(Variable_LDe);
            end  
        end
    end
end

% example usage 
% TimepointData.MA_StepDur = TimepointDataGrab(MA_StepDur,Phase,Timepoints,'Stride');

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
    %   d. linearCIs: CIs based on the assumption of linearization.
    %   e. aic: AIC value for the fitted model.
    %
    %   Toolboxes required: Optimization Toolbox, Global Optimization Toolbox, Particle Swarm Optimization Toolbox - Brian Birge
    %
    %   Copyright (c) <2018> <Usman Rashid>
    %   Licensed under the MIT License. See LICENSE in the project root for
    %   license information.
    %
    %   https://github.com/GallVp/knkTools/tree/master/exponentialModels
    %   Article on Model:  Rashid U, Kumari N, Signal N, Taylor D, Vandal AC. On Nonlinear Regression for Trends in Split-Belt Treadmill Training. Brain Sciences. 2020; 10(10):737. https://doi.org/10.3390/brainsci10100737
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
    
    %% Set up particle swarm algorithm for optimization
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
    
    %% Run optimization twice and chose the solution with the smaller cost
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
    subplot(2, 2, 1)
    plot(1:length(dataSeries), dataSeries, 'k.');
    hold
    plot(1:length(fitted), fitted, 'r-', 'LineWidth', 1.5);
    title('Single exp. model')
    box off;
    ylabel('Symmetry')
    ax = axis;
    axis([1 length(dataSeries) ax(3) ax(4)]);
    xlabel('Stride No.')
    
    subplot(2, 2, 2)
    plot(1:length(dataSeries), resids, 'k.');
    box off;
    ax = axis;
    axis([1 length(dataSeries) ax(3) ax(4)]);
    xlabel('Stride No.')
    
    subplot(2, 2, 3)
    histogram(resids, 10, 'FaceColor', [1 1 1], 'LineWidth', 1.5);
    box off;
    xlabel('Symmetry');
    ylabel('Frequency');
    
    subplot(2, 2, 4)
    qh = qqplot(resids);
    title('');
    ylabel('Quantiles of residuals');
    xlabel('Standard normal quantiles');
    qh(1).MarkerEdgeColor = [0 0 0];

end

%% Create a data structure for a variable (stance, swing, x,y,z)
function DataStruct = CreateDataStruct(Variable, HeelStrikes, ToeOffs)
    GaitCycles = length(HeelStrikes) - 1;
     % Set WindowLength to the length of the first SwingWindow 
            MaxWindowLengthSwing = length(Variable(1:3, ToeOffs(1):HeelStrikes(1)));
            MaxWindowLengthStance = length(Variable(1:3,HeelStrikes(1):ToeOffs(2)));
    
            DataStruct.Swing = struct('X', nan(GaitCycles, MaxWindowLengthSwing), 'Y', nan(GaitCycles, MaxWindowLengthSwing), 'Z', nan(GaitCycles, MaxWindowLengthSwing));
            DataStruct.Stance = struct('X', nan(GaitCycles, MaxWindowLengthStance), 'Y', nan(GaitCycles, MaxWindowLengthStance), 'Z', nan(GaitCycles, MaxWindowLengthStance));

        for i = 1:GaitCycles
            SwingWindow = Variable(1:3, ToeOffs(i):HeelStrikes(i));
            StanceWindow = Variable(1:3, HeelStrikes(i):ToeOffs(i+1));
    
            DataStruct.Swing.X(i, 1:length(SwingWindow)) = SwingWindow(1, :);
            DataStruct.Swing.Y(i, 1:length(SwingWindow)) = SwingWindow(2, :);
            DataStruct.Swing.Z(i, 1:length(SwingWindow)) = SwingWindow(3, :);
    
            DataStruct.Stance.X(i, 1:length(StanceWindow)) = StanceWindow(1, :);
            DataStruct.Stance.Y(i, 1:length(StanceWindow)) = StanceWindow(2, :);
            DataStruct.Stance.Z(i, 1:length(StanceWindow)) = StanceWindow(3, :);
        end
   
end

% Usage example for More Affected Ankle Angles
%MA_AnkleAng = CreateJointStruct(MA_AnkleAng, MA_HeelStrikes, MA_ToeOffs);

%% Function to interpolate nested structure 
function InterpolatedStruct = InterpolateNestedStruct(RawStruct, Phase, MaxColumns)
    FieldNames = fieldnames(RawStruct);
    InterpolatedStruct = struct();

    for i = 1:numel(FieldNames)
        FieldName = FieldNames{i};
        FieldValue = RawStruct.(FieldName);
        
        if isstruct(FieldValue)
            % If the field is a structure, recursively interpolate it
            InterpolatedStruct.(FieldName) = InterpolateNestedStruct(FieldValue, Phase, MaxColumns);
        else
             % If the field is numeric data, interpolate each row individually
            [NumRows, ~] = size(FieldValue);
            
            for Row = 1:NumRows
                CurrentRowData = FieldValue(Row, :);
                NaNIndices = isnan(CurrentRowData);
                NonNaNIndices = ~NaNIndices;
    
                 if contains(FieldName, 'X') || contains(FieldName, 'Y') || contains(FieldName, 'Z')
                        if sum(NonNaNIndices) > 0
                            % Interpolate or extrapolate only the non-NaN values
                            X = find(NonNaNIndices);
                            XInterpolated = linspace(X(1), X(end), MaxColumns);
                            InterpolatedRowData = interp1(X, CurrentRowData(NonNaNIndices), XInterpolated, 'linear', 'extrap');
                            InterpolatedData(Row, :) = InterpolatedRowData;
                        else
                            % If all values are NaN, keep them as NaN
                            InterpolatedData(Row, :) = NaN;
                        end
                 else % Interpolated calculated variables to number of data points based on phase
                     if strcmp(Phase,'Baseline')
                            DataPoints = 100;
                     elseif strcmp(Phase,'Adaptation') || strcmp(Phase,'Deadaptation') 
                            DataPoints = 500;
                     end

                        X = find(NonNaNIndices);
                            XInterpolated = linspace(X(1), X(end), DataPoints);
                            InterpolatedRowData = interp1(X, CurrentRowData(NonNaNIndices), XInterpolated, 'linear', 'extrap');
                            InterpolatedData(Row, :) = InterpolatedRowData;
                 end

            end
            InterpolatedStruct.(FieldName) = InterpolatedData;
        end
    end
end

%Example usage         
% SBNT_Data_Interp.(SubNum).(Phase) = InterpolateNestedStruct((SBNT_Data.(SubNum).(Visit).(Phase)), 100);
