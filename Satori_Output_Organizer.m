%% fNIRS Data Organizer for Satori Output
% Originating Author: Andrew Hagen
% Last Revised: 5/9/2024

% This script imports the raw Satori GLM output xlsx and organizes the Oxy Results sheet to have only beta weights, with conditions and participants labeled across each row, and channels across each column

% Prompts 
    % 1. Select raw (unedited) Satori GLM output file
    % 2. How many participants in analysis - crucial for column deletion and participant labels
    % 3. How many conditions (event names)
    % 4. List each condition name in separate prompts
    % 5. Y or N for exporting
    % 6. Y or N for if you want the organized data to be exported as a separate file in the same folder as the original file or if you want to add the organized data as a new sheet in the original file

clc
clear
close all

% Prompt the user to select an Excel file
[filename, filepath] = uigetfile('*.xlsx', 'Select an Excel file');
fullFilepath = fullfile(filepath, filename);

% Prompt the user for the number of participants
numParticipants = inputdlg('Enter the number of participants: ','Participants');
numParticipants = str2double(numParticipants{1}); % Convert the string to a number

% Prompt the user for condition names
numConditions = inputdlg('Enter the number of conditions: ', 'Conditions');
numConditions = str2double(numConditions{1}); % Convert the string to a number
conditionNames = cell(1,numConditions);
for i = 1:numConditions
    conditionNames{i} = inputdlg(['Enter name for condition ' num2str(i) ': '], 'Condition Names');
end

%Import 'MULTI-STUDY GLM RESULTS' sheet as cell array
   participantsArray = readcell(fullFilepath,"FileType","spreadsheet","Sheet",1);
   participantsArray = participantsArray(6:(5 + numParticipants),:);
%Import 'Oxy Results' sheet as cell array
    dataArray = readcell(fullFilepath,"FileType","spreadsheet","Sheet",2); 

%% Delete all columns after column 3 except for every fourth column (This is the beta weight column
numCols = size(dataArray, 2);
columnsToKeep = [1:3, 7:4:numCols]; % Columns 1, 2, 3, and every fourth column after 3
dataArray = dataArray(:, columnsToKeep);

%% Fill in missing condition names for each participant
colIndex = 3; % Starting column index for condition names
for Participant = 1:numParticipants
    for condition = 1:numel(conditionNames)
        dataArray{1, colIndex} = conditionNames{condition};
        colIndex = colIndex + 1;
    end
end

% Delete extra columns containing residuals - There is group of 4 columns for each participant
dataArray(:,colIndex:end) = [];

%% Insert the participant IDs starting from column 3
% Calculate the total number of columns needed for participant IDs
totalColumns = 2 + numParticipants * numConditions;

% Create a new row to insert at the top of dataArray
participantRow = cell(1, totalColumns);

% Insert the participant IDs starting from column 3
for i = 1:numParticipants
    participantID = participantsArray{i,i}; % Get participant ID
    if isnumeric(participantID)
        participantID = num2str(participantID); % Convert numeric ID to string
    end
    startCol = 2 + (i - 1) * numConditions + 1;
    endCol = startCol + numConditions - 1;
    participantRow(1, startCol:endCol) = repmat({participantID}, 1, numConditions);
end

% Insert the participant row at the top of dataArray
dataArray = [participantRow; dataArray];

%% Transpose dataArray and label
dataArray = dataArray';
dataArray{1,1} = 'ParticipantID';
dataArray{1,2} = 'Condition';
dataArray{2,2} = [];
dataArray(:,3) = [];

% Convert all cell contents to strings so we can export
for i = 1:numel(dataArray)
    if iscell(dataArray{i})
        dataArray{i} = cell2mat(dataArray{i});
    end
end

%% Export dataArray based on users preference
exportQ1 = inputdlg('Would you like to export the data? Y or N:', 'Export');
    if strcmp(exportQ1,'Y')
        exportQ2 = questdlg('Would you like to save the data to the orignal xlsx on a new sheet or save the data onto a new xlsx?:', 'Export Location', 'Save to original xlsx', 'Save to new xlsx', 'Cancel', 'Cancel');
            if strcmp(exportQ2,'Save to original xlsx') 
                writecell(dataArray, (fullFilepath), 'Sheet','Oxy_Betas_Organized');
                disp('Data exported successfully :)')
            elseif strcmp(exportQ2,'Save to new xlsx') 
                writecell(dataArray, [filepath 'Oxy_Betas_Organized.xlsx']);
                disp('Data exported successfully :)')
            end
    end
