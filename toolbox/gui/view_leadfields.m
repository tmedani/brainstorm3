function hFig = view_leadfields(HeadmodelFiles)
% VIEW_LEADFIELDS: Show all the leadfield vectors from a "Gain matrix" of the forward model.
% 
% USAGE:  hFig = view_leadfields(HeadmodelFiles)
%
% DOCUMENTATION:
%    - https://www.researchgate.net/publication/260603026_Biomagnetism
%    - http://www.bem.fi/book/11/11x/1119x.htm   

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2020 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: John Mosher, Takfarinas Medani, Francois Tadel, 2020


%% ===== PARSE INPUTS =====
if ischar(HeadmodelFiles)
    HeadmodelFiles = {HeadmodelFiles};
end

%% ===== GET DATA =====
LF_finale = [];
HeadmodelMat = cell(1, length(HeadmodelFiles));
SubjectName = cell(1, length(HeadmodelFiles));
CortexSurface = [];
ChannelMat = [];
Channels = [];
iChannels = [];
allModalities = {'EEG', 'MEG', 'SEEG', 'ECOG'};
markersLocs = [];
bst_progress('start', 'View leadfields', 'Loading headmodels...');
for iFile = 1:length(HeadmodelFiles)
    % Get study description
    [sStudy, iStudy, iHeadModel] = bst_get('HeadModelFile', HeadmodelFiles{iFile});
    % Load lead field matrix
    HeadmodelMat{iFile} = in_bst_headmodel(HeadmodelFiles{iFile});
    % Check the dimensions of the leadfield
    if (iFile >= 2)
        if (size(HeadmodelMat{iFile}.Gain,1) ~= size(HeadmodelMat{1}.Gain,1))
            error('The files have different numbers of sensors.');
        elseif (size(HeadmodelMat{iFile}.Gain,2) ~= size(HeadmodelMat{1}.Gain,2))
            warning(['The overlay of models with different source spaces is not recommended.' 10 ...
                'Only the source space of the first file will be displayed.']);
        end
    end

    % Get the modalities used in this study
    if isempty(sStudy.HeadModel(iHeadModel).EEGMethod)
        allModalities = setdiff(allModalities, 'EEG');
    end
    if isempty(sStudy.HeadModel(iHeadModel).MEGMethod)
        allModalities = setdiff(allModalities, 'MEG');
    end
    if isempty(sStudy.HeadModel(iHeadModel).SEEGMethod)
        allModalities = setdiff(allModalities, 'SEEG');
    end
    if isempty(sStudy.HeadModel(iHeadModel).ECOGMethod)
        allModalities = setdiff(allModalities, 'ECOG');
    end

    % Load channel file
    if (iFile == 1)
        ChannelMat = in_bst_channel(sStudy.Channel.FileName, 'Channel');
    else
        newChanMat = in_bst_channel(sStudy.Channel.FileName, 'Channel');
        if ~isequal({ChannelMat.Channel.Name}, {newChanMat.Channel.Name})
            error('The files have different lists of channels.');
        end
    end
    SubjectName{iFile} = bst_fileparts(bst_fileparts(sStudy.Channel.FileName));

    % Get surface file to display in the figure
    if isempty(CortexSurface) && ~isempty(HeadmodelMat{iFile}.SurfaceFile)
        CortexSurface = HeadmodelMat{iFile}.SurfaceFile;
    end
    % Load the source space
    if isempty(HeadmodelMat{iFile}.GridLoc)
        TessMat = in_tess_bst(HeadmodelMat{iFile}.SurfaceFile);
        HeadmodelMat{iFile}.GridLoc = TessMat.Vertices;
    end
end
if isempty(allModalities)
    error('No modality available for all the files.');
end


%% ===== SELECT MODALITY/REFERENCE =====
% Ask modality + EEG reference
selectedModality = [];
isAvgRef = 1; 
iRef = [];
if ~SelectModality
    bst_progress('stop');
    return;
end
% Update current lead fields
GetLeadField();


%% ===== DISPLAY =====
% Display cortex surface
hFig = [];
SurfAlpha = 0.5;
SurfColor = [0.5 0.5 0.5] ;
hFig = view_surface(CortexSurface, SurfAlpha, SurfColor, hFig);
if isempty(hFig)
    error('No reference surface available');
end
hold on
% Set orientation: left
figure_3d('SetStandardView', hFig, 'left');
% Update figure name
set(hFig, 'Name', ['Leadfield: ' HeadmodelFiles{1}]);
% Get axes handles
hAxes = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
ColorOrder = get(hAxes, 'ColorOrder');
% Initialize list of quiver objects
hQuiver = zeros(length(HeadmodelFiles),1);
% Hack keyboard callback
KeyPressFcn_bak = get(hFig, 'KeyPressFcn');
set(hFig, 'KeyPressFcn', @KeyPress_Callback);

% Create legend
hLabel = uicontrol('Style', 'text', ...
    'String',              '...', ...
    'Units',               'Pixels', ...
    'Position',            [6 20 400 18], ...
    'HorizontalAlignment', 'left', ...
    'FontUnits',           'points', ...
    'FontSize',            bst_get('FigFont'), ...
    'ForegroundColor',     [.3 1 .3], ...
    'BackgroundColor',     [0 0 0], ...
    'Parent',              hFig);

%% ===== DISPLAY LEADFIELD =====
% Current sensor
useLogScale = false;
useLogScaleLegendMsg = 'Off';
iChannel = 1;
iQuiverSize = 1;
iQuiverWidth = 1;
iThresholdLF = 1;
iBalance = 0; % sense of the threshold
DrawArrows();
bst_progress('stop');



%% =================================================================================
%  === INTERNAL CALLBACKS ==========================================================
%  =================================================================================

%% ===== KEYBOARD CALLBACK =====
    function KeyPress_Callback(hFig, keyEvent)
        switch (keyEvent.Key)
            % === LEFT, RIGHT, PAGEUP, PAGEDOWN : Processed by TimeWindow  ===
            case {'leftarrow',}
                if ismember('shift', keyEvent.Modifier)
                    iQuiverSize = iQuiverSize /1.2;
                elseif ismember('control', keyEvent.Modifier)
                    iQuiverWidth = iQuiverWidth / 1.2;
                elseif ismember('alt', keyEvent.Modifier)
                    iThresholdLF = iThresholdLF - 0.01;
                else
                    iChannel = iChannel - 1; 
                end
            case {'rightarrow'}
                 if ismember('shift', keyEvent.Modifier)
                    iQuiverSize = iQuiverSize * 1.2;
                 elseif ismember('control', keyEvent.Modifier)
                     iQuiverWidth = iQuiverWidth * 1.2;
                elseif ismember('alt', keyEvent.Modifier)
                    iThresholdLF = iThresholdLF + 0.01;
                 else
                     iChannel = iChannel + 1;
                 end
            case 'uparrow'
                if ismember('shift', keyEvent.Modifier)
                    iQuiverSize = iQuiverSize * 1.2;
                elseif ismember('control', keyEvent.Modifier)
                    iQuiverWidth = iQuiverWidth * 1.2;
                elseif ismember('alt', keyEvent.Modifier)
                    iThresholdLF = iThresholdLF + 0.01;
                 else                   
                    if ~isempty(iRef)
                        iRef = iRef + 1;
                    end
                end
            case 'downarrow'
                if ismember('shift', keyEvent.Modifier)
                    iQuiverSize = iQuiverSize / 1.2;
                 elseif ismember('control', keyEvent.Modifier)
                    iQuiverWidth = iQuiverWidth / 1.2;
                elseif ismember('alt', keyEvent.Modifier)
                    iThresholdLF = iThresholdLF - 0.01;
                 else
                    if ~isempty(iRef)
                        iRef = iRef - 1;
                    end
                end
            case 'r' %% not for MEG
                SelectReference();
            case 't' %% not for MEG
                SelectTarget();
            case 's'
                if ~isempty(findobj(hFig, 'Tag', 'SetVertices'))
                    delete(findobj(hFig, 'Tag', 'SetVertices'))
                else
                    hold on;
                    plot3(HeadmodelMat{1}.GridLoc(:,1), HeadmodelMat{1}.GridLoc(:,2), HeadmodelMat{1}.GridLoc(:,3), 'r.', ...
                        'Parent', hAxes, ...
                        'Tag', 'SetVertices');
                end
            case 'e'
                hold on;
                % Plot sensors
                if ~isempty(findobj(hAxes, 'Tag', 'allChannel'))
                    delete(findobj(hAxes, 'Tag', 'allChannel'))
                else
                    if length(Channels) > 10
                        hSensors = figure_3d('PlotSensorsNet', hAxes, markersLocs, 0, 0);
                        set(hSensors, 'LineWidth', 1, 'MarkerSize', 5,'Tag','allChannel');
                    end
                end
            case 'm'
                if ~SelectModality()
                    return;
                end
                GetLeadField();
                if ~isempty(findobj(hAxes, 'Tag', 'allChannel'))
                    delete(findobj(hAxes, 'Tag', 'allChannel'))
                end
           
            case 'l'
                if ismember('shift', keyEvent.Modifier)
                    useLogScale = ~useLogScale;
                    if (useLogScale)
                        useLogScaleLegendMsg = 'On';
                    else
                        useLogScaleLegendMsg = 'Off';
                    end
                end
                
             case 'return'
                if ismember('shift', keyEvent.Modifier)
                    useLogScale = ~useLogScale;
                    if (useLogScale)
                        useLogScaleLegendMsg = 'On';
                    else
                        useLogScaleLegendMsg = 'Off';
                    end
                elseif ismember('alt', keyEvent.Modifier)
                     iBalance = ~iBalance;
                else
                    return;    
                end
            case 'h'
                java_dialog('msgbox', ['<HTML><TABLE>' ...
                    '<TR><TD><B>Left arrow</B></TD><TD>Previous target channel (red color)</TD></TR>' ....
                    '<TR><TD><B>Right arrow</B></TD><TD>Next target channel (red color)</TD></TR>'....
                    '<TR><TD><B>Up arrow</B></TD><TD>Previous ref channel (green color)</TD></TR>'....
                    '<TR><TD><B>Down arrow</B></TD><TD>Next ref channel (green color)</TD></TR>'....
                    '<TR><TD><B>Shift + uparrow</B></TD><TD>Increase the vector length</TD></TR>'...
                    '<TR><TD><B>Shift + downarrow</B></TD><TD>Decrease the vector length</TD></TR>'...                   
                    '<TR><TD><B>Control + uparrow</B></TD><TD>Increase the vector width</TD></TR>'...
                    '<TR><TD><B>Control + downarrow</B></TD><TD>Decrease the vector width</TD></TR>'... 
                    '<TR><TD><B>Shift + L</B></TD><TD>Toggle on/off logarithmic scale</TD></TR>'...
                    '<TR><TD><B>M</B></TD><TD>Change the <B>M</B>odality (MEG, EEG, SEEG, ECOG)</TD></TR>'....
                    '<TR><TD><B>R</B></TD><TD>Select the <B>R</B>eference channel</TD></TR>'....
                    '<TR><TD><B>T</B></TD><TD>Select the <B>T</B>arget channel</TD></TR>'....
                    '<TR><TD><B>S</B></TD><TD>Show/hide the source grid</TR>'....
                    '<TR><TD><B>E</B></TD><TD>Show/hide the sensors</TD></TR>'...
                    '<TR><TD><B>0 to 9</B></TD><TD>Change view</TD></TR>'...
                    '</TABLE>'], 'Keyboard shortcuts');
           
                  
                % === NUMBERS : VIEW SHORTCUTS ===
            case '0'
                SetStandardView(hFig, {'left', 'right', 'top'});
            case '1'
                SetStandardView(hFig, 'left');
                
            case '2'
                SetStandardView(hFig, 'bottom');
            case '3'
                SetStandardView(hFig, 'right');
                
            case '4'
                SetStandardView(hFig, 'front');
                
            case '5'
                SetStandardView(hFig, 'top');
                
            case '6'
                SetStandardView(hFig, 'back');
                
            case '7'
                SetStandardView(hFig, {'left', 'right'});
                
            case '8'
                SetStandardView(hFig, {'bottom', 'top'});
                
            case '9'
                SetStandardView(hFig, {'front', 'back'});
                
            otherwise
                KeyPressFcn_bak(hQuiver, keyEvent);
                return;
        end
        % Redraw arrows
        if (iChannel <= 0)
            iChannel = length(Channels);
        end
        if (iChannel > length(Channels))
            iChannel = 1;
        end
        % Redraw arrows
        if (iRef <= 0)
            iRef = length(Channels);
        end
        if (iRef > length(Channels))
            iRef = 1;
        end
        DrawArrows();
    end

%% ===== DRAW CURRENT CHANNEL =====
    function DrawArrows()
        % Delete previous Channels and sensors
        delete(findobj(hAxes, '-depth', 1, 'Tag', 'lfArrows'));
        if isprop(hAxes, 'ColorOrderIndex')
            set(hAxes, 'ColorOrderIndex', 1);
        else
            set(hAxes, 'NextPlot', 'new');
        end
        % Draw current LF
        for iLF = 1:length(LF_finale)
            % EEG
            if ismember(selectedModality, {'EEG','ECOG','SEEG'})
                if isAvgRef
                    LeadField = LF_finale{iLF}(iChannel,:) - mean(LF_finale{iLF},1);
                    strRef = ' / Ref: Avg';
                else
                    LeadField = LF_finale{iLF}(iChannel,:) - LF_finale{iLF}(iRef,:);
                    strRef = Channels(iRef).Name;
                end
            else % MEG
                LeadField = LF_finale{iLF}(iChannel,:);
                strRef = '';
                if ~isempty(findobj(hAxes, '-depth', 1, 'Tag', 'RefChannel'))
                    delete(findobj(hAxes, '-depth', 1, 'Tag', 'RefChannel'));
                end
            end
            % Display arrows
            LeadField = reshape(LeadField,3,[])'; % each column is a vector
            if(useLogScale)
                LeadField = logScaledLeadField(LeadField);
            end
            % thresholding
            normLF = sqrt((LeadField(:,1) ).^2 +(LeadField(:,2) ).^2 + (LeadField(:,3)).^2);
            [col1, ind] = sort(normLF, 'ascend');
            LeadFieldRordered = LeadField(ind,:);
            cdf = cumsum(col1); % Compute cdf
            cdf = cdf/cdf(end); % Normalize
            % Find index bellow or above the thresholding
            if iBalance == 0 % 0 ==> inferior and 1 is superior 
                index = find(cdf <= iThresholdLF);
                iSymbole = '<=';
            else
                index = find(cdf > iThresholdLF);   
                iSymbole = '>';
            end
            dataValue = zeros(size(LeadFieldRordered));
            dataValue(index,:) = LeadFieldRordered(index,:);

            hQuiver(iLF) = quiver3(...
                ...HeadmodelMat{iLF}.GridLoc(:,1), HeadmodelMat{iLF}.GridLoc(:,2), HeadmodelMat{iLF}.GridLoc(:,3), ...
                ...LeadField(:,1), LeadField(:,2), LeadField(:,3), ...
                HeadmodelMat{iLF}.GridLoc(ind,1), HeadmodelMat{iLF}.GridLoc(ind,2), HeadmodelMat{iLF}.GridLoc(ind,3), ...
                dataValue(:,1), dataValue(:,2), dataValue(:,3), ...
                iQuiverSize, ...
                'Parent',    hAxes, ...
                'LineWidth', iQuiverWidth, ...
                'Color',     ColorOrder(mod(iLF-1, length(ColorOrder)) + 1, :), ...
                'Tag',       'lfArrows');
            % Arrow legends
            strLegend{iLF} = [SubjectName{iLF} ' : ' selectedModality  ' ' HeadmodelMat{iLF}.Comment];
            hold on
        end

        % Remove previous selected sensor
        delete(findobj(hAxes, '-depth', 1, 'Tag', 'SelChannel'));
        % Plot selected sensor
        if ~isempty(Channels(iChannel).Loc) && ~ismember(Channels(iChannel).Name, {'EEG','MEG','MEG MAG', 'MEG GRAD'})
            line(Channels(iChannel).Loc(1,1), Channels(iChannel).Loc(2,1), Channels(iChannel).Loc(3,1), ...
                'Parent',          hAxes, ...
                'LineWidth',       2, ...
                'LineStyle',       'none', ...
                'Marker',          'o', ...
                'MarkerFaceColor', [1 0 0], ...
                'MarkerEdgeColor', [.4 .4 .4], ...
                'MarkerSize',      8, ...
                'Tag',             'SelChannel');
        end
        if ~strcmp(selectedModality,'MEG')
            % Remove previous selected reference
            delete(findobj(hAxes, '-depth', 1, 'Tag', 'RefChannel'));
            % Plot the reference electrode
            if ~isAvgRef && ~isempty(Channels(iRef).Loc)
                line(Channels(iRef).Loc(1,1), Channels(iRef).Loc(2,1), Channels(iRef).Loc(3,1), ...
                    'Parent',          hAxes, ...
                    'LineWidth',       2, ...
                    'LineStyle',       'none', ...
                    'Marker',          '*', ...
                    'MarkerFaceColor', [0 1 1], ...
                    'MarkerEdgeColor', [.4 .8 .4], ...
                    'MarkerSize',      8, ...
                    'Tag',             'RefChannel');
            end
            % Title bar (channel name)
            if isAvgRef
                strTitle = sprintf('Target channel(red) #%d/%d  (%s) | %s Ref Channel(green) = AvgRef  | iThresholdLF %s %s %%| Log. scale %s', iChannel, length(Channels), Channels(iChannel).Name,selectedModality, iSymbole,  num2str(iThresholdLF*100),useLogScaleLegendMsg);
            else
                strTitle = sprintf('Target channel(red) #%d/%d  (%s) | %s Ref Channel(green) = %s| iThresholdLF %s %s %%| Log. scale %s', iChannel, length(Channels), Channels(iChannel).Name,selectedModality,Channels(iRef).Name, iSymbole, num2str(iThresholdLF*100),useLogScaleLegendMsg);
            end
        else
            strTitle = sprintf('Target channel (red) #%d/%d  (%s) | iThresholdLF %s %s %%|Log. scale %s', iChannel, length(Channels), Channels(iChannel).Name, iSymbole,num2str(iThresholdLF*100),useLogScaleLegendMsg);
        end
        
        if (iChannel == 1) && (length(Channels) > 1)
            strTitle = [strTitle, '       [Press arrows for next/previous channel (or H for help)]'];
        end
        set(hLabel, 'String', strTitle, 'Position', [10 1 1500 35],'ForegroundColor', [1 1 1]);
        % Arrows legend
        legend(hQuiver, strLegend, ...
            'TextColor',   'w', ...
            'fontsize',    bst_get('FigFont'), ...
            'Interpreter', 'None', ...
            'Location',    'NorthEast', ...
            'Tag',         'LegendOverlay');
        legend('boxoff');
    end


%% ===== SET MODALITY =====
    function isOk = SelectModality()
        % Ask modality if there is more than one
        if (length(allModalities) > 1)
            selectedModality = java_dialog('question', 'Select the modality:', ...
                'Display the Lead Field', [], allModalities, allModalities{1});
            if isempty(selectedModality)
                isOk = 0;
                return;
            end
        elseif (length(allModalities) == 1)
            selectedModality = allModalities{1};
        end
        % Get the corresponding channels
        iChannels = good_channel(ChannelMat.Channel, [], selectedModality);
        Channels = ChannelMat.Channel(iChannels);
        % Get channels locations
        if length(Channels) > 10
            markersLocs = cell2mat(cellfun(@(c)c(:,1), {Channels.Loc}, 'UniformOutput', 0))';
        end
        % Ask to select reference
        isOk = SelectReference();
    end


%% ===== SET REFERENCE =====
    function isOk = SelectReference()
        isOk = 1;
        if ~strcmp(selectedModality,'MEG')
%             [isAvgRef, isCancel] = java_dialog('confirm', ...
%                 ['<HTML>Do you want to use the <B>average refence</B> for the ' selectedModality ' ?<BR>'...
%                 'Otherwise you will choose one reference electrode.'], [selectedModality ' average reference'], [], ...
%                 {'Yes, use average reference'}, 1);
%             if isCancel
%                 isOk = 0;
%                 return;
%             end
%             if ~isAvgRef
%                 % Ask for the reference electrode
%                 refChan = java_dialog('combo', '<HTML>Select the reference channel:<BR><BR>', [selectedModality ' reference'], [], {Channels.Name});
%                 if isempty(refChan)
%                     isOk = 0;
%                     return;
%                 end
%                 iRef = find(strcmpi({Channels.Name}, refChan));
%             end
                % Ask for the reference electrode
                refChan = java_dialog('combo', '<HTML>Select the reference channel (green color):<BR><BR>', [selectedModality ' reference'], [], {'Average Ref', Channels.Name});
                if isempty(refChan)
                    isOk = 0;
                    return;
                end                
                iRef = find(strcmpi({Channels.Name}, refChan));
                if isempty(iRef)
                    isAvgRef = 1;
                else
                    isAvgRef = 0;
                end                   
        end
    end

%% ===== SET TARGET =====
    function isOk = SelectTarget()
        isOk = 1;
        % Ask for the target electrode
        trgChan = java_dialog('combo', '<HTML>Select the target channel (red color):<BR><BR>', [selectedModality ' Target'], [], {Channels.Name});
        if isempty(trgChan)
            isOk = 0;
            return;
        end
        iChannel = find(strcmpi({Channels.Name}, trgChan));
    end

    %% ===== GET LEADFIELD =====
    function GetLeadField       
        % Update the LF according to the selected channels only
        for iLF = 1:length(HeadmodelFiles)
            LF_finale{iLF} = HeadmodelMat{iLF}.Gain(iChannels,:);
        end
    end

%% ===== LEADFIELD TO LOG SPACE =====
    function lf_log = logScaledLeadField(lf)
        lf_2 = lf.^2;
        r = sqrt(sum(lf_2,2));
        rho = sqrt(lf_2(:,1) + lf_2(:,2));
        t = atan2(rho,lf(:,3));
        f = atan2(lf(:,2),lf(:,1));
        lf_log = [ log10(r) .* sin(t) .* cos(f) ...
                   log10(r) .* sin(t) .* sin(f) ...
                   log10(r) .* cos(t)];
    end

    %= = From this line, these functions are part of the figure3d
    %% ===== SET STANDARD VIEW =====
    function SetStandardView(hFig, viewNames)
        % Make sure that viewNames is a cell array
        if ischar(viewNames)
            viewNames = {viewNames};
        end
        % Get Axes handle
        hAxes = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
        % Get the data types displayed in this figure
        ColormapInfo = getappdata(hFig, 'Colormap');
        % Get surface information
        TessInfo = getappdata(hFig, 'Surface');

        % ===== ANATOMY ORIENTATION =====
        % If MRI displayed in the figure, use the orientation of the slices, instead of the orientation of the axes
        R = eye(3);
        % Get the mri surface
        Ranat = [];
        if ismember('anatomy', ColormapInfo.AllTypes)
            iTess = find(strcmpi({TessInfo.Name}, 'Anatomy'));
            if ~isempty(iTess)
                % Get the subject MRI structure in memory
                sMri = bst_memory('GetMri', TessInfo(iTess).SurfaceFile);
                % Calculate transformation: SCS => MRI  (inverse MRI => SCS)
                Ranat = pinv(sMri.SCS.R);
            end
        end
        % Displaying a surface: Load the SCS field from the MRI
        if isempty(Ranat) && ~isempty(TessInfo) && ~isempty(TessInfo(1).SurfaceFile)
            % Get subject
            sSubject = bst_get('SurfaceFile', TessInfo(1).SurfaceFile);
            % If there is an MRI associated with it
            if ~isempty(sSubject) && ~isempty(sSubject.Anatomy) && ~isempty(sSubject.Anatomy(sSubject.iAnatomy).FileName)
                % Load the SCS+MNI transformation from this file
                sMri = load(file_fullpath(sSubject.Anatomy(sSubject.iAnatomy).FileName), 'NCS', 'SCS', 'Comment');
                if isfield(sMri, 'NCS') && isfield(sMri.NCS, 'R') && ~isempty(sMri.NCS.R) && isfield(sMri, 'SCS') && isfield(sMri.SCS, 'R') && ~isempty(sMri.SCS.R)
                    % Calculate the SCS => MNI rotation   (inverse(MRI=>SCS) * MRI=>MNI)
                    Ranat = sMri.NCS.R * pinv(sMri.SCS.R);
                end
            end
        end
        % Get the rotation to change orientation
        if ~isempty(Ranat)
            R = [0 1 0;-1 0 0; 0 0 1] * Ranat;
        end    

        % ===== MOVE CAMERA =====
        % Apply the first orientation to the target figure
        switch lower(viewNames{1})
            case {'left', 'right_intern'}
                newView = [0,1,0];
                newCamup = [0 0 1];
            case {'right', 'left_intern'}
                newView = [0,-1,0];
                newCamup = [0 0 1];
            case 'back'
                newView = [-1,0,0];
                newCamup = [0 0 1];
            case 'front'
                newView = [1,0,0];
                newCamup = [0 0 1];
            case 'bottom'
                newView = [0,0,-1];
                newCamup = [1 0 0];
            case 'top'
                newView = [0,0,1];
                newCamup = [1 0 0];
        end
        % Update camera position
        view(hAxes, newView * R);
        camup(hAxes, double(newCamup * R));
        % Update head light position
        camlight(findobj(hAxes, '-depth', 1, 'Tag', 'FrontLight'), 'headlight');
        % Select only one hemisphere
        if any(ismember(viewNames, {'right_intern', 'left_intern'}))
            bst_figures('SetCurrentFigure', hFig, '3D');
            drawnow;
            if strcmpi(viewNames{1}, 'right_intern')
                panel_surface('SelectHemispheres', 'right');
            elseif strcmpi(viewNames{1}, 'left_intern')
                panel_surface('SelectHemispheres', 'left');
            else
                panel_surface('SelectHemispheres', 'none');
            end
        end

        % ===== OTHER FIGURES =====
        % If there are other view to represent
        if (length(viewNames) > 1)
            hClones = bst_figures('GetClones', hFig);
            % Process the other required views
            for i = 2:length(viewNames)
                if ~isempty(hClones)
                    % Use an already cloned figure
                    hNewFig = hClones(1);
                    hClones(1) = [];
                else
                    % Clone figure
                    hNewFig = bst_figures('CloneFigure', hFig);
                end
                % Set orientation
                SetStandardView(hNewFig, viewNames(i));
            end
            % If there are some cloned figures left : close them
            if ~isempty(hClones)
                close(hClones);
                % Update figures layout
                gui_layout('Update');
            end
        end
    end
end
