classdef DiffusionModeAnalysisProcess < PostTrackingProcess
    % A concrete class for analyzing track diffusion modes (instead of
    % traditional diffusion types)
    
    % Khuloud Jaqaman, August 2017
 %
% Copyright (C) 2020, Jaqaman Lab - UTSouthwestern
%
% This file is part of diffModesF2F.
%
% diffModesF2F is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% diffModesF2F is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with diffModesF2F.  If not, see <http://www.gnu.org/licenses/>.
%
%
   
    methods (Access = public)
        function obj = DiffusionModeAnalysisProcess(owner, varargin)
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                super_args{1} = owner;
                super_args{2} = DiffusionModeAnalysisProcess.getName;
                super_args{3} = @analyzeMovieDiffusionModes;
                if isempty(funParams)  % Default funParams
                    funParams = DiffusionModeAnalysisProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@PostTrackingProcess(super_args{:});
        end

        function h=draw(obj,iChan,varargin)
            h = obj.draw@PostTrackingProcess(iChan,varargin{:},'useCache',true);
        end

        function varargout = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            outputList = {'diffModeAnalysisRes', 'diffModeDecomposition', 'tracks'};
            ip =inputParser;
            ip.addRequired('iChan',@(x) isscalar(x) && obj.checkChanNum(x));
            ip.addOptional('iFrame',[],@(x) isempty(x) || isscalar(x) && obj.checkFrameNum(x));
            ip.addParamValue('iZ',[], @(x)ismember(x,1:obj.owner_.zSize_));
            ip.addParamValue('useCache',false,@islogical);
            ip.addParamValue('output',outputList,@(x) all(ismember(x,outputList)));
            ip.parse(iChan,varargin{:})
            iFrame = ip.Results.iFrame;
            output = ip.Results.output;
            if ischar(output),output={output}; end
            nOutput = numel(output);
            
            % Data loading

            % load outFilePaths_{1,iChan}
            s = cached.load(obj.outFilePaths_{1,iChan}, '-useCache', ip.Results.useCache, output{:});
            
            varargout = cell(nOutput);
            for i = 1:nOutput
                switch output{i}
                    case 'tracks'
                        tracksFinal = s.(output{i});
                        if ~isempty(iFrame)
                            % Filter tracks existing in input frame
                            trackSEL=getTrackSEL(tracksFinal);
                            validTracks = (iFrame>=trackSEL(:,1) &iFrame<=trackSEL(:,2));
                            [tracksFinal(~validTracks).tracksCoordAmpCG]=deal([]);
                            
                            for j=find(validTracks)'
                                tracksFinal(j).tracksCoordAmpCG = tracksFinal(j).tracksCoordAmpCG(:,1:8*(iFrame-trackSEL(j,1)+1));
                            end
                            varargout{i} = tracksFinal;
                        else
                            varargout{i} = tracksFinal;
                        end
                    case 'diffModeAnalysisRes'
                        varargout{i} = s.(output{i});
                    case 'diffModeDecomposition'
                        varargout{i} = s.(output{i});
                end
            end
        end
        
        function output = getDrawableOutput(obj)
            types = DiffusionModeAnalysisProcess.getTrackTypes();
            colors = vertcat(types.color);
            output(1).name='Classified tracks';
            output(1).var='tracks';
            output(1).formatData=@DiffusionModeAnalysisProcess.formatTracks;
            output(1).type='overlay';
            output(1).defaultDisplayMethod=@(x)TracksDisplay('Color', colors);
        end
        
    end
    
    methods (Static)
        
        function name = getName()
            name = 'Diffusion Mode Analysis';
        end
        
        %         function h = GUI()
        %             h = @motionAnalysisProcessGUI;
        %         end
        %
        %         function alpha = getAlphaValues()
        %             alpha=[0.01 0.05 0.1 0.2];
        %         end
        %
        %         function methods = getConfinementRadiusMethods()
        %             methods(1).type = 0;
        %             methods(1).name = 'Mean positional standard deviation';
        %             methods(2).type = 1;
        %             methods(2).name = 'Minimum positional standard deviation';
        %             methods(3).type = 2;
        %             methods(3).name = 'Rectangle approximation';
        %         end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1 : numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'DiffusionModeAnalysis'];

            if owner.is3D
                funParams.probDim = 3;
            else
                funParams.probDim = 2;
            end
            
            funParams.minLength = 5;
            funParams.alpha = 0.01;
            funParams.showPlot = 1;
            funParams.maxNumMode = 10;
            funParams.binStrategy = 2;
            funParams.plotName = 'Figure';
            funParams.subSampSize = [];
            funParams.doControl = 1;
            funParams.diffModeDividerStruct = [];
            funParams.forceDecompose = 1;
            
        end
        
        function displayTracks = formatTracks(tracks)
            % Format classified tracks into structure for display
            
            % Read track types and classification matrix
            types = DiffusionModeAnalysisProcess.getTrackTypes();
            % Column vector of the number of track segments
            track_class = vertcat(tracks.classification);

            % Number of labels per track needed
            nLabels = cellfun('size',{tracks.classification},1);
            % Labels is a cell array of indices corresponding to track segment types,
            %   see getTrackTypes
            % Initialize all labels as unlabeled
            labels = arrayfun(@(x) ones(x,1)*numel(types),nLabels,'UniformOutput',false);
            % Map indicates position of last label for each track
            map = cumsum(nLabels);
            
            % logical array per track index of if there is more than one label per track
            nLabels_gt_1 = nLabels > 1;
            % labels index where the labels for each track starts
            %  if there is more than one label per track
            start = map(nLabels_gt_1)-nLabels(nLabels_gt_1)+1;
            % labels index where the labels for each track ends
            %  if there is more than one label per track
            finish = map(nLabels_gt_1);

            % Set labels as needed
            for i = 1 : numel(types) - 1
                % idx is a logical array if each track _segment_ belongs to types(i)
                idx = types(i).f(track_class);
                % idx2 is a cell array of logical arrays
                %  the number of cells corresponds to each track index
                %  the index of the logical array in each cell refers to each track segment
                % Setup logical arrays. This works for when nLabels == 1
                idx2 = num2cell(idx(map));
                % deal with nLabels > 1 separately, grab range of labels for each track segment
                %  corresponding to each track        
                idx2(nLabels_gt_1) = arrayfun(@(s,e) idx(s:e),start,finish,'UniformOutput',false);
                % idx is now a cell array of logical arrays marking for each track
                %  if each track segment belongs to types(i)
                idx = idx2;
                % Select only the indices where at least one segment belongs to type(i)
                any_idx = cellfun(@any,idx);
                % Assign label as i for each track segment belonging to types(i)
                % Merge with previous labels
                labels(any_idx) = cellfun(@(idx_i,labels_i) labels_i.*~idx_i + i.*idx_i, ...
                    idx(any_idx), labels(any_idx)','UniformOutput',false);
            end
            % Assign labels to each track
            [tracks.label] = deal(labels{:});

            % Format tracks using TrackingProcess utility function
            displayTracks = TrackingProcess.formatTracks(tracks);
        end
        
        function types = getTrackTypes()
            % Get the color map for classified tracks
            %
            % Mode 1: brownish red
            types(1).name = 'Mode 1';
            types(1).f = @(x) x(:, 1) == 1;
            types(1).color = [0.8 0.2 0];
            % Mode 2: blue
            types(2).name = 'Mode 2';
            types(2).f = @(x) x(:, 1) == 2;
            types(2).color = [0 0 1];
            % Mode 3: cyan
            types(3).name = 'Mode 3';
            types(3).f = @(x) x(:, 1) == 3;
            types(3).color = [0 1 1];
            % Mode 4: magenta
            types(4).name = 'Mode 4';
            types(4).f = @(x) x(:, 1) == 4;
            types(4).color = [1 0 1];
            % Mode 5: yellow
            types(5).name = 'Mode 5';
            types(5).f = @(x) x(:, 1) == 5;
            types(5).color = [1 0.7 0];
            % Unclassified: gray
            types(6).name = 'Too short for classification';
            types(6).f = @(x) isnan(x(:,1));
            types(6).color = [0.7 0.7 0.7];
        end
        
        function hIm = showTrackTypes(newFigure)
            % DiffusionModeAnalysisProcess.showTypes: Show types and their colors
            
            if(nargin < 1)
                newFigure = false;
            end
            if(newFigure)
                figure;
            end
            
            types = DiffusionModeAnalysisProcess.getTrackTypes;
            
            % Create a grey background
            hIm = imshow(ones(300)*0.3);
            % Display type names in their color from top to bottom
            for i=1:length(types)
                text(1,length(types)-i+1,types(i).name, ...
                    'Color',types(i).color, ...
                    'Units','characters');
            end
        end
        
    end
end
