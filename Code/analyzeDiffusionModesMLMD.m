function analyzeDiffusionModesMLMD(MLMD,minLength,alpha,showPlot,maxNumMode,...
    binStrategy,plotName,subSampSize,doControl,diffModeDividerStruct,forceDecompose,channel2analyze)
%analyzeDiffusionModesMLMD performs diffusion mode analysis on single particle tracks
%
%SYNOPSIS analyzeDiffusionModesMLMD(MLMD,minLength,alpha,showPlot,maxNumMode,...
%    binStrategy,plotName,subSampSize,doControl,diffModeDividerStruct,forceDecompose)
%
%INPUT 
%   Mandatory
%       MLMD: MovieList or MovieData object for movie(s) to be analyzed
%   Optional
%       - forceDecompose: 1 to force run diffusion mode decomposition, 0 to
%       use old decomposition results if available. Default: 1.
%       - channel2analyze: Index of channel to analyze. Default: 1.
%       - For all other parameters, see functions "getDiffModes" and
%       "trackDiffModeAnalysis" for parameter information.
%
%OUPUT Output is saved in directory DiffusionModeAnalysis belonging to each
%      analyzed movie. See functions "getDiffModes" and
%      "trackDiffModeAnalysis" for detailed description of saved output.
%
%Khuloud Jaqaman, August 2017
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
%% Input

%check whether correct number of input arguments was used
if nargin < 1
    disp('--analyzeDiffusionModesMLMD: Function needs at least 1 input argument!');
    return
end

%assign default values of optional input variables
%check rest of input

if nargin < 2 || isempty(minLength)
    minLength = 5;
end

if nargin < 3 || isempty(alpha)
    alpha = 0.01;
end

if nargin < 4 || isempty(showPlot)
    showPlot = 1;
end

if nargin < 5 || isempty(maxNumMode)
    maxNumMode = 10;
end

if nargin < 6 || isempty(binStrategy)
    binStrategy = 2;
end

if nargin < 7 || isempty(plotName)
    plotName = 'Figure';
end

if nargin < 8 || isempty(subSampSize)
    subSampSize = [];
end
if isempty(subSampSize)
    numSamp = 1;
else
    numSamp = 10;
    showPlot = 0;
end

if nargin < 9 || isempty(doControl)
    doControl = 1;
end

if nargin < 10 || isempty(diffModeDividerStruct)
    diffModeDividerStruct = [];
end

if nargin < 11 || isempty(forceDecompose)
    forceDecompose = 1;
end

if nargin < 12 || isempty(channel2analyze)
    channel2analyze = 1;
end

%% Analysis

%determine if input is a MovieList or MovieData object
if isa(MLMD,'MovieList') %if it is a movieDist
    
    listFlag = 1;
    
    %rename to ML
    ML = MLMD;
    clear MLMD
    
    %get number of movies
    numMovies = length(ML.movieDataFile_);
    
else %if it is a movieData
    
    listFlag = 0;
    
    %rename to MD
    MD = MLMD;
    clear MLMD
    
    %assign number of movies to 1
    numMovies = 1;
    
end

%go over all movies and run diffusion mode analysis
for iM = 1 : numMovies
    
    %get movieData of current movie
    if listFlag == 1
        MD = MovieData.load(ML.movieDataFile_{iM});
    end
    
    %add diffusion mode analysis process if never run before
    iProc = MD.getProcessIndex('DiffusionModeAnalysisProcess',1,0);
    if isempty(iProc)
        iProc=numel(MD.processes_)+1;
        MD.addProcess(DiffusionModeAnalysisProcess(MD));
    end
    
    %define function parameters
    
    %general
    funParams = MD.getProcess(iProc).funParams_;
    funParams.ChannelIndex = channel2analyze;
    
    %function-specific
    funParams.minLength = minLength;
    funParams.alpha = alpha;
    funParams.showPlot = showPlot;
    funParams.maxNumMode = maxNumMode;
    funParams.binStrategy = binStrategy;
    funParams.plotName = plotName;
    funParams.subSampSize = subSampSize;
    funParams.doControl = doControl;
    funParams.diffModeDividerStruct = diffModeDividerStruct;
    funParams.forceDecompose = forceDecompose;
    
    %general again
    parseProcessParams(MD.getProcess(iProc),funParams);
    
    %% Run analysis processes
    cellfun(@run,MD.processes_(iProc));
    
end

