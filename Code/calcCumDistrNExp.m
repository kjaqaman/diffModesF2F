function cumDistrNExp = calcCumDistrNExp(param,abscissa)
%CALCCUMDISTRNEXP calculates the cumulative distribution of N exponentials
%
%SYNOPSIS cumDistrNExp = calcCumDistrNExp(param,abscissa)
%
%INPUT  param         : Vector of parameters indicating the means and
%                       amplitude of the exponentials.
%       abscissa      : Abscissa values at which the cumulative
%                       distribution is calculated.
%
%OUTPUT cumDistrNExp  : Values of the resulting cumulative distribution
%                       given the input abscissa values.
%       jacobianMat   : Jacobian matrix.

%Khuloud Jaqaman, June 2012
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

%% Output

cumDistrNExp = [];

%% Input

%check whether correct number of input arguments was used
if nargin < 2
    disp('--calcCumDistrNExp: Incorrect number of input arguments!');
    return
end

%% Calculating the cumulative distribution

%get number of exponentials
numExp = length(param)/2;

%get their means and amplitudes
expMean = param(1:numExp);
expAmp  = param(numExp+1:end);

%get number of data points
numData = length(abscissa);

%calculate the cumulative distribution and the jacobian matrix
cumDistrNExp = zeros(numData,1);
% jacobianMat = zeros(numData,2*numExp);
for iExp = 1 : numExp
    cumDistrNExp = cumDistrNExp + expAmp(iExp)*(1 - exp(-abscissa/expMean(iExp)));
    %     jacobianMat(:,iExp) = -(expAmp(iExp)/expMean(iExp)^2)*exp(-abscissa/expMean(iExp));
    %     jacobianMat(:,numExp+iExp) = 1 - exp(-abscissa/expMean(iExp));
end
    
%% ~~ the end ~~
    
