function [ DSout ] = mergeChannels( sourceCh , targetCh , varargin )
% Add additional channels to existing dataset in drEEMLC format
%   
% USEAGE:
%           [ DSout ] = mergeChannels( sourceCh , targetCh , varargin )
%
% INPUTS
%            sourceCh: old index of channel that is to be merged
%            targetCh: new index of channel that is to be merged. This might be reorganized
%                      when empty channels are found.
%            varargin: datasets to be merged, first input will be basis of merging (targetCh),
%                      all other datasets will be sources (sourceCh)
%
% OUTPUTS
%               DS: data structure containing HPLC data
%
% NOTES:    In case nSample differs between datasets (eg 20 2D scans are merged with 1 absorbance run),
%           the first sample will be transfered and the rest filled up with NaNs
% Examples:
%       DSmerged=mergeChannels(7,8,DS2dScan,DSpda);
%
% Notice:
% This mfile is part of the LC-addon for drEEM.
% % Check https://github.com/urbanwuensch/LC-addon-for-drEEM for the latest version
%
%
%
% mergeChannels.m: Copyright (C) 2017 Urban J Wünsch
% Technical University of Denmark
% National Institute of Aquatic Resources
% Section for Oceans and Arctic
% Kemitorvet
% 2800 Kgs. Lyngby, Denmark
% urbw@aqua.dtu.dk
%
%
% Version 1, January 2016 First version


if nargin==0
    help mergeChannels
    return
end

nDS= length(varargin);
DSout=varargin{1};

if size(sourceCh)~=size(targetCh)
    error('Size mismatch between targetCh and SourceCh')
end
if nDS-1~=size(sourceCh,1)
    error('Indicies of channels do not match number of source datasets.')
end

for i_DS=1:size(sourceCh,1)
    for i_sCh=1:size(sourceCh,2)
        sourceCh(i_DS,i_sCh);
        targetCh(i_DS,i_sCh);
        DSout.ChData(targetCh(i_DS,i_sCh))=varargin{i_DS+1}.ChData(sourceCh(i_DS,i_sCh));
        if DSout.nSample~=varargin{i_DS+1}.nSample
            if ndims(varargin{i_DS+1}.ChData(sourceCh(i_DS,i_sCh)))==3
                DSout.ChData(targetCh(i_DS-1,i_sCh)).data(varargin{i_DS+1}.nSample+1:DSout.nSample,:,:)=0;
            elseif ndims(varargin{i_DS+1}.ChData(sourceCh(i_DS,i_sCh)))==2
                DSout.ChData(targetCh(i_DS,i_sCh)).data(varargin{i_DS+1}.nSample+1:DSout.nSample,:)=0;
            end
        end
    end
end
empty_elems = arrayfun(@(s) isempty(s.data) & isempty(s.Rtime),DSout.ChData);
DSout.ChData(empty_elems) = [];
% New addition since change in blc code:
DSout.blc.type='none';
DSout.blc.func=cell(numel(DSout.ChData),1)
for n=1:numel(DSout.ChData)
    if ndims(DSout.ChData(n).data)==3
        DSout.blc.func{n}=cell(size(DSout.ChData(n).data,1),size(DSout.ChData(n).data,3))
    elseif ndims(DSout.ChData(n).data)==2
        DSout.blc.func{n}=cell(size(DSout.ChData(n).data,1))
    end
end

checkIntegrity(DSout);
disp(' ')
disp('Merged dataset has changed in organization:')
disp('* * * * * * * *')
for n=1:length(DSout.ChData)

    disp(['Channel ',num2str(n,'%02i'),': ',DSout.ChData(n).ident]);
end
disp('* * * * * * * *')