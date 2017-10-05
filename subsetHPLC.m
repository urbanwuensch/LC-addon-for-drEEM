function [ DS ] = subsetHPLC(DS , delSamples , deltR, Xname )
% delete samples and / or retention time from HPLC datasets
%   
% USEAGE:
%           [ DS ] = subsetHPLC(DS , delSamples , deltR )
%
% INPUTS
%            DS:            Dataset do be processed
%            delSamples:    Samples to be deleted
%            deltR:         Retention times to be removed. 
%                           E.g. [5 6;25 65] will remove minutes 5-6 and 25-65
%                           Since knnsearch is used, the retiontimes can be out of bounds,
%                           if so, the nearest neighbour will be found and used
%
% OUTPUTS
%               DS: data structure containing the subset HPLC data
%
% Examples:
%       DSout=subsetHPLC(DS,[1:2],[5 6;25 65])
%
% Notice:
% This mfile is part of the LC-addon for drEEM.
% % Check https://github.com/urbanwuensch/LC-addon-for-drEEM for the latest version
%
%
%
% subsetHPLC.m: Copyright (C) 2017 Urban J Wünsch
% Technical University of Denmark
% National Institute of Aquatic Resources
% Section for Oceans and Arctic
% Kemitorvet
% 2800 Kgs. Lyngby, Denmark
% urbw@aqua.dtu.dk
%
% Version 1, June 2016 First version


%% Function init
if nargin==0
    help interplXaxis
    return
elseif nargin<4
    error('Not enough input arguments.')
end


if ~isempty(deltR)
    if ~isfield(DS.ChData,Xname)
        if ~isfield(DS.ChData,Xname)
            error('Xaxis not found. Take a look in your DS and see if the axis you want to subset really exists')
        end
    end
end
% [ ~,XaxisString ] = chooseXaxis( DS,1,'Rtime'); % Channel shouldn't matter
% 
% 
% if ~strcmp(XaxisString,Xname)
%     disp(XaxisString)
%     disp(Xname)
%     error('Xaxis not found. Take a look in your DS and see if the axis you want to subset really exists')
% end


DS.nSample=DS.nSample-numel(delSamples);
DS.filelist(delSamples)=[];
% Deletion of sampels
for nCh=1:length(DS.ChData)
    DS.ChData(nCh).data(delSamples,:,:)=[];
end

%% Deletion of Xaxis
if ~isempty(deltR)
    for nCh=1:length(DS.ChData)
        [ XaxisVal,XaxisStr ] = chooseXaxis( DS,nCh,Xname);
        for ndeltR=1:size(deltR,1)
            idx1=knnsearch(XaxisVal,deltR(ndeltR,1));
            idx2=knnsearch(XaxisVal,deltR(ndeltR,2));
            sample=DS.ChData(nCh).data;
            
            if ndims(sample)==3
                sample(:,idx1:idx2,:)=[];
            elseif ndims(sample)==2
                sample(:,idx1:idx2)=[];
            end
            XaxisVal(idx1:idx2)=[];
            
            DS.ChData(nCh)=setfield(DS.ChData(nCh),'data',sample);
            DS.ChData(nCh)=setfield(DS.ChData(nCh),XaxisStr,XaxisVal);
            
            if strcmp(XaxisStr,'Evol')&&isfield(DS.ChData(nCh),'Rtime')
                Rtime=getfield(DS.ChData(nCh),'Rtime');
                Rtime(idx1:idx2)=[];
                DS.ChData(nCh)=setfield(DS.ChData(nCh),'Rtime',Rtime);
            elseif strcmp(XaxisStr,'Rtime')&&isfield(DS.ChData(nCh),'Evol')
                Evol=getfield(DS.ChData(nCh),'Evol');
                Evol(idx1:idx2)=[];
                DS.ChData(nCh)=setfield(DS.ChData(nCh),'Evol',Evol);
            end
        end
        clearvars XaxisVal XaxisStr Evol Rtime
    end
end
checkIntegrity(DS);

end

function [ XaxisValues,XaxisString ] = chooseXaxis( DS,ChSel,pref )
% Little function that will find which x-axes are present in your DS and return one to use

XaxisValues=[];
XaxisString=[];

potXaxisStrings{1}='Rtime';   potXaxisUnit{1}='Rtime';
potXaxisStrings{2}='Evol';    potXaxisUnit{2}='Evol';

dataDims=size(DS.ChData(ChSel).data); % dataDims(2) is the Xaxis
try
    if strcmp(pref,'Rtime')
        order=numel(potXaxisStrings):-1:1;
    elseif strcmp(pref,'Evol')
        order=1:numel(potXaxisStrings);
    end
catch
    order=1:numel(potXaxisStrings);
end
for n=order
    if isfield(DS.ChData,potXaxisStrings{n})
        if numel(getfield(DS.ChData(ChSel),potXaxisStrings{n}))==dataDims(2)
            XaxisValues=getfield(DS.ChData(ChSel),potXaxisStrings{n});
            XaxisString=potXaxisUnit{n};
        end
    end
end

if isempty(XaxisValues)
    error(['None of the Xaxis found. Channel idx #', num2str(ChSel)])
end

end
