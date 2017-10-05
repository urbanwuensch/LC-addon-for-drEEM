function [ DS ] = Rtime2Evol( DS, tR , tflow )
% Use Retention times and supplied flow profile to calculate elution
% volumes
%   
% USEAGE:
%           [ DS ] = Rtime2Evol( DS, tR , tflow )
%
% INPUTS
%            DS:    Dataset containing channel data
%            tR:    Retention times of profile (like specified in HPLC software)
%            tflow: Respective flow values belonging to tR values.
%
% OUTPUTS
%               DS: data structure containing elution volume in channel
%               data
%             
% Examples:
%       DS=Rtime2Evol(DS,[0 32.5 36 43 44.5 62.5],[0.2 0.2 0.25 0.25 0.2 0.2]);
%
% Notice:
% This mfile is part of the LC-addon for drEEM.
% % Check https://github.com/urbanwuensch/LC-addon-for-drEEM for the latest version
%
%
%
% Rtime2Evol.m: Copyright (C) 2017 Urban J Wünsch
% Technical University of Denmark
% National Institute of Aquatic Resources
% Section for Oceans and Arctic
% Kemitorvet
% 2800 Kgs. Lyngby, Denmark
% urbw@aqua.dtu.dk
%
%
% Version 1, September 2016 First version

for nCh=1:numel(DS.ChData)
    Rtime=DS.ChData(nCh).Rtime;
    flow=nan(size(Rtime));
    
    for ntR=2:1:numel(tR)

        flfit=fit([Rtime(knnsearch(Rtime,tR(ntR-1))) Rtime(knnsearch(Rtime,tR(ntR)))]', [tflow(ntR-1) tflow(ntR)]', 'pchip');
        flline=feval(flfit, Rtime(knnsearch(Rtime,tR(ntR-1)):knnsearch(Rtime,tR(ntR))));
        flow(knnsearch(Rtime,tR(ntR-1)):knnsearch(Rtime,tR(ntR)))=flline;
    end
    Evol(1,1)=Rtime(1,1)*flow(1,1);
    for n=2:numel(Rtime)
        Evol(n,1)=Evol(n-1)+(Rtime(n,1)-Rtime(n-1,1))*flow(n,1);
    end
    DS.ChData(nCh).Evol=Evol;
    clearvars flow flfit flline Evol Rtime
end

checkIntegrity(DS);