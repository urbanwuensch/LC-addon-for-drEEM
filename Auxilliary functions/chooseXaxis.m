function [ XaxisValues,XaxisString ] = chooseXaxis( DS,ChSel,pref )
% Find which x-axes are present in your DS and return one to use
% Notice:
% This mfile is part of the LC-addon for drEEM.
% % Check https://github.com/urbanwuensch/LC-addon-for-drEEM for the latest version
%
%
%
% chooseXaxis.m: Copyright (C) 2017 Urban J Wünsch
% Technical University of Denmark
% National Institute of Aquatic Resources
% Section for Oceans and Arctic
% Kemitorvet
% 2800 Kgs. Lyngby, Denmark
% urbw@aqua.dtu.dk
%
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
