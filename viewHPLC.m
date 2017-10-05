function [] = viewHPLC( DS,ChSel,SampleSel )
% function for the viewing of 2D HPLC data ~=PDA or Flu scans!
%   
% USEAGE:
%           viewHPLC( DS , ChSel, SampleSel )
%
% INPUTS
%            DS:              HPSEC-toolbox formatted dataset
%            ChSel:           numeric or character, e.g. 1 or 'B-Ch1'
%                             if ndims(ChSel)>2, user will be asked for input
%                             to specify wavelength selection
%            SampleSel:       optional, select samples to plot, default: all
%
% OUTPUTS
%            none. Well, the plots ;)
%             
% Examples:
%       	viewHPLC( DS , 'B-Ch1', [1:DS.nSample] )
%
% Notice:
% This mfile is NOT part of anything (yet).
%
% NOTES:    All files have to be identical in dimension. Please exclude optional data
%           such as peak tables that could alter files between runs! 
%
%
% Notice:
% This mfile is part of the LC-addon for drEEM.
% % Check https://github.com/urbanwuensch/LC-addon-for-drEEM for the latest version
%
%
%
% viewHPLC.m: Copyright (C) 2017 Urban J Wünsch
% Technical University of Denmark
% National Institute of Aquatic Resources
% Section for Oceans and Arctic
% Kemitorvet
% 2800 Kgs. Lyngby, Denmark
% urbw@aqua.dtu.dk
%
% Version 2, Dec 2016 added functionality to select channels with string
% Version 1, Feb 2016 initial draft

%% Function init
if nargin==0
    help viewHPLC
    return
elseif nargin<2
    error('Not enough input arguments.')
end

if nargin==2
    SampleSel=[1:DS.nSample];
end

if ~isnumeric(ChSel)
    ChSel=find(~cellfun(@isempty,(arrayfun(@(n) strfind(DS.ChData(n).ident, ChSel), 1:numel(DS.ChData),'UniformOutput',false))))
    if numel(ChSel)>1
        disp('Which channel did you mean? Type number!')
        for n=1:numel(ChSel)
            disp([num2str(ChSel(n)),': ',DS.ChData(ChSel(n)).ident])
        end
    prompt = '';
    ChSel= input(prompt);
    end
end
[ XaxisValues,XaxisString ] = chooseXaxis( DS,ChSel );

% Check if user selected 3D-Ch
if ndims(DS.ChData(ChSel).data)>2
    disp('Selected Channel has two modes, please specify wavelength:')
    prompt = ['[',num2str(min(DS.ChData(ChSel).wave)),' - ',num2str(max(DS.ChData(ChSel).wave)),']'];
    UserWave= input(prompt);
    
    f1=figure('InvertHardcopy','off','Color',[1 1 1]);
    set(gcf, 'units', 'normalized', 'pos', [0.1 0.2 0.8 0.4])

    plot(XaxisValues,squeeze(DS.ChData(ChSel).data(SampleSel,:,knnsearch(DS.ChData(ChSel).wave',UserWave))))
    xlabel(XaxisString);
    title([DS.ChData(ChSel).ident,' at '...
        ,num2str(DS.ChData(ChSel).wave(knnsearch(DS.ChData(ChSel).wave',UserWave))),...
       'nm'])

    if numel(SampleSel)<10
        legend(DS.filelist(SampleSel),'Interpreter', 'none','location','eastoutside')
    else
        disp('Too many lines, legend omitted.')
    end
    axis tight
    
else

    f1=figure('InvertHardcopy','off','Color',[1 1 1]);
    set(gcf, 'units', 'normalized', 'pos', [0.1 0.2 0.8 0.4])

    plot(XaxisValues,DS.ChData(ChSel).data(SampleSel,:))
    xlabel(XaxisString);
    title(DS.ChData(ChSel).ident)

    if numel(SampleSel)<10
        legend(DS.filelist(SampleSel),'Interpreter', 'none','location','eastoutside')
    else
        disp('Too many lines, legend omitted.')
    end
    axis tight
end

end

function [ XaxisValues,XaxisString ] = chooseXaxis( DS,ChSel )
% Little function that will find which x-axes are present in your DS and return one to use

potXaxisStrings{1}='Rtime';potXaxisUnit{1}='Retention time [min]';
potXaxisStrings{2}='Evol';potXaxisUnit{2}='Elution volume [mL]';

dataDims=size(DS.ChData(ChSel).data); % dataDims(2) is the Xaxis
for n=1:numel(potXaxisStrings)
    if isfield(DS.ChData,potXaxisStrings{n})
        if numel(getfield(DS.ChData,potXaxisStrings{n}))==dataDims(2)
            XaxisValues=getfield(DS.ChData,potXaxisStrings{n});
            XaxisString=potXaxisUnit{n};
        end
    end
end
end