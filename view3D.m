function [] = view3D(DS,varargin)
% Plot PDA data
%   
% USEAGE:
%           pdaview(DS,ncontours,ContLine)
%
% INPUTS
%               DS:         Structure containing variables .X, .Abs_wave, and .Rtime
%               ncontours (opt): Step height for Z axis in contour plot. smaller = faster, but less detail
%               ContLine (opt):  LineStyle property
%
%
%
% Notice:
% This mfile is part of the LC-addon for drEEM.
% % Check https://github.com/urbanwuensch/LC-addon-for-drEEM for the latest version
%
%
%
% view3D.m: Copyright (C) 2017 Urban J Wünsch
% Technical University of Denmark
% National Institute of Aquatic Resources
% Section for Oceans and Arctic
% Kemitorvet
% 2800 Kgs. Lyngby, Denmark
% urbw@aqua.dtu.dk
%

if nargin==0
    help view3D
    return
end

disp(' ')
disp(' ')
disp('view3D.m')
disp('----------')
disp('This script will let you plot a dataset of PDA chromatograms')
disp('Press Ctrl+C to cancel...')
disp('----------')
disp('After plotting each contour plot you have the following options:')
disp('1. Press ''next'' to continue with the next sample')
disp('2. Press ''Chromatogram'' to view the absorbance spectrum at the selected retention time')
disp('3. Press ''Spectrum'' to view the chromatogram at the selected wavelength')
disp('4. Press ''Both'' to view the chromatogram at the selected wavelength and the absorbance spectrum at the selected retention time')
disp('5. Use the slider to adjust CLim for the contour plot.')
disp('----------')

% Init
% Argin checks and defs
narginchk(1,3);
numvarargs = length(varargin);

% set defaults for optional inputs
optargs = {60  '-'};

% Overwrite defaults in case user has preference
optargs(1:numvarargs) = varargin;
[steps,ContLine] = optargs{:};


% Initial conversion of Wavelength to index in DS.X
for iCh=1:numel(DS.ChData)
    if ndims(DS.ChData(iCh).data)==3
        idx3D(iCh)=iCh;
    else
        idx3D(iCh)=nan;
    end
end

if numel(find(~isnan(idx3D')))==1
    ChIdx=find(~isnan(idx3D'));
else
    i_count=1;
    for i=idx3D(~isnan(idx3D))
        disp([char(num2str(i)),':',char(DS.ChData(i).ident)]);i_count=i_count+1;
    end
    prompt = 'Select 3D channel to display: ';
    ChIdx= input(prompt);
end
X=getfield(DS.ChData(ChIdx),'data'); % Field extraction

if isfield(DS.ChData(ChIdx),'Rtime')
    Xname='Rtime';
elseif isfield(DS.ChData(ChIdx),'Evol')
    Xname='Evol';
else
    error('Cannot determine Xaxis. Evol or Rtime!?')
end
x=getfield(DS.ChData(ChIdx),Xname);
y=getfield(DS.ChData(ChIdx),'wave');

fig1=figure('name',char(DS.ChData(ChIdx).ident));set(fig1,'InvertHardcopy','off');
for i=1:DS.nSample
    if ~ishandle(fig1); return; end; % Ends function when plot is closed by user
    try
        pda=maxnorm(squeeze(X(i,:,:))');contours=linspace(nanmin(nanmin(pda)),nanmax(nanmax(pda)),steps);
        warning('OFF','MATLAB:contourf:EmptyV6OutputArgument')
       subplot(20,1,1:18);[~,h,~] = contourf(x,y,pda,contours);
      
       warning('ON','MATLAB:contourf:EmptyV6OutputArgument')
       % Formatting
       colorbar;set(h,'LineStyle',ContLine);set(gca,'fontsize',8,'FontName','times new roman')
       xlabel(Xname);ylabel('Wavelength');title(sprintf('%s_%d',char(DS.filelist(i))), 'Interpreter', 'none');
       % uicontrol objects
       uicontrol('Style','text','Position',[100 20 100 20],'String','Further options: ',...
           'Units','normalized','Position', [0.02 0.04 0.15 0.05]);
       uicontrol('Style', 'pushbutton','String','Chromatogram','Callback',{@chromatogram,x,y,pda'},...
           'Units','normalized','Position', [0.18 0.05 0.15 0.05]);
       uicontrol('Style', 'pushbutton','String','Spectrum','Callback',{@spectrum,x,y,pda'},...
           'Units','normalized','Position', [0.35 0.05 0.15 0.05]);
       uicontrol('Style', 'pushbutton','String','Both','Callback',{@ChromaSpectrum,DS,i,getfield(DS.ChData(ChIdx),'wave'),ChIdx},...
           'Units','normalized','Position', [0.55 0.05 0.15 0.05]);
       uicontrol('Style', 'pushbutton','String','Next','Callback',{@b4_callback,DS,i},...
           'Units','normalized','Position', [0.77 0.05 0.15 0.05]);
       
       %ConLim might be a problem here
       ConLim=get(gca,'CLim');
       uicontrol('Style', 'slider','String','CLim',...
       'Min',ConLim(1)+0.01,'Max',ConLim(2)*2,'Value',ConLim(2),...
       'Units','normalized','Position', [0.92 0.2 0.04 0.73],...
       'Callback', {@ContourCLim,ConLim});
       %wait for uiinput on fig1 
       uiwait(gcf)
   catch ME
       disp(' ');rethrow(ME); % Repeats last error message
   end
   
   % In case the last sample is reached, this will just tell the user that the end was reached.
   if i==DS.nSample; close all;end
end
disp(' ')
disp('Done.')
end


%% Functions that get called by uicontrol buttons
function chromatogram(~,~,x,y,X)
            [~,ySel]=ginput(1);
            fig2=figure(2);
            set(gcf, 'units', 'normalized', 'pos', [0.2 0.2 0.4 0.4]);
            set(gcf,'numbertitle','off','name','Chromatogram');
            [h]=plot(x,maxnorm(X(:,knnsearch(y',ySel))))
            hold on
            disp('Close when ready.')
            waitfor(fig2)
end

function spectrum(~,~,x,y,X)
            [xSel,~]=ginput(1);
            fig3=figure(3);
            set(gcf, 'units', 'normalized', 'pos', [0.2 0.2 0.4 0.4]);
            set(gcf,'numbertitle','off','name','PDA spectrum');
            [h]=plot(y,maxnorm(X(knnsearch(x,xSel),:)))
            hold on
            disp('Close when ready.')
            waitfor(fig3);
end


function ChromaSpectrum(~,~,DS,i,wave,ChIdx,Xname)
            [x,y]=ginput(1);
            fig4=figure(4);
            set(gcf, 'units', 'normalized', 'pos', [0.2 0.2 0.7 0.5]);
            set(gcf,'numbertitle','off','name','Chromatogram + PDA spectrum');
            subplot(1,2,1);
            [h]=plot(x,maxnorm(X(:,knnsearch(y',ySel))))
            subplot(1,2,2);
            [h]=plot(y,maxnorm(X(knnsearch(x,xSel),:)))
            disp('Close when ready.')
            waitfor(fig4);
end

function b4_callback(~,~,~,~)
uiresume
end

function ContourCLim(source,~,ConLim)
    val = source.Value;
    set(gca,'CLim',[ConLim(1) val])
end

function [data]=maxnorm(data)
data=data./nanmax(nanmax(data));
end
function [data]=minmaxnorm(data)
data=(data/nanmin(data))/(nanmax(data)-nanmin(data));
end