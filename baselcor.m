function [DS ] = baselcor( DS,channel,varargin)
% Baseline correction function. Corrections are done based on retention
% time information provided. Interpolation method will then be used to
% estimate and mitigate baseline drift
% 
%USEAGE: 
%       [ DS ] = baselcor( DS,channel,XaxisBase,type,plot_opt,OpMode)
%       [ DS ] = baselcor( DS,varargin )
%
%INPUTS: 
%        DS:         data structure containing EEMs to be modelled in data.X.
%        channel:    Number of components in the model to be fitted.
%        XaxisBase:  List of X-axis values that will be used to fit
%                    baseline function
%        type:       fitting type: 'poly1' (linear), or 'smoothingspline' (piecewise
%                    interpolation). See 'fittype' in fit.m for more
%                    options
%        plot_opt:   If 1, plot will be shown and function paused until user
%                    hits enter
%        OpMode:       "apply" if nothing is specified, alternative: 'reverse'.
%
%OUTPUTS:
%       DS: Corrected dataset
%
%    Additionally, a plot will be shown of the sum of squared errors and 
%    number of iterations before the model converged, for each run. 
%    The run number with the least squares solution will be identified.
%
%Examples
%   DSm=baselcor(DSm,1,[5 39],0)

% Notice:
% This mfile is part of the LC-addon for drEEM.
% % Check https://github.com/urbanwuensch/LC-addon-for-drEEM for the latest version
%
%
%
% baselcor.m: Copyright (C) 2017 Urban J Wünsch
% Technical University of Denmark
% National Institute of Aquatic Resources
% Section for Oceans and Arctic
% Kemitorvet
% 2800 Kgs. Lyngby, Denmark
% urbw@aqua.dtu.dk

%% Function initialization

if nargin==0
    help baselcor
    return
end


narginchk(2,6);
% Xaxis is either Rtime or Evol
try 
    if size(DS.ChData(channel).data,2)==size(DS.ChData(channel).Evol,1)
        Xaxis=DS.ChData(channel).Evol;
    end
end;
try
    if size(DS.ChData(channel).data,2)==size(DS.ChData(channel).Rtime,1)
        Xaxis=DS.ChData(channel).Rtime;
    end
end

if nargin==2
    warning('No baseline fitting parameters supplied. Default values will be used');
    XaxisBase=[min(Xaxis) max(Xaxis)];
    type='poly1';
    plot_opt=0;
    OpMode= 'apply';
elseif nargin==3
    XaxisBase=varargin{1};
    type='poly1';
    plot_opt=0;
    OpMode= 'apply';
elseif nargin== 4
    XaxisBase=varargin{1};
    type=varargin{2};
    plot_opt=0;
    OpMode= 'apply';
elseif nargin == 5
    XaxisBase=varargin{1};
    type=varargin{2};
    plot_opt=varargin{3};
    OpMode= 'apply';
elseif nargin == 6
    XaxisBase=varargin{1};
    type=varargin{2};
    plot_opt=varargin{3};
    OpMode= varargin{4};
end


% Data transfer into scheme with internal format
X.X=DS.ChData(channel).data;
X.Xaxis=Xaxis;
X.nX=numel(X.Xaxis);
X.nSample=DS.nSample;
X.nAbs=size(X.X,3);

selection_l=knnsearch(X.Xaxis,XaxisBase');
selection_r=knnsearch(X.Xaxis,XaxisBase'+0.25);

try blc.func=DS.blc.func; 
catch
    warning('Field ''blc'' not found. That might not be good...');
end
%Waitbar init
h = waitbar(0,'Wait you must...');total=X.nSample*X.nAbs;i=1;
if strcmp(OpMode,'apply')
    % START SAMPLE LOOP
    for s=1:X.nSample
        if all(~isnan(squeeze(X.X(s,:,1))))
        % START WAVELENGTH LOOP
        for n=1:X.nAbs
            % Definition of variables
            y_selected=squeeze(X.X(s,:,n));
            y=nan(X.nX,1);
            x=X.Xaxis;
            % Actual baseline correction
            [y_corrected,func]=baseDetermination(X,y_selected,y,x,selection_l,selection_r,XaxisBase,n,s,plot_opt,type);
            blc.func{channel}{s,n}=func;
            X.X(s,:,n)=y_corrected;
            % Waitbar update
            waitbar(i/total,h,'Hang on. The minions are working...');i=i+1;
            % END WAVELENGTH LOOP
        end
        % END SAMPLE LOOP
        else
            disp(['Sample ',num2str(s),' seems to not contain any information. It was skipped.'])
        end
    end
    close(h)
    DS.ChData(channel).data=X.X;
    DS.blc=blc;
    DS.blc.type=type;
elseif strcmp(OpMode,'reverse')
    for s=1:X.nSample
        for n=1:X.nAbs
            type=DS.blc.type;
            if ndims(X.X)==3
                func=DS.blc.func{channel}{s,n};
                y_selected=squeeze(X.X(s,:,n));
                y=nan(X.nXaxis,1);
                x=X.Xaxis;
                baseline=feval(func, X.Xaxis);
                y=y_selected+baseline';
                X.X(s,:,n)=y';
            end
            
            
            if ~plot_opt==0
                figure(1);hold off
                p1=plot(X.Xaxis,y_selected);hold on;axis tight
                CYLim=get(gca,'YLim');ylim([0.5*CYLim(1) 0.5*CYLim(2)])
                p2=plot(X.Xaxis,y);ylim([0.9*CYLim(1) 0.1*CYLim(2)])
                hold on;legend([p1 p2],'Cor','Reversed (Org)');
                l=line([min(X.Xaxis) max(X.Xaxis)],[0 0]);l.Color='k';l.LineWidth=1;
                if isnumeric(plot_opt)
                    pause(plot_opt);
                else
                    pause
                end
            end
            
            
            
            waitbar(i/total,h,'Hang on. The minions are working...');i=i+1;
        end
    end
    DS.ChData(channel).data=X.X;
    DS.blc.type=['Reversed with ''',type,''' fit'];
end
delete(h)
%checkIntegrity(DS)
disp('Baseline correction performed successfully.')
end


%% Baseline fitting function
function[y_corrected,func]=baseDetermination(X,y_selected,y,x,selection_l,selection_r,XaxisBase,n,s,plot_opt,type)
    % extraction of y-data from chromatogram at selected retention times
    for m=1:size(XaxisBase,2)
        y(selection_l(m),1)=nanmean(y_selected(1,selection_l(m):selection_r(m)));
    end
    % Reduction of data to selected retention times
    x(isnan(y(:,1)),:)=[];
    y(isnan(y(:,1)),:)=[];
    % Actual baseline fit based on selected retention times
    Basefit=fit(x, y, type);
    func=Basefit;
    % Calculation of baseline values across the cromatogram
    baseline=feval(Basefit, X.Xaxis);
    % Baseline-corrected Chromatogram
    y_corrected=y_selected-baseline';
    % If you want to check: Here are the results:
    if ~plot_opt==0
        figure(1);hold off
        p1=plot(X.Xaxis,y_selected);hold on;axis tight
        CYLim=get(gca,'YLim');ylim([0.5*CYLim(1) 0.5*CYLim(2)])
        plot(Basefit,x,y);ylim([0.5*CYLim(1) 0.5*CYLim(2)])
        p2=plot(X.Xaxis,y_corrected);ylim([0.9*CYLim(1) 0.1*CYLim(2)])
        hold on;legend([p1 p2],'Org','Cor');
        l=line([min(X.Xaxis) max(X.Xaxis)],[0 0]);l.Color='k';l.LineWidth=1;
        if isnumeric(plot_opt)
            pause(plot_opt);
        else
            pause
        end
    end
end