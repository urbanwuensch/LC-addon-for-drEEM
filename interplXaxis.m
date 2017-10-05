function [ DSout ] = interplXaxis( DS,Xaxis,XaxisName,OpMode,alignCh )
% 
% DS:           dataset
% Xaxis:        Xaxis values to be assigned to Xaxis in DS
% XaxisName:    Name of Xaxis in DS 
% OpMode:         'interpolation' or 'alignment', default is 'interpolation'
% alignCh:      if OpMode is 'alignment', you can specify if you want to only align a selected channel (default: all)
%
%
% USAGE:        [ DSout ] = interplXaxis( DS,Xaxis,XaxisName,OpMode )
% NOTE on OpMode: 'interpolation' will use Xaxis and interpolate data so it fits Xaxis
%               'alignment' will assume Xaxis are the true values for DS and interpolate to fit DS' xaxis
%                This is useful, if you want to shift values by a fixed amount and then interpolate to old values
%               'alignment' needs equal Xaxis dimensions across all channels (you could run interpolation first ;) )
%
% Notice:
% This mfile is part of the LC-addon for drEEM.
% % Check https://github.com/urbanwuensch/LC-addon-for-drEEM for the latest version
%
%
%
% interplXaxis.m: Copyright (C) 2017 Urban J Wünsch
% Technical University of Denmark
% National Institute of Aquatic Resources
% Section for Oceans and Arctic
% Kemitorvet
% 2800 Kgs. Lyngby, Denmark
% urbw@aqua.dtu.dk
%


%% Function init
if nargin==0
    help interplXaxis
    return
elseif nargin<3
    error('Not enough input arguments.')
end

if nargin==3
    OpMode= 'interpolation';
end

if ~isfield(DS.ChData(1), XaxisName)
    error('XaxisName does not exist in DS.')
end

try
    checkIntegrity(DS)
catch
    error('DS is not correctly formated. Run checkIntegrity(DS) for details')
end

DSout=DS; %
if size(Xaxis,1)>size(Xaxis,2);Xaxis=Xaxis';end

if strcmp(OpMode,'alignment')&&nargin==4
    alignCh=1:numel(DS.ChData);
end
disp(['Mode: ',OpMode]);

%% Interpolations
if strcmp(OpMode,'interpolation')
    for n=1:numel(DS.ChData)
        if ndims(DS.ChData(n).data)==2
            x=getfield(DS.ChData(n), XaxisName);
            yIt=nan(DS.nSample,size(Xaxis,2));
            for i=1:DS.nSample
                Sample=DS.ChData(n).data(i,:);
                yIt(i,:)=interp1(x,Sample,Xaxis,'pchip');
                clearvars Sample
            end
            DSout.ChData(n).data=yIt;        
            clearvars x yIt

        elseif ndims(DS.ChData(n).data)==3
             Y=getfield(DS.ChData(n), XaxisName);
             X=DS.ChData(n).wave;
             Yq=Xaxis;
             Xq=DS.ChData(n).wave;
             
             [X,Y]=meshgrid(X,Y);
             [Xq,Yq]=meshgrid(Xq,Yq);
             Sq=nan(DS.nSample,size(Xaxis,2),size(Y,2));
            for i=1:DS.nSample
                S=squeeze(DS.ChData(n).data(i,:,:));
                Sq(i,:,:)=interp2(X,Y,S,Xq,Yq,'spline');
                clearvars S
            end
            DSout.ChData(n).data=Sq;
            clearvars X Y yIt Xq Yq
        end

        if size(Xaxis,1)>size(Xaxis,2)
            DSout.ChData(n)=setfield(DSout.ChData(n),XaxisName,Xaxis);
        else
            DSout.ChData(n)=setfield(DSout.ChData(n),XaxisName,Xaxis');
        end
    end
elseif strcmp(OpMode,'alignment')
        for n=alignCh
            if ndims(DS.ChData(n).data)==2
                x=getfield(DS.ChData(n), XaxisName);
                yIt=nan(DS.nSample,size(Xaxis,2));
                for i=1:DS.nSample
                    Sample=DS.ChData(n).data(i,:);
                    yIt(i,:)=interp1(Xaxis,Sample,x,'pchip');
                    clearvars Sample
                end
                DSout.ChData(n).data=yIt;        
                clearvars x yIt

            elseif ndims(DS.ChData(n).data)==3
                 Y=Xaxis;
                 X=DS.ChData(n).wave;
                 Yq=getfield(DS.ChData(n), XaxisName);
                 Xq=DS.ChData(n).wave;

                 [X,Y]=meshgrid(X,Y);
                 [Xq,Yq]=meshgrid(Xq,Yq);
                 Sq=nan(DS.nSample,size(Xaxis,2),size(Y,2));
                for i=1:DS.nSample
                    S=squeeze(DS.ChData(n).data(i,:,:));
                    Sq(i,:,:)=interp2(X,Y,S,Xq,Yq,'spline');
                    clearvars S
                end
                DSout.ChData(n).data=Sq;
                clearvars X Y yIt Xq Yq
            end
        end
end

% Some final things, remove other Xaxis if present, since they are simply incorrect now:
if strcmp(XaxisName,'Rtime')
    try DSout.ChData=rmfield(DSout.ChData,'Evol');end
elseif strcmp(XaxisName,'Evol')
    try DSout.ChData=rmfield(DSout.ChData,'Rtime');end
end
% Final check of dataset integrity
try
    checkIntegrity(DSout)
catch
    error('result dataset is not consistent. This is weird')
end