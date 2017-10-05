function [ drEEMout ] = hplc2eem(DS,Ex,absopt)
% Convert HPLC 2D emission scans and absorbance data (drEEMLC format) into drEEM compatible EEMs
%   
% USEAGE:
%           [ drEEMout ] = HPLC2drEEM(DS,Ex,absopt)
%
% INPUTS
%               DS:       data structure containing the 2D scans and absorbance data
%               Ex:       Excitation wavelengths, e.g. [240:10:450]
%               absopt:   if 1, absorbance will be included, if zero, absorbance will not be included
%
% OUTPUTS
%               drEEMout: drEEM compatible dataset
%
% NOTES
%           1. The function assumes that the absorbance data is located in in the proper channel
%              and is located in the first sample of that channel (DS.ChData.data(1,:,:))
%           2. Please make sure that every DS.filelist entry contains a reference to the ex wavelength.
%              This is to make sure that the correct Ex wavelengths are assigned to a file
%             
% Examples:
%       [ Xin ] = HPLC2drEEM(DS,[250:10:450],1);
%
% Notice:
% This mfile is part of the LC-addon for drEEM.
% % Check https://github.com/urbanwuensch/LC-addon-for-drEEM for the latest version
%
%
%
% hplc2eem.m: Copyright (C) 2017 Urban J Wünsch
% Technical University of Denmark
% National Institute of Aquatic Resources
% Section for Oceans and Arctic
% Kemitorvet
% 2800 Kgs. Lyngby, Denmark
% urbw@aqua.dtu.dk
%
% Version 1, January 2016 First version
%%
if nargin==0
    help hplc2eem
    return
end

disp(' ')
disp('Conversion of HPLC chromatogram data into drEEM scheme')
disp('* * * * * * ')

%% Locate Abs and Flu channel
for i=1:numel(DS.ChData)
    if ~isempty(strfind(DS.ChData(i).ident,'2D Emission'))
        EmScanCh=i;
    end
    if ~isempty(strfind(DS.ChData(i).ident,'[PDA 3D]'))&&absopt~=0
        AbsCh=i;
    end
end

%% Rtime and / or Evol there?
if isfield(DS.ChData(EmScanCh),'Rtime')
    Xname='Rtime';
elseif isfield(DS.ChData(EmScanCh),'Evol')
    Xname='Evol';
else
    error('Can''t find the Filename axis. It should be ''Rtime'' or ''Evol''')
end
    

Xaxis=getfield(DS.ChData(EmScanCh),Xname);
if absopt~=0
    if numel(Xaxis)~=numel(getfield(DS.ChData(AbsCh),Xname))
        error('size of Xaxis unequal (Fl and Abs)')
    end
end

%% Tranfer of misc data
nSample=size(DS.ChData(EmScanCh).data,2);
% prompt = 'Type filename prefix: ';
% filename_pre = input(prompt,'s');
for i=1:nSample
    %filelist{i,1}= [filename_pre,char(num2str(Rtime(i)))];
    filelist{i,1}= [char(num2str(Xaxis(i)))];
end


%% Transfer of fluorescence data into drEEM scheme
nEm=size(DS.ChData(EmScanCh).data,3);
nEx=size(Ex,2);
Em=DS.ChData(EmScanCh).wave';
X=nan(nSample,nEm,nEx);

for Ex_i=1:size(Ex,2)
    for Em_i=1:size(DS.ChData(EmScanCh).data,3)
        X(:,Em_i,Ex_i)=DS.ChData(EmScanCh).data(Ex_i,:,Em_i);
    end
end

%% Transfer of absorbance data into drEEM scheme
if absopt==1
    nAbs=size(DS.ChData(AbsCh).data,3);
    Abs_wave=DS.ChData(AbsCh).wave;



    nAbs=size(DS.ChData(AbsCh).data,3);
    Abs_wave=DS.ChData(AbsCh).wave;

    Abs=squeeze(DS.ChData(AbsCh).data(1,:,:));

    Abs_wave_i=[2.*round(ceil(nanmin(Abs_wave))/2):round(mean(diff(Abs_wave))):floor(nanmax(Abs_wave))];
    nAbs_i=size(Abs_wave_i,2);

    for i=1:nSample
        Abs_i(i,:)=interp1(Abs_wave',Abs(i,:)',Abs_wave_i')';
    end

    Abs_wave=Abs_wave_i;
    Abs=Abs_i;

    Abs_ife=nan(size(Abs));
    Abs_ife(1,:)=Abs_wave;
    Abs_ife(2:nSample+1,:)=Abs;

    Abs_ife(isnan(Abs_ife)) = 0 ;
    Abs(isnan(Abs)) = 0 ;
else
    disp('Absorbance not imported.')
end
%% Definition of final export variable drEEMout

drEEMout.nSample=nSample;
drEEMout.filelist=filelist;
drEEMout.i=[1:i]';
drEEMout=setfield(drEEMout,Xname,Xaxis);

%X=X./10E5;
%disp('All intensity values have been divided by 10E5.')
disp(['The max. FL intensity is ',num2str(max(max(max(X))))])
disp('If this number is >1,000, consider dividing by a factor')
drEEMout.X=X;
drEEMout.Em=Em;
drEEMout.Ex=Ex';
drEEMout.nEm=nEm;
drEEMout.nEx=nEx;

if absopt==1
    drEEMout.Abs=Abs;
    drEEMout.Abs_ife=Abs_ife;
    drEEMout.Abs_wave=Abs_wave;
end


%% Safety measure: Assume Ex is mentioned in the filenames and check that things are in the correct order
for i=1:DS.nSample
    k{i}=strfind(DS.filelist{i},num2str(Ex(i)));
    if isempty(k{i})
        error('It appears your list is not properly sorted! Check filenames and try again')
    end
end
if all(cell2mat(k))
    disp('All Ex wavelength files seem to be sorted in the correct order.')
end


drEEMout = orderfields(drEEMout);
disp(' ')
disp('Conversion complete')
disp('* * * * * * ')
