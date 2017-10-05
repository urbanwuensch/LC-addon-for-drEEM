function [DS]=importrun(ident)
% Import HPLC data exported by LabSolutions v5.80
%   
% USEAGE:
%           [DS]=aqy(DS,mode,method,qy_Std,qy_std_choice)
%
% INPUTS
%            ident: identifier for files to be imported
%
% OUTPUTS
%               DS: data structure containing HPLC data
%             
% Examples:
%       DS=importHPLC('*.txt');
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
% importrun.m: Copyright (C) 2017 Urban J Wünsch
% Technical University of Denmark
% National Institute of Aquatic Resources
% Section for Oceans and Arctic
% Kemitorvet
% 2800 Kgs. Lyngby, Denmark
% urbw@aqua.dtu.dk
%
% Version 1, January 2016 First version

% Welcome messages
disp(' ')
disp(' ')
disp(' ')
disp('importHPLC.m')
disp('------------')
disp('Reading in metadata...')
if exist('str2doubleq');opt='fast';else;opt='slow';end

% Scan for importable files
Folder = pwd; % Folder will be current working directory
files   = dir(fullfile(Folder, ident)); % scan for all files containing ident in name
if size(files,1)==0;error('Nothing to import...');end % error in case no files are found

%% Identification of importable data + fetching of metadata
mdImport = textread(fullfile(Folder,files(1).name),'%s','delimiter','\n') ;
emptyCells=find(cellfun('isempty',mdImport));mdImport=mdImport(~cellfun('isempty',mdImport));

EmptyCount=0;
for k=1:numel(mdImport)
    if isempty(mdImport{k})
        metadata{k}='';
        EmptyCount=EmptyCount+1;
    else
    mdProcessed{k,:}=textscan(mdImport{k},'%s%s%[^\n\r]', 'Delimiter', ',','ReturnOnError', false,'EmptyValue',NaN);
    metadata{k,1}=mdProcessed{k,1}{1,1}{1,1};metadata{k,2}=mdProcessed{k,1}{1,2};
    end
end;EmptyCount=EmptyCount-1;
clearvars mdImport mdProcessed

% New addition: Fetch channel list automatically for import! Leave no one behind ;)
hit=1;
for n=1:size(metadata,1)
    if ~isempty(strmatch('[PD',metadata{n,1}))|| ~isempty(strmatch('[L',metadata{n,1}))
        Ch{hit,1}=metadata{n,1};hit=hit+1;
    end
end


for n=1:length(Ch)
    Channels(n,:)=ChIdent(metadata,Ch{n,1});
    % Extract wavelengths in case of 3D data
    if strcmp(Ch(n,:),'[PDA 3D]')&&~isempty(Channels(n).md)
        % Extraction of Wavelenths
        %offset=sum(Channels(n,:).idx(1,1)>emptyCells);
        %wavemat=dlmread(fullfile(Folder,files(1).name),',',[Channels(n,:).idx(1,1)-2+offset 1, Channels(n,:).idx(1,1)-2+offset, Channels(n,:).idx(2,2)]);
        %wavemat=wavemat/100;Channels(n,:).wave=wavemat;clearvars wavemat
    end
end
clearvars Ch ident n k

% 
UserInput=0;
for n=1:size(Channels)
    if strmatch('[PDA 3D]',Channels(n).ident)
        cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
        % Based on the number of wavelengths this will guess if '[PDA 3D]' tag is PDA or FLD data
        % FLD data usually has <100 wavelength points, PDA has >200
        if Channels(n).md.Dvalues(find(cellfun(cellfind('# of Wavelength Axis Points'),Channels(n).md.Dtype)))<100
            UserInput=2;
            disp('It was assumed that the files contain 3D-FLD data, NOT 3D-PDA data.')
        else
            UserInput=1;
            disp('It was assumed that the files contain 3D-PDA data, NOT 3D-FLD data.')
        end
        idx3d=n;
    end
end

disp(' ')
disp(['Importing ',num2str(length(files)),' samples...'])
%% Data import.
for k = 1:length(files)
    disp(files(k).name)
    
    % Read in data
    % This function is slow as hell for lots of lines (Channel data!), look for a faster alternative!
    [ dataAll,data2D ] = readHPSECdata( files(k).name);
    for i=1:length(Channels)
    % Import of 2D data
    if ~strcmp(Channels(i).ident,'[PDA 3D]')&&~isempty(Channels(i).idx)
        if strcmp(opt,'fast')
           RtimeVal=real(str2doubleq(data2D(Channels(i).idx(1,1):Channels(i).idx(1,2),1)));
           tempData=real(str2doubleq(data2D(Channels(i).idx(1,1):Channels(i).idx(1,2),2)));
        elseif strcmp(opt,'slow')
           RtimeVal=str2double(data2D(Channels(i).idx(1,1):Channels(i).idx(1,2),1));
           tempData=str2double(data2D(Channels(i).idx(1,1):Channels(i).idx(1,2),2));
        end
        if  all(round(Channels(i).Rtime,4)==round((RtimeVal),4))
            Channels(i).data(k,:)=tempData;clearvars temp
        else
            csvwrite('ErrorLog_Rtime_.csv',[Channels(i).Rtime RtimeVal])
            error(['Retention time mismatch between 1st and ',num2str(k),'th sample: ','(',files(k).name,')'])
        end
    % Import of 3D data
    elseif strcmp(Channels(i).ident,'[PDA 3D]')&&~isempty(Channels(i).idx)
        if all(Channels(i).Rtime==str2double(data2D(Channels(i).idx(1,1):Channels(i).idx(1,2),1)))
            RowIdx=[Channels(i).idx(1,1):Channels(i).idx(1,2)];
            X=nan(Channels(i).idx(1,2)-Channels(i).idx(1,1),numel([Channels(i).idx(2,1):Channels(i).idx(2,2)])+1);
            wave=cellfun(@str2double,dataAll{RowIdx(1)-1})/100;
            wave=wave(Channels(i).idx(2,1):end);
            if k==1
                Channels(i).wave=wave;
            elseif ~all(Channels(i).wave==wave)
                error('Wavelength dimension problem during import!')
            end
            if strcmp(opt,'fast')
                for m=1:numel(RowIdx)
                    temp=real(cellfun(@str2doubleq,dataAll{RowIdx(m)})); 
                    X(m,:)=temp(:,Channels(i).idx(2,1):end);
                end
            elseif strcmp(opt,'slow')
                if k==1;warning('Download str2doubleq to speed this code up!');end
                for m=1:numel(RowIdx)
                    temp=cellfun(@str2double,dataAll{RowIdx(m)});
                    X(m,:)=temp(:,Channels(i).idx(2,1):end);
                end
            end
            Channels(i).data(k,:,:)=X;
        else
            error(['Retention time mismatch between 1st and ',num2str(i),'th sample.'])
        end
    end
    end
    filename=strsplit(char(data2D(find(strcmp(data2D(:,1),'Data File Name')),2)),'\');
    filename=filename(1,end);
    filelist{k,:}=filename{1,1};
    clearvars data2D pdamat
end



%% Final allocation of DS

% Shimadzu is dumb, so the 2D Emission profiles are named wrong. This is the fix
if UserInput==2
   Channels(idx3d).ident='[2D Emission Scan(Detector A-3D)]';
end

DS.filelist=filelist;
DS.nSample=k;

if UserInput==1
    if Channels(idx3d).md{find(strcmp(Channels(idx3d).md{:,1},'# of Wavelength Axis Points')),2}==size(squeeze(Channels(idx3d).data(1,1,:)),1)
        DS.nAbs=size(squeeze(Channels(idx3d).data(1,1,:)),1);
    end
end
DS.ChData=Channels;

for iCh=1:numel(DS.ChData)
    DS.ChData(iCh).data(isnan( DS.ChData(iCh).data)) = 0;
end

DS.ChData=rmfield(DS.ChData,'idx');
DS = orderfields(DS);
checkIntegrity(DS);


fclose('all');
disp(' ')
disp('----------------')
disp('Import complete.')
end

%% END OF importHPLC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ ChOut ] = ChIdent( metadata,channel )

if ~isempty(find(strcmp(metadata(:,1),channel)));
    disp([channel,' found']);
    ChOut.ident=channel;
    ChOut.md=MDextract(metadata,channel);
    ChOut.idx = indexQuery( metadata,ChOut.md,channel );
    ChOut.Rtime=str2double(metadata(ChOut.idx(1,1):ChOut.idx(1,2),1));
    ChOut.wave=[];
    ChOut.Evol=nan(size(ChOut.Rtime));
else
    ChOut.ident=channel;
    ChOut.md=[];
    ChOut.idx =[];
    ChOut.Rtime =[];
    ChOut.wave=[];
    ChOut.Evol=nan(size(ChOut.Rtime));
end
end
%% End of function

function [idxOut] = indexQuery( data,metadata,target )
% Look for row indicies

start=find(strcmp(data(:,1),target));
for i=start:start+25
    if ~isnan(str2double(data(i,1)))
    lastIdx=i-1;
    break
    end
end
if ~strcmp(target,'[PDA 3D]')
    Vlast=lastIdx+metadata{find(strcmp(metadata{:,1},'# of Points')),2};
    C(1,1)=2;C(1,2)=2;
elseif strcmp(target,'[PDA 3D]')
    Vlast=lastIdx+metadata{find(strcmp(metadata{:,1},'# of Time Axis Points')),2};
    C(1,1)=2;
    C(1,2)=metadata{find(strcmp(metadata{:,1},'# of Wavelength Axis Points')),2};
end


if ~round(str2double(data{Vlast,1}),3)==round(metadata{find(strcmp(metadata{:,1},'End Time(min)')),2},3)
    error(['Error during index determination of import of ',target,'. Missmatch at end!'])
end
if ~round(str2double(data{lastIdx+1,1}),3)==round(metadata{find(strcmp(metadata{:,1},'Start Time(min)')),2},3)
    error(['Error during index determination of import of ',target,'. Missmatch at start!'])
end


R(1,1)=lastIdx+1;
R(1,2)=Vlast;

idxOut=[R;C];

end

%% End of function

function [ dataOut ] = MDextract( data,target )
%Extraction of metadata into data given target




start=find(strcmp(data(:,1),target));

for i=start:start+25
    if ~isnan(str2double(data(i,1)))
    last=i-1;
    break
    end
end

md=data(start:last,:);

hit=0;
for n=1:size(md,1)
    if ~isnan(str2double(md{n,2}))
        hit=hit+1;
    Dtype(hit,1)=md(n,1);
    Dvalues(hit,1)=str2double(md{n,2});
    end
end
dataOut=table(Dtype,Dvalues);


end
% End of function
%% Function for reading all ASCII-data of HPSEC chromatograms
function [ dataAll,data2D ] = readHPSECdata( filename )
    baseFileName = filename; 
    
    Nrows = numel(textread(fullfile(pwd,baseFileName),'%1c%*[^\n]'));
    dataAll=cell(Nrows,1);
    theLine=cell(Nrows+1,1);
    fid = fopen(fullfile(pwd,baseFileName),'r');

   
    if fid==-1
        error(['Couldn''t open ',filename]);
    end
    clearvars theLine
    theLine{1,1} = fgetl(fid);
    n=2;
    while ~feof(fid)
        theLine{n,1} = fgetl(fid);
        if ~isempty(theLine{n,1})
            dataAll{n,1}=regexp(theLine{n,1}, regexptranslate('escape', ','), 'split');
            n=n+1;
        end
    end
    fclose(fid);
    
    % Preps for 2D data
    hplc2D=cell(1,2);hplc2D{1}=cell(Nrows,1);hplc2D{2}=cell(Nrows,1);
    data2D=cell(Nrows,2);
    for i=2:numel(dataAll)
        hplc2D{1}{i,1}=dataAll{i}{1};
        try hplc2D{2}{i,1}=dataAll{i}{2};end % This needs to be handeled in case it's an empty cell in the second row
    end
    data2D(:,1)=hplc2D{1,1};
    data2D(:,2)=hplc2D{1,2};
end
% END