function [  ] = checkIntegrity( DS )
% Check integrity of drEEMLC datasets with respect to ndim and size
% USEAGE:
%           [  ] = checkIntegrity( DS )
%
% INPUTS
%            ident: DS in the format of drEEMLC
%
% OUTPUTS
%             none. However, if checks fail, the function produces an error
%             
% Examples:
%       checkIntegrity( DS );
%
% Notice:
% This mfile is part of the LC-addon for drEEM.
% % Check https://github.com/urbanwuensch/LC-addon-for-drEEM for the latest version
%
%
%
% checkIntegrity.m: Copyright (C) 2017 Urban J Wünsch
% Technical University of Denmark
% National Institute of Aquatic Resources
% Section for Oceans and Arctic
% Kemitorvet
% 2800 Kgs. Lyngby, Denmark
% urbw@aqua.dtu.dk
%
% Version 1, January 2016 First version

% Reset warning messages to an empty message
warning('')

%% Perform ndim and size checks for each of the channels individually
for nCh=1:numel(DS.ChData)
    ChData = DS.ChData(nCh); % Create a temporary Ch variable for easy access
    % Check if Channel is either completly populated or  empty
    if strcmp(checkConents(ChData),'FAIL');
        warning('Abs / pres: FAIL');
    end
    if strcmp(chechNaN(ChData),'FAIL');
        warning('NaNs present: FAIL');
    end
    results{nCh,4}=['Abs / pres: ',checkConents(ChData)];
    results{nCh,5}=['NaN: ',chechNaN(ChData)];
    % Rest of checks: nSample, nRtime, nWave, nWave (3rd dimension)
    if ~isempty(ChData.data)
         % nSample Check
        if strcmp(nSampleCheck(ChData,DS.nSample),'PASS')&&strcmp(nSampleCheck(ChData,size(DS.filelist,1)),'PASS')
            results{nCh,1}=([char(num2str(nCh)), ' nSample: PASS...']);
        elseif strcmp(nSampleCheck(ChData,DS.nSample),'FAIL')&&strcmp(nSampleCheck(ChData,size(DS.filelist,1)),'PASS')
            warning([ChData.ident, ' nSample: FAIL...'])
            results{nCh,1}=([char(num2str(nCh)), ' nSample: FAIL...']);
        end
        
        % nRtime Check
        if isfield(DS.ChData,'Rtime')
%             if isfield(DS,'alignment')
%                 nRtime=size(DS.Rtime,1);
%             else 
                nRtime=size(ChData.Rtime,1);
            %end
            if strcmp(nRtimeCheck(ChData,nRtime),'PASS')
                results{nCh,2}=([char(num2str(nCh)), ' nRtime: PASS...']);
            elseif strcmp(nRtimeCheck(ChData,nRtime),'FAIL')&&isempty(strfind(ChData.ident,'[LC Status Trace'))
                warning([ChData.ident, ' nRtime: FAIL...'])
                results{nCh,2}=([char(num2str(nCh)), ' nRtime: FAIL...']);
            end
        end
        % nEvol Check
        if isfield(DS.ChData(1),'Evol')
            if isfield(DS,'alignment')
                nEvol=size(DS.Evol,1);
            else 
                nEvol=size(ChData.Evol,1);
            end
            if strcmp(nRtimeCheck(ChData,nEvol),'PASS')
                results{nCh,2}=([char(num2str(nCh)), ' nEvol: PASS...']);
            elseif strcmp(nRtimeCheck(ChData,nEvol),'FAIL')&&isempty(strfind(ChData.ident,'[LC Status Trace'))
                warning([ChData.ident, ' nEvol: FAIL...'])
                results{nCh,2}=([char(num2str(nCh)), ' nEvol: FAIL...']);
            end
        end
        
        % n3rd Check
        if ndims(ChData.data)==3
            if strcmp(n3rdCheck(ChData,size(ChData.wave,2)),'PASS')
                results{nCh,3}=([char(num2str(nCh)), ' n3rd dim.: PASS...']);
            elseif strcmp(n3rdCheck(ChData,size(ChData.wave,2)),'FAIL')
                warning([ChData.ident, ' n3rd dim.: FAIL...'])
                results{nCh,3}=([char(num2str(nCh)), ' n3rd dim.: FAIL...']);
            end
        else
            results{nCh,3}=[];
        end
    end

end
% Throw error if warning has been issued
if isempty(strfind(lastwarn,'FAIL'))
    %disp('(Integrity of dataset intact)')
elseif ~isempty(strfind(lastwarn,'FAIL'))
    disp(' ')
    disp(' ')
    disp('---------')
    results( all(cellfun(@isempty,results),2), : ) = [];
    disp(results)
    error('Dataset integrity check failed.')
end

end
% END OF CheckIntegrity.m
%% Auxilliary functions

% Function for nSample check
function [result] = nSampleCheck (ChData,nSample)

    if ndims(ChData.data)==3
        dim=1;
    elseif ndims(ChData.data)==2
        dim=1;
    end

    if size(ChData.data,dim)==nSample
        result='PASS';
    else 
        result='FAIL';
    end

end


% Function for nRtime check
function [result] = nRtimeCheck (ChData,nRtime)

    if ndims(ChData.data)==3
        dim=2;
    elseif ndims(ChData.data)==2
        dim=2;
    end

    if size(ChData.data,dim)==nRtime
        result='PASS';
    else 
        result='FAIL';
    end

end



% Function for n3rd dimension check
function [result] = n3rdCheck (ChData,n3rd)

    if ndims(ChData.data)==3
        dim=3;
    elseif ndims(ChData.data)==2
        result='PASS';
        return
    end

    if size(ChData.data,dim)==n3rd
        result='PASS';
    else 
        result='FAIL';
    end

end

% Function for contents empty/populated check
function [result] = checkConents (ChData)
    try;contents(1,1)=~isempty(ChData.md);catch warning('Field md not found');end
    try;contents(2,1)=~isempty(ChData.Rtime);end
    try;contents(2,1)=~isempty(ChData.Evol);end
    contents(3,1)=~isempty(ChData.data);
    contents(4,1)=~isempty(ChData.wave);

    if sum(contents)==4&&~isempty(ChData.wave)
        result='PASS';
    elseif sum(contents)==3&&isempty(ChData.wave)
        result='PASS';
    elseif sum(contents)==0
        result='PASS';
    else
        result='FAIL';
    end

end

% Function for contents empty/populated check
function [result] = chechNaN (ChData)
    if any(isnan(reshape(ChData.data,[],1)));
        result='FAIL';
    else
        result='PASS';
    end
end


%% EOF