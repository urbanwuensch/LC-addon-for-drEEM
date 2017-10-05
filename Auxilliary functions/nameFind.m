function [idx ] = nameFind( DS,string )
%Simply returns the indices of stringmatches in DS.filelist
%
% Notice:
% This mfile is part of the LC-addon for drEEM.
% % Check https://github.com/urbanwuensch/LC-addon-for-drEEM for the latest version
%
%
%
% nameFind.m: Copyright (C) 2017 Urban J Wünsch
% Technical University of Denmark
% National Institute of Aquatic Resources
% Section for Oceans and Arctic
% Kemitorvet
% 2800 Kgs. Lyngby, Denmark
% urbw@aqua.dtu.dk
%


idx=find(~cellfun('isempty',strfind(DS.filelist,string)))';
end