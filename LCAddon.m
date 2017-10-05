% The LC-addon for drEEM for MATLAB (ver. 0.1)
% Version 0.1 October-2017
%
% Copyright (C) 2017  Urban J. Wuensch
% Technical University of Denmark, Kemitorvet, 2800 Kgs. Lyngby
% urbw@aqua.dtu.dk
%%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation version 3 of the License <http://www.gnu.org/licenses/>
%
% This program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with 
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin 
% Street, Fifth Floor, Boston, MA  02110-1301, USA.
%
%-----------------------------------------------
%DATA IMPORT/EXPORT
%
%importrun                 Import ASCII data from Shimadzu LabSolutions
%hplc2eem                  Transfer HPSEC data scheme into drEEM format
%
%-----------------------------------------------
%DATA MANIPULATION
%
%baselcor                  Correct baseline drifts by subtraction of baseline function
%mergeChannels             Merge Channels stored in separate variables belonging to one sample
%Rtime2Evol                Convert retention time to elution volume with the flow diagram
%subsetHPLC                Remove samples and / or retention times from datasets
%interplXaxis              Interpolate retention time or elution volume to fit given values & Eliminate elution volume and retention time difference between detectors

%-----------------------------------------------
%PLOTTING
%
%view3D                    view chromatograms of 3D channels
%viewHPLC                  view chromatograms of 2D channels
%
%-----------------------------------------------
%AUXILLARY FUNCTIONS
%
%checkIntegrity            Confirm the integrity of a dataset. Matching dimensions etc.
%chooseXaxis               Finds x-axis present in datasets (retention time or elution volume) and returns values
%InterDetectorVol          Determines interdetector volume
%nameFind                  Returns index of string. Useful when trying to plot specific sample.