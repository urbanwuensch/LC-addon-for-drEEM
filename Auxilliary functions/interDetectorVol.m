function [ DVout ] = interDetectorVol( DS,flowrate,idxPDA,idxFL,idxRI,PDAwave)
% [ DV ] = interDetectorVol( DS.flowrate,idxPDA,idxFL,idxRI,PDAwave)
% Notice:
% This mfile is part of the LC-addon for drEEM.
% % Check https://github.com/urbanwuensch/LC-addon-for-drEEM for the latest version
%
%
%
% interDectectorVol.m: Copyright (C) 2017 Urban J Wünsch
% Technical University of Denmark
% National Institute of Aquatic Resources
% Section for Oceans and Arctic
% Kemitorvet
% 2800 Kgs. Lyngby, Denmark
% urbw@aqua.dtu.dk
%
figure('InvertHardcopy','off','Color',[1 1 1],'units','normalized','outerposition',[0.25 0.4 0.5 0.5]);

for i=1:DS.nSample
    %%Fetch data for sample
    A.data=FeatureScaling(DS.ChData(idxFL).data(i,:));
    A.Rtime=DS.ChData(idxFL).Rtime;
   
    B.data=FeatureScaling(DS.ChData(idxRI).data(i,:));
    B.Rtime=DS.ChData(idxRI).Rtime;

    PDA.data=FeatureScaling(DS.ChData(idxPDA).data(i,:,knnsearch(DS.ChData(idxPDA).wave',PDAwave)));
    PDA.Rtime=DS.ChData(idxPDA).Rtime;
    
    PDA.data=interp1(PDA.Rtime',PDA.data',A.Rtime')';
    PDA.Rtime=A.Rtime;

    lwr=min(find(PDA.data>0.05));upr=max(find(B.data>0.05));
    
    A.Rtime=A.Rtime(lwr:upr);
    A.data=A.data(lwr:upr);
    B.Rtime=B.Rtime(lwr:upr);
    PDA.Rtime=PDA.Rtime(lwr:upr);
    B.data=B.data(lwr:upr);
    PDA.data=PDA.data(lwr:upr);
    
    % Modify data

    A.data(isnan(A.data)) = 0 ;B.data(isnan(B.data)) = 0 ;PDA.data(isnan(PDA.data)) = 0 ;
    f1 = fit(A.Rtime, A.data', 'smoothingspline','SmoothingParam', 0.99999);
    f2 = fit(B.Rtime, B.data', 'smoothingspline','SmoothingParam', 0.99999);
    f3 = fit(PDA.Rtime, PDA.data, 'smoothingspline','SmoothingParam', 0.99999);
    subplot(2,3,1);plot(f1,A.Rtime,A.data);legend off;hold on;title('FL');xlabel('Rtime');ylabel('Normalized intensity');ylim([0 1.1])
    subplot(2,3,2);plot(f2,B.Rtime,B.data);legend off;hold on;title('RI');xlabel('Rtime');ylim([0 1.1]);ylabel(' ')
    subplot(2,3,3);plot(f3,PDA.Rtime,PDA.data);legend off;hold on;title('PDA');xlabel('Rtime');ylim([0 1.1]);ylabel(' ')
    A.data=feval(f1,A.Rtime);
    B.data=feval(f2,A.Rtime);
    PDA.data=feval(f3,PDA.Rtime);

    
    % Calculate deadvolume
    Rvol=A.Rtime*flowrate;
    [~,A.x]=findpeaks(FeatureScaling(A.data),Rvol,'MinPeakHeight',0.8,'NPeaks',1,'SortStr','descend');
    [~,B.x]=findpeaks(FeatureScaling(B.data),Rvol,'MinPeakHeight',0.8,'NPeaks',1,'SortStr','descend');
    [~,PDA.x]=findpeaks(FeatureScaling(PDA.data),Rvol,'MinPeakHeight',0.8,'NPeaks',1,'SortStr','descend');

    disp(PDA.x/flowrate)
    disp(A.x/flowrate)
    disp(B.x/flowrate)
    
    DV(1,i)=A.x-PDA.x;
    DV(2,i)=B.x-PDA.x;

end

DVout(1,1)=nanmean(DV(1,:));
DVout(2,1)=nanmean(DV(2,:));
disp(' ')
disp(' ')
disp('---------------------')
disp(['COV for FL:',char(num2str(std(DV(1,:))/nanmean(DV(1,:))*100)),'%'])
disp(['COV for RI:',char(num2str(std(DV(2,:))/nanmean(DV(2,:))*100)),'%'])
disp('---------------------')
subplot(2,3,4:6);
plot(DS.ChData(idxFL).Rtime-DVout(1,1)/flowrate,FeatureScaling(DS.ChData(idxFL).data(1,:)));hold on
plot(DS.ChData(idxRI).Rtime-DVout(2,1)/flowrate,FeatureScaling(DS.ChData(idxRI).data(1,:)))
plot(DS.ChData(idxPDA).Rtime,FeatureScaling(DS.ChData(idxPDA).data(i,:,knnsearch(DS.ChData(idxPDA).wave',PDAwave))));title('Endresult example')
xlim([0 inf]);ylim([0 1.1]);xlabel('Retention time (aligned)');ylabel('Normalized intensity')
legend('FL','RI','PDA')

%% Tailing factors
% FL
disp(' ')
[~,c]=findpeaks(FeatureScaling(A.data),Rvol,'MinPeakHeight',0.8,'NPeaks',1,'SortStr','descend');
a=c-Rvol(knnsearch(FeatureScaling(A.data(Rvol<c)),0.1));
b=c-c+Rvol(knnsearch(FeatureScaling(A.data(Rvol>c)),0.1));
disp(['FL  tailing 10%: ',char(num2str(b/a))])

% RI
[~,c]=findpeaks(FeatureScaling(B.data),Rvol,'MinPeakHeight',0.8,'NPeaks',1,'SortStr','descend');
a=c-Rvol(knnsearch(FeatureScaling(B.data(Rvol<c)),0.1));
b=c-c+Rvol(knnsearch(FeatureScaling(B.data(Rvol>c)),0.1));
disp(['RI  tailing 10%: ',char(num2str(b/a))])


% PDA
% FL
[~,c]=findpeaks(FeatureScaling(PDA.data),Rvol,'MinPeakHeight',0.8,'NPeaks',1,'SortStr','descend');
a=c-Rvol(knnsearch(FeatureScaling(PDA.data(Rvol<c)),0.1));
b=c-c+Rvol(knnsearch(FeatureScaling(PDA.data(Rvol>c)),0.1));
disp(['PDA tailing 10%: ',char(num2str(b/a))])
disp('---------------------')
disp(' ')
disp(['PDA-FL volume: ',char(num2str(nanmean(DV(1,:))))]);
disp(['PDA-RI volume: ',char(num2str(nanmean(DV(2,:))))])
disp(' ')
end

