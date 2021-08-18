%Correlating the mean distances in eye tracking data to ISPC

%Eyetracking data begins at movie
%Movie begins in fMRI data at 10s in (5 TRs)

method=2; %method of downsampling eyetracking data: 1=interpolation 2=binning
offset=7; %HRF lag offset
stereo = 1; %1 = 3D, 2 = 2D
dir = '/'; %set directory here

if stereo==1
    data_rh=importdata([dir,'Movies_3D_rh.mat']);
    data_lh=importdata([dir,'Movies_3D_lh.mat']);
else
    data_rh=importdata([dir,'Movies_2D_rh.mat']);
    data_lh=importdata([dir,'Movies_2D_lh.mat']);
end

data_lh(:,[1:offset,121:121+offset-1])=[]; data_rh(:,[1:offset,121:121+offset-1])=[];

eyetracking_data1= importdata([dir,'meandistances_movie12D.mat']);  
eyetracking_data2= importdata([dir,'meandistances_movie22D.mat']);

%Downsample mean distances down to 120
if method==1
    eyetrack1 = interp1(eyetracking_data1,1:60:6780,'PCHIP');
    eyetrack2 = interp1(eyetracking_data2,1:60:6780,'PCHIP');
else
    eyetrack1 = zeros(1,113);
    eyetrack2 = zeros(1,113);
    for i = 1:113
        t=eyetracking_data1(1,i*60-9:i*60+10);
        eyetrack1(1,i)=mean(t);
        t2=eyetracking_data2(1,i*60-9:i*60+10);
        eyetrack2(1,i)=mean(t2);
    end
end
clear('t','t2','i','j','eyetracking_data1','eyetracking_data2','data_lh','data_rh')

eyetrack = cat(2,eyetrack1,eyetrack2);
%flip distances by making all numbers negative then adding 2*mean
m1 = max(eyetrack)+min(eyetrack);
eyetrack = -(eyetrack);
eyetrack = eyetrack + m1;
%make it small
eyetrack = eyetrack/10000; %100000
clear('eyetrack1','eyetrack2');

%% per visual area correlation
cd([dir,'\consolidated_areas_std60']);
pvals = zeros(22,3);
rvals = zeros(22,3);

%left
va_results_lh=zeros(22,2);
for i=1:22
    VA_ROIs = dlmread(['roi',num2str(i),'_lh.1D']);
    temp = data_lh(find(any(VA_ROIs,2)),:);
    temp = squeeze(nanmean(temp,1));
    temp = smooth(temp,0.04,'moving');
    eye = smooth(eyetrack,0.04,'moving');
    [r,p,rlo,rup]=corrcoef(eye,temp);
    va_results_lh(i,1)=r(1,2);
    va_results_lh(i,2)=p(1,2);
end
pvals(:,1) = mafdr(va_results_lh(:,2),'BHFDR','true')
rvals(:,1) = va_results_lh(:,1);

%right
va_results_rh=zeros(22,2);
for i=1:22
    VA_ROIs = dlmread(['roi',num2str(i),'_rh.1D']);
    temp = data_rh(find(any(VA_ROIs,2)),:);
    temp = squeeze(nanmean(temp,1));
    temp = smooth(temp,0.04,'moving');
    eye = smooth(eyetrack,0.04,'moving');
    [r,p,rlo,rup]=corrcoef(eye,temp);
    va_results_rh(i,1)=r(1,2);
    va_results_rh(i,2)=p(1,2);
end
pvals(:,2) = mafdr(va_results_rh(:,2),'BHFDR','true')
rvals(:,2) = va_results_rh(:,1);

%both
va_results=zeros(22,4);
for i=1:22
    VA_ROI_lh = dlmread(['roi',num2str(i),'_lh.1D']);
    VA_ROI_rh = dlmread(['roi',num2str(i),'_rh.1D']);
    temp_lh = data_lh(find(any(VA_ROI_lh,2)),:);
    temp_rh = data_rh(find(any(VA_ROI_rh,2)),:);
    temp = [temp_lh',temp_rh'];
    temp = mean(temp,2)';
    eye = smooth(eyetrack,0.04,'moving');
    temp = smooth(temp,0.04,'moving')
    [r,p,rlo,rup]=corrcoef(eye,temp)
    va_results(i,1)=r(1,2);
    va_results(i,2)=p(1,2);
    va_results(i,3)=rlo(1,2);
    va_results(i,4)=rup(1,2);
end
pvals(:,3) = mafdr(va_results(:,2),'BHFDR','true')
rvals(:,3) = va_results(:,1);
ci = zeros(22,2);
ci(:,1) = va_results(:,3);
ci(:,2) = va_results(:,4);


save('pvals_movie24_eyeclip2D', 'pvals', '-v7.3');
save('Rvals_movie24_eyeclip2D', 'rvals', '-v7.3');
save('ci_movie24_eyeclip2D', 'ci', '-v7.3');


%%
%per node correlation
%LEFT
protoeye_corr = zeros(36002,2);
for i = 1:36002
   eye = smooth(eyetrack,0.04,'moving'); %0.08 or 0.04 moving
   temp = smooth(data_lh(i,:),0.04,'moving'); %0.05 rloess or 0.04 moving
   [R,P] = corrcoef(eye,temp);
   protoeye_corr(i,1) = R(1,2);
   protoeye_corr(i,2) = P(1,2);
   i
end
protoeye_corr_fdr = protoeye_corr(:,2)*36002;
results = zeros(36002,2);
results(protoeye_corr_fdr<0.05,1)= protoeye_corr(protoeye_corr_fdr<0.05,1);
results(protoeye_corr_fdr<0.05,2)= protoeye_corr(protoeye_corr_fdr<0.05,2);

dlmwrite('3D_rvals_lh.1D', results(:,1));
dlmwrite('3D_pvals_lh.1D', results(:,2));

%RIGHT
protoeye_corr = zeros(36002,2);
for i = 1:36002
   eye = smooth(eyetrack,0.04,'moving'); %0.08 or 0.04 moving
   temp = smooth(data_rh(i,:),0.04,'moving'); %0.05 rloess or 0.04 moving
   [R,P] = corrcoef(eye,temp);
   protoeye_corr(i,1) = R(1,2);
   protoeye_corr(i,2) = P(1,2);
   i
end
protoeye_corr_fdr = protoeye_corr(:,2)*36002;

results = zeros(36002,2);
results(protoeye_corr_fdr<0.05,1)= protoeye_corr(protoeye_corr_fdr<0.05,1);
results(protoeye_corr_fdr<0.05,2)= protoeye_corr(protoeye_corr_fdr<0.05,2);

dlmwrite('3D_rvals_rh.1D', results(:,1));
dlmwrite('3D_pval_rh.1D', results(:,2));