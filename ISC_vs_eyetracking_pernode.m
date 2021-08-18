%ISC correlation with saliency at each node across cortex

%Eyetracking data begins at movie
%Movie begins in fMRI data at 10s in (5 TRs)

method=2; %method of downsampling eyetracking data: 1=interpolation 2=binning
offset=7; %HRF lag offset
stereo = 1; %1 = 3D, 2 = 2D
dir = '/'; %set directory here

if stereo==1
    data_rh(:,1:120)=importdata([dir,'Movies_3D_1_rh.mat']);
    data_lh(:,1:120)=importdata([dir,'Movies_3D_1_lh.mat']);
    data_rh(:,121:240)=importdata([dir,'Movies_3D_2_rh.mat']);
    data_lh(:,121:240)=importdata([dir,'Movies_3D_2_lh.mat']);    
else
    data_rh(:,1:120)=importdata([dir,'Movies_2D_1_rh.mat']);
    data_lh(:,1:120)=importdata([dir,'Movies_2D_1_lh.mat']);
    data_rh(:,121:240)=importdata([dir,'Movies_2D_2_rh.mat']);
    data_lh(:,121:240)=importdata([dir,'Movies_2D_2_lh.mat']);    
end

%apply offset
data_lh(:,[1:offset,121:121+offset-1],:)=[]; data_rh(:,[1:offset,121:121+offset-1],:)=[];

if stereo == 1
    eyetracking_data1= importdata([dir,'\meandistances_movie13D.mat']);
    eyetracking_data2= importdata([dir,'\meandistances_movie23D.mat']);
else
    eyetracking_data1= importdata([dir,'\meandistances_movie12D.mat']);
    eyetracking_data2= importdata([dir,'\meandistances_movie22D.mat']);
end

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
clear('eyetrack1','eyetrack2');

%convolve with HRF
t = 1:1:30; % MEASUREMENTS
h = gampdf(t,6) + -.5*gampdf(t,16); % HRF MODEL
h = h/max(h); % SCALE HRF TO HAVE MAX AMPLITUDE OF 1
plot(h);
eyetrack2 = conv(h,eyetrack);
eye = eyetrack2;
eye(:,227:end)=[];
%% per node correlation
results_lh = zeros(36002,2);
results_rh = zeros(36002,2);
ISC_lh = zeros(1,226); ISC_rh = zeros(1,226);
for node = 2:36002
    temp_lh = squeeze(data_lh(node,:,:));
    temp_rh = squeeze(data_rh(node,:,:));
   %sliding window ISC
    for time = 5:222
        part_lh = temp_lh(time-4:time+4,:); part_rh = temp_lh(time-4:time+4,:);
        part_lh = atanh(corrcoef(part_lh)); part_rh = atanh(corrcoef(part_rh));
        part_lh(isinf(part_lh))=0; part_rh(isinf(part_rh))=0;
        part_lh = tanh(nanmean(nonzeros(triu(part_lh)))); part_rh = tanh(nanmean(nonzeros(triu(part_rh))));
        ISC_lh(1,time)=part_lh; ISC_rh(1,time)=part_rh;
    end
    ISC_lh(1,1:4) = [ISC_lh(1,5),ISC_lh(1,5),ISC_lh(1,5),ISC_lh(1,5)];
    ISC_rh(1,1:4) = [ISC_rh(1,5),ISC_rh(1,5),ISC_rh(1,5),ISC_rh(1,5)];
    ISC_lh(1,223:226) = [ISC_lh(1,222),ISC_lh(1,222),ISC_lh(1,222),ISC_lh(1,222)];
    ISC_rh(1,223:226) = [ISC_rh(1,222),ISC_rh(1,222),ISC_rh(1,222),ISC_rh(1,222)];

    [r,p]=corrcoef(eye,ISC_lh);
    results_lh(node,1)=r(1,2);
    results_lh(node,2)=p(1,2);
    [r,p]=corrcoef(eye,ISC_rh);
    results_rh(node,1)=r(1,2);
    results_rh(node,2)=p(1,2);
    
end
%bonferonni
pvals_lh = results_lh(:,2)*36002;
pvals_rh = results_rh(:,2)*36002;

finalresults_lh = zeros(36002,1);
finalresults_rh = zeros(36002,1);

finalresults_lh(pvals_lh<0.05,1)= results_lh(pvals_lh<0.05,1);
finalresults_rh(pvals_rh<0.05,1)= results_rh(pvals_rh<0.05,1);

if stereo == 1
    dlmwrite('3D_rvals_lh.1D', finalresults_lh);
    dlmwrite('3D_rvals_rh.1D', finalresults_rh);  
else
    dlmwrite('2D_rvals_lh.1D', finalresults_lh);
    dlmwrite('2D_rvals_rh.1D', finalresults_rh);
end