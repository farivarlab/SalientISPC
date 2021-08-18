%function quadrant_analysis()
dir = '/'; %set directory here

load([dir,'Subs_all_binneddata_movie1.mat']); %only 1-2 column of dim=2
load([dir,'Subs_all_binneddata_movie2.mat']); %only 3-4 column of dim=2

movie1_3D = squeeze(binneddata1(:,1,:,:));
movie1_2D = squeeze(binneddata1(:,2,:,:));
movie2_3D = squeeze(binneddata2(:,3,:,:));
movie2_2D = squeeze(binneddata2(:,4,:,:)); 

%Downsample mean distances down to 120/clip and concatenate clips
eyetrack1 = zeros(14,113,2); %clip1
eyetrack2 = zeros(14,113,2); %clip2
for i = 1:113
    t=movie1_3D(:,i*60-9:i*60+10,:);
    eyetrack1(:,i,:)=nanmean(t,2);
    t2=movie2_3D(:,i*60-9:i*60+10,:);
    eyetrack2(:,i,:)=nanmean(t2,2);
end
eyetrack_3D = cat(2,eyetrack1,eyetrack2);
eyetrack1 = zeros(14,113,2); %clip1
eyetrack2 = zeros(14,113,2); %clip2
for i = 1:113
    t=movie1_2D(:,i*60-9:i*60+10,:);
    eyetrack1(:,i,:)=nanmean(t,2);
    t2=movie2_2D(:,i*60-9:i*60+10,:);
    eyetrack2(:,i,:)=nanmean(t2,2);
end
eyetrack_2D = cat(2,eyetrack1,eyetrack2);
clear('eyetrack1','eyetrack2','movie1_3D','movie1_2D','movie2_3D','movie2_2D','binneddata1','binneddata2');

%load protocorr data
data_rh_2D=importdata([dir,'Movies_2D_rh.mat');
data_lh_2D=importdata([dir,'Movies_2D_lh.mat');
data_rh_3D=importdata([dir,'Movies_3D_rh.mat');
data_lh_3D=importdata([dir,'Movies_3D_lh.mat');
%add 7 TR offset by deleting first 7 TRs of data
data_lh_3D(:,[1:7,121:121+7-1])=[]; data_rh_3D(:,[1:7,121:121+7-1])=[];
data_lh_2D(:,[1:7,121:121+7-1])=[]; data_rh_2D(:,[1:7,121:121+7-1])=[];
%get data from V1-V3
cd([dir,'converted_to_std60/']);
roi_data_lh_3D = cell(1,6);roi_data_rh_3D = cell(1,6);
roi_data_lh_2D = cell(1,6);roi_data_rh_2D = cell(1,6);
% 01 - V1v	    
% 02 - V1d	   
% 03 - V2v	   
% 04 - V2d	   
% 05 - V3v	    
% 06 - V3d
for i=1:6
    templh = dlmread(['roi',num2str(i),'_lh.1D']);
    temprh = dlmread(['roi',num2str(i),'_rh.1D']);
    roi_data_lh_3D{1,i}=data_lh_3D(find(templh),:);
    roi_data_rh_3D{1,i}=data_rh_3D(find(temprh),:);
    roi_data_lh_2D{1,i}=data_lh_2D(find(templh),:);
    roi_data_rh_2D{1,i}=data_rh_2D(find(temprh),:);
end
clear('data_lh_3D','data_lh_2D','data_rh_3D','data_rh_2D');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
quads3D=zeros(4,226); quads2D=zeros(4,226);
%num subs in each quadrant
for tr = 1:226
    temp = squeeze(eyetrack_2D(:,tr,:));
    for sub = 1:14
        if (temp(sub,1)<920 && temp(sub,2)>550)
            quads2D(1,tr)= quads2D(1,tr)+1;
        elseif (temp(sub,1)>920 && temp(sub,2)>550)
            quads2D(2,tr)= quads2D(2,tr)+1;
        elseif (temp(sub,1)>920 && temp(sub,2)<550)
            quads2D(3,tr)= quads2D(3,tr)+1;
        elseif (temp(sub,1)<920 && temp(sub,2)<550)
            quads2D(4,tr)= quads2D(4,tr)+1;
        end
    end
    temp = squeeze(eyetrack_3D(:,tr,:));
    for sub = 1:14
        if (temp(sub,1)<920 && temp(sub,2)>550)
            quads3D(1,tr)= quads3D(1,tr)+1;
        elseif (temp(sub,1)>920 && temp(sub,2)>550)
            quads3D(2,tr)= quads3D(2,tr)+1;
        elseif (temp(sub,1)>920 && temp(sub,2)<550)
            quads3D(3,tr)= quads3D(3,tr)+1;
        elseif (temp(sub,1)<920 && temp(sub,2)<550)
            quads3D(4,tr)= quads3D(4,tr)+1;
        end
    end
end
%smooth
for i = 1:4
    conv_quads3D(i,:) = smooth(quads3D(i,:),0.04,'moving');
    conv_quads2D(i,:) = smooth(quads2D(i,:),0.04,'moving');
end

%compare to ISPC
ventral = [1,3,5];
dorsal = [2,4,6];
results = zeros(3,4,2,2); %V1-V3,quadrant,3D-2D,[r,p]
permrvals = zeros(1000,3,4,2); %perms, visual area, quadrant, 3D/2D
roiperms = zeros(226,1000);
for k = 1:3
    %left ventral 3D / quadrant 2
    tic
    roitemp = roi_data_lh_3D{1,ventral(k)};
    roitemp = nanmean(roitemp,1);
    roitemp = smooth(roitemp,0.04,'moving');
    [r,p] = corr(roitemp,conv_quads3D(2,:)');
    results(k,2,1,:) = [r,p];
    toc
    %left dorsal 3D / quadrant 3
    roitemp = roi_data_lh_3D{1,dorsal(k)};
    roitemp = nanmean(roitemp,1);
    roitemp = smooth(roitemp,0.04,'moving');
    [r,p] = corr(roitemp,conv_quads3D(3,:)');
    results(k,3,1,:) = [r,p];
    
    %right ventral 3D / quadrant 1
    roitemp = roi_data_rh_3D{1,ventral(k)};
    roitemp = nanmean(roitemp,1);
    roitemp = smooth(roitemp,0.04,'moving');
    [r,p] = corr(roitemp,conv_quads3D(1,:)');
    results(k,1,1,:) = [r,p];
   
    %right dorsal 3D / quadrant 4
    roitemp = roi_data_rh_3D{1,dorsal(k)};
    roitemp = nanmean(roitemp,1);
    roitemp = smooth(roitemp,0.04,'moving');
    [r,p] = corr(roitemp,conv_quads3D(4,:)');
    results(k,4,1,:) = [r,p];
    
    %left ventral 2D / quadrant 2
    roitemp = roi_data_lh_2D{1,ventral(k)};
    roitemp = nanmean(roitemp,1);
    roitemp = smooth(roitemp,0.04,'moving');
    [r,p] = corr(roitemp,conv_quads2D(2,:)');
    results(k,2,2,:) = [r,p];
    
    %left dorsal 3D / quadrant 3
    roitemp = roi_data_lh_2D{1,dorsal(k)};
    roitemp = nanmean(roitemp,1);
    roitemp = smooth(roitemp,0.04,'moving');
    [r,p] = corr(roitemp,conv_quads2D(3,:)');
    results(k,3,2,:) = [r,p];
    
    %right ventral 3D / quadrant 1
    roitemp = roi_data_rh_2D{1,ventral(k)};
    roitemp = nanmean(roitemp,1);
    roitemp = smooth(roitemp,0.04,'moving');
    [r,p] = corr(roitemp,conv_quads2D(1,:)');
    results(k,1,2,:) = [r,p];
    
    %right dorsal 3D / quadrant 4
    roitemp = roi_data_rh_2D{1,dorsal(k)};
    roitemp = nanmean(roitemp,1);
    roitemp = smooth(roitemp,0.04,'moving');
    [r,p] = corr(roitemp,conv_quads2D(4,:)');
    results(k,4,2,:) = [r,p];
end
presults3D = squeeze(results(:,:,1,2));
presults2D = squeeze(results(:,:,2,2));
rresults3D = squeeze(results(:,:,1,1));
rresults2D = squeeze(results(:,:,2,1));
test = reshape(presults3D, [12,1]); test = mafdr(test,'BHFDR','true'); presults3D = reshape(test,[3,4]);
test = reshape(presults2D, [12,1]); test = mafdr(test,'BHFDR','true'); presults2D = reshape(test,[3,4]);





