%% load data
% Fix bugs in V2 to make sure the Otsu mask is calculated from eTE0
% with minimal motion in whole image


clear

load('E:\BabyMotion\Minimal_To_Mask & Hct_by_sex\BIOCARD\Ctrl_Label_Diff_Ref_all.mat');
load('E:\BabyMotion\BIOCARD\MergedData\ROI_retained_Timepoint1.mat');
load('E:\BabyMotion\Minimal_To_Mask & Hct_by_sex\BIOCARD\hct_by_sex.mat');
load('E:\BabyMotion\BIOCARD\MergedData\Manual_exclusion_Timepoint1.mat');

if ~exist('scan_discard_flag')
    scan_discard_flag = zeros(size(control_all,4),1);
end

if ~exist('manual_exclusion')
    manual_exclusion = zeros(12,size(control_all,4));
end

if ~exist('hct_all')
    hct_all = zeros(size(control_all,4),1);
    hct_all(:) = 0.41;
end

fprintf('Use same mask for images in a scan to get L1\n');
%% parameters
rep_num = 3; eTE_num = 4;
ete_seq = [1 40 80 160];

%% reorganize retained data
retained_scan_num = length(find(scan_discard_flag==0));
control_all_retained = zeros(64,64,12,retained_scan_num);
diff_all_retained = zeros(64,64,12,retained_scan_num);
label_all_retained = zeros(64,64,12,retained_scan_num);
manual_exclusion_retained = zeros(12,retained_scan_num);
scan_discard_flag_retained = zeros(retained_scan_num,1);
retained_num = length(scan_discard_flag)-sum(scan_discard_flag);
ii=1;

for i = 1:length(scan_discard_flag)
    if ~scan_discard_flag(i)
        control_all_retained(:,:,:,ii) = control_all(:,:,:,i);
        diff_all_retained(:,:,:,ii) = diff_all(:,:,:,i);
        label_all_retained(:,:,:,ii) = label_all(:,:,:,i);
        manual_exclusion_retained(:,ii) = manual_exclusion(:,i);
        scan_discard_flag_retained(:,:,:,ii) = scan_discard_flag(i);
        ii = ii+1;
    end
end

%% get rid of SSS signal
dilateSize = 2;
SNR_selected = 15*((220/160)^2); % BIOCARD
fprintf('Using %d as SNR\n',SNR_selected);

[L1_bottom_noSSS_retained, mask_bottom_retained_noSSS] = CalL1NoSSS_threshV2(label_all_retained,diff_all_retained,dilateSize,SNR_selected);

L1_bottom_noSSS_Meanlabel_noSSS_retained = zeros(size(L1_bottom_noSSS_retained));
label_mean_noSSS_bottom = zeros(12,retained_scan_num);
average_label_mean_noSSS_bottom = zeros(12,retained_scan_num);

for irun = 1:retained_scan_num
    for idyn = 1:12
        cnt_label = label_all_retained(:,:,idyn,irun);
        label_mean_noSSS_bottom(idyn,irun) = mean(cnt_label(logical(mask_bottom_retained_noSSS(:,:,irun))));
    end
    
    for ieTE = 1:4
        cnt_run_label = label_all_retained(:,:,:,irun);
        range = (3*(ieTE-1)+1):3*ieTE;
        average_label_mean_noSSS_bottom(ieTE,irun) = mean(label_mean_noSSS_bottom(range,irun));
    end
end

for i = 1:retained_scan_num
    L1_bottom_noSSS_Meanlabel_noSSS_retained(:,i) = L1_bottom_noSSS_retained(:,i)./average_label_mean_noSSS_bottom(1,i);
end


%% find optimal theshold for median norm and min norm ratio L1
Range = [0.015 0.04];
step_size = 0.001;
nlist = 4;

% [opt_threshold,opt_T2,opt_Yv,opt_deltar2,discard_index] = OptimalDice_adult(temp_all,L1_bottom_noSSS_Meanlabel_noSSS_retained,diff_all_retained,control_all_retained,hct_all,Range,step_size,manual_exclusion_retained,nlist);
[opt_threshold,opt_T2,opt_Yv,opt_deltar2,discard_index,TurePositive,FalsePositive,AUC_value] = OptimalROC(temp_all,L1_bottom_noSSS_Meanlabel_noSSS_retained,diff_all_retained,control_all_retained,hct_all,manual_exclusion_retained,nlist);

opt_ind = int16((opt_threshold - Range(1))/step_size + 1);
opt_discard_result = cell(size(control_all_retained,4),1);
for i = 1:size(control_all_retained,4)
    opt_discard_result{i} = discard_index{opt_ind,i};
end

%% test on specific threshold

selected_threshold = 0.0099;
nlist = 4;

test_result = L1_bottom_noSSS_Meanlabel_noSSS_retained>selected_threshold;
% test_result = zeros(12,223);
[T2,Yv,ci,deltar2] = CalKeyIndV2_adult(diff_all,control_all,test_result,temp_all,hct_all,nlist);

discard_result = cell(size(control_all_retained,4),1);
for i = 1:size(control_all_retained,4)
    discard_result{i} = find(test_result(:,i));
end

%% Calculate Youden Index
Num_excluded = sum(manual_exclusion_retained(:));
TotalImgNum = length(manual_exclusion_retained(:));
TruePositiveMatrix = manual_exclusion_retained&test_result;
FalsePositiveMatrix = (~manual_exclusion_retained)&test_result;
TurePositive = sum(TruePositiveMatrix(:))/Num_excluded; 
FalsePositive = sum(FalsePositiveMatrix(:))/(TotalImgNum-Num_excluded);

YoudenInd = TurePositive-FalsePositive;
%% check images one by one
for i = 1:size(diff_all,4)
    show_imgs_sc(control_all(:,:,:,i),1:12,0,1000);title([num2str(i) 'th scan'])
    pause(3);
    close all
end