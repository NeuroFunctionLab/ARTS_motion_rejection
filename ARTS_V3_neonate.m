%% load data
% Fix bugs in V2 to make sure the Otsu mask is calculated from eTE0
% with minimal motion in whole image

clear
folderName = 'E:\BabyMotion\Control_baby_data_20221215';
load([folderName '\Ctrl_Label_Diff_Ref_all.mat']);
load([folderName '\Manual_Exclusion_all.mat']);
load([folderName '\Motion_vec_all.mat']);

% folderName = 'E:\BabyMotion\HIE_Baby_20230118';
% load([folderName '\Ctrl_Label_Diff_Ref_all.mat']);
% load([folderName '\Manual_Exclusion.mat']);


if ~exist('scan_discard_flag')
    scan_discard_flag = zeros(size(control_all,4),1);
end

if ~exist('manual_exclusion')
    manual_exclusion = zeros(12,size(control_all,4));
end

bBottomAlready = 0;
fprintf('Use same mask for images in a scan to get L1\n');
%% parameters
rep_num = 3; eTE_num = 4;
ete_seq = [0.44 40 80 160];
retained_scan_num = length(find(scan_discard_flag==0));

%% reorganize retained data
control_all_retained = zeros(64,64,12,retained_scan_num);
diff_all_retained = zeros(64,64,12,retained_scan_num);
label_all_retained = zeros(64,64,12,retained_scan_num);
manual_exclusion_retained = zeros(12,retained_scan_num);
scan_discard_flag_retained = zeros(retained_scan_num,1);
% motion_vec_retained = cell(retained_scan_num,4);
retained_num = length(scan_discard_flag)-sum(scan_discard_flag);
ref_retained = zeros(4,retained_scan_num);

ii=0;
for i = 1:length(scan_discard_flag)
    if ~scan_discard_flag(i)
        ii = ii+1;
        control_all_retained(:,:,:,ii) = control_all(:,:,:,i);
        diff_all_retained(:,:,:,ii) = diff_all(:,:,:,i);
        label_all_retained(:,:,:,ii) = label_all(:,:,:,i);
        manual_exclusion_retained(:,ii) = manual_exclusion(:,i);
        ref_retained(:,ii) = ref_all(:,i);
        scan_discard_flag_retained(:,:,:,ii) = scan_discard_flag(i);
%         for ieTE = 1:eTE_num
%             motion_vec_retained{ii,ieTE} = zeros(6,6);
%             motion_vec_retained{ii,ieTE} = motion_vec{i,ieTE};
%         end
    end
end

%% get rid of SSS signal
load([folderName '\ROI_retained.mat']);
load([folderName '\Ctrl_Label_Diff_Ref_retained.mat']);
load([folderName '\Manual_Exclusion_retained.mat']);
dilateSize = 2;
SNR_selected = 15; % baby
fprintf('Using %d as SNR\n',SNR_selected);

[L1_bottom_noSSS_retained, mask_bottom_retained_noSSS] = CalL1NoSSS_threshV2(label_all_retained,diff_all_retained,dilateSize,SNR_selected);


L1_bottom_noSSS_Meanlabel_noSSS_retained = zeros(size(L1_bottom_noSSS_retained));
L1_bottom_noSSS_MedianTissue_noSSS_retained = zeros(size(L1_bottom_noSSS_retained));
label_mean_noSSS_bottom = zeros(12,retained_scan_num);
average_label_mean_noSSS_bottom = zeros(eTE_num,retained_scan_num);

for irun = 1:retained_scan_num
    for idyn = 1:12
        cnt_label = label_all_retained(:,:,idyn,irun);
        label_mean_noSSS_bottom(idyn,irun) = mean(cnt_label(logical(mask_bottom_retained_noSSS(:,:,irun))));
    end
    
    for ieTE = 1:4
        range = (3*(ieTE-1)+1):3*ieTE;
        average_label_mean_noSSS_bottom(ieTE,irun) = mean(label_mean_noSSS_bottom(range,irun));
    end
end

for i = 1:retained_scan_num
    L1_bottom_noSSS_Meanlabel_noSSS_retained(:,i) = L1_bottom_noSSS_retained(:,i)./average_label_mean_noSSS_bottom(1,i);
    L1_bottom_noSSS_MedianTissue_noSSS_retained(:,i) = L1_bottom_noSSS_retained(:,i)./median(L1_bottom_noSSS_retained(:,i));
end

%% find optimal theshold for median norm and min norm ratio L1
Range = [0.028 0.07];
step_size = 0.001;

% Range = [1 1.5];
% step_size = 0.01;

nlist = 4;
load([folderName '\hct_retained_by_sex_average.mat']);
% OptimizeMethod = 'Mean';

% [opt_threshold,opt_T2,opt_Yv,opt_deltar2,discard_index] = OptimalDice(temp_all,L1_bottom_noSSS_Meanlabel_noSSS_retained,diff_all_retained,control_all_retained,hct_all,Range,step_size,manual_exclusion_retained,nlist);
% [opt_threshold,opt_T2,opt_Yv,opt_deltar2,discard_index] = OptimalDeltaR2(temp_all,L1_bottom_noSSS_Meanlabel_noSSS_retained,diff_all_retained,control_all_retained,hct_all,Range,step_size,nlist,OptimizeMethod);
[opt_threshold,opt_T2,opt_Yv,opt_deltar2,discard_index,TurePositive,FalsePositive,AUC_value] = OptimalROC(temp_all,L1_bottom_noSSS_Meanlabel_noSSS_retained,diff_all_retained,control_all_retained,hct_all,manual_exclusion_retained,nlist);

% [opt_threshold,opt_T2,opt_Yv,opt_deltar2,discard_index] = OptimalDice(temp_all,L1_bottom_noSSS_MedianTissue_noSSS_retained,diff_all_retained,control_all_retained,hct_all,Range,step_size,manual_exclusion_retained,nlist);
% [opt_threshold,opt_T2,opt_Yv,opt_deltar2,discard_index] = OptimalDeltaR2(temp_all,L1_bottom_noSSS_MedianTissue_noSSS_retained,diff_all_retained,control_all_retained,hct_all,Range,step_size,nlist);


opt_ind = int16((opt_threshold - Range(1))/step_size + 1);
opt_discard_result = cell(size(control_all_retained,4),1);
for i = 1:size(control_all_retained,4)
    opt_discard_result{i} = discard_index{opt_ind,i};
end
test_result = L1_bottom_noSSS_Meanlabel_noSSS_retained>opt_threshold(1);
dice_result = dice(test_result,logical(manual_exclusion_retained))

%% test on specific threshold
% load('E:\BabyMotion\Control_baby_data_20220913\Previous_Manual_Result\hct_all.mat');
load([folderName '\hct_retained_by_sex_average.mat']);


selected_threshold = 0.0353;
% selected_threshold = 1.33;
nlist = 4;

test_result = L1_bottom_noSSS_Meanlabel_noSSS_retained>selected_threshold;
% test_result = L1_bottom_noSSS_MedianTissue_noSSS_retained>selected_threshold;
[T2,Yv,ci,deltar2] = CalKeyIndV2(diff_all_retained,control_all_retained,test_result,temp_all,hct_all,nlist);

discard_result = cell(size(control_all_retained,4),1);
for i = 1:size(control_all_retained,4)
    discard_result{i} = find(test_result(:,i));
end
dice_result = dice(test_result,logical(manual_exclusion_retained))

%% Test adaptive threshold
load([folderName '\hct_retained_by_sex.mat']);
nlist = 4;
range = [0.03 0.07];
step_size = 0.001;

[T2,Yv,thresh,deltar2,test_result] = AdaptiveThreshold(diff_all_retained,control_all_retained,L1_bottom_noSSS_Meanlabel_noSSS_retained,range,step_size,temp_all,hct_all,nlist);

discard_result = cell(size(control_all_retained,4),1);
for i = 1:size(control_all_retained,4)
    discard_result{i} = find(test_result(:,i));
end
dice_result = dice(test_result,manual_exclusion_retained)