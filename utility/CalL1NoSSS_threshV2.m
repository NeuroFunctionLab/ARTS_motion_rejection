function [L1_retained, mask_noSSS_noise] = CalL1NoSSS_threshV2(label_all,diff_all,dilateSize,SNR_selected,dynnum)
    % Fix bugs to make sure brain mask is calculated from minimal motion
    % eTE0 image
    rownum = size(diff_all,1);
    columnnum = size(diff_all,2);
    Total_dyn_num = size(diff_all,3);

    L1_retained = zeros(Total_dyn_num,1);
    mask_noSSS_noise = zeros(rownum,columnnum);
    
    diff_all = abs(diff_all(:,:,:));

    L1_eTE0 = zeros(dynnum,1);
    for idyn = 1:dynnum
        SelectedDif = diff_all(:,:,idyn);
        L1_eTE0(idyn) = mean(SelectedDif(:));
    end
    minL1Ind = find(L1_eTE0==min(L1_eTE0));
    anat=diff_all(:,:,minL1Ind);
    selectedLabel = label_all(:,:,minL1Ind);
    % calculate mask
    NormalizedLabel = (selectedLabel-min(selectedLabel(:)))./(max(selectedLabel(:)-min(selectedLabel(:))));
    OtsuThreshold = graythresh(NormalizedLabel);
    OtsuWholeMask = NormalizedLabel>OtsuThreshold;
    OtsuWholeMask = imfill(OtsuWholeMask, 'holes');
    bottomOtsuMask = OtsuWholeMask;
    [rowInd, ~] = find(bottomOtsuMask);
    centroid_row = round(mean(rowInd));
    bottomOtsuMask(1:centroid_row, :) = 0;
    cnt_mask = bottomOtsuMask;

    % SNR test
    im_bottom = anat(cnt_mask);
    [values, ~] = sort(im_bottom(:), 'descend');
    mean_noise = mean(values(ceil(length(values)/2):end));
    theshold_selected_noise = mean_noise*SNR_selected;
    SSSmask_noise = cnt_mask.*(anat > theshold_selected_noise);
    Dilated_SSSmask_noise = imdilate(SSSmask_noise, strel('disk',dilateSize));
    mask_noSSS_noise(:,:) = cnt_mask.*(~Dilated_SSSmask_noise);

    for idyn = 1:Total_dyn_num
        cnt_dif_dyn = diff_all(:,:,idyn);
        L1_retained(idyn) = mean(cnt_dif_dyn(logical(mask_noSSS_noise(:,:))));
    end
end