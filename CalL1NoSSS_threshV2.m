function [L1_retained, mask_noSSS_noise] = CalL1NoSSS_threshV2(label_all,diff_all,dilateSize,SNR_selected)
    % Fix bugs to make sure brain mask is calculated from minimal motion
    % eTE0 image
    rownum = size(diff_all,1);
    columnnum = size(diff_all,2);
    dyn_num = size(diff_all,3);
    run_num = size(diff_all,4);

    L1_retained = zeros(dyn_num,run_num);
    mask_noSSS_noise = zeros(rownum,columnnum,run_num);
    
    for irun = 1:run_num
%         cnt_mask = logical(mask_all(:,:,irun));
        cnt_dif = abs(diff_all(:,:,:,irun));
        cnt_label = label_all(:,:,:,irun);
        
        L1_eTE0 = zeros(3,1);
        for idyn = 1:3
            SelectedDif = cnt_dif(:,:,idyn);
            L1_eTE0(idyn) = mean(SelectedDif(:));
        end
        minL1Ind = find(L1_eTE0==min(L1_eTE0));
        anat=cnt_dif(:,:,minL1Ind);
        selectedLabel = cnt_label(:,:,minL1Ind);
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
        mask_noSSS_noise(:,:,irun) = cnt_mask.*(~Dilated_SSSmask_noise);

%         if irun == 24   
%             figure();
%             subplot(2,2,1);imshow(anat,[0 100]);title('diff image')
%             subplot(2,2,2);imshow(SSSmask_noise);title('SSS mask')
%             subplot(2,2,3);imshow(Dilated_SSSmask_noise);title('Dilated SSS mask')
%             subplot(2,2,4);imshow(mask_noSSS_noise(:,:,irun));title('mask without SSS')
%         end


        for idyn = 1:dyn_num
            cnt_dif_dyn = cnt_dif(:,:,idyn);
            L1_retained(idyn,irun) = mean(cnt_dif_dyn(logical(mask_noSSS_noise(:,:,irun))));
        end
    end
end