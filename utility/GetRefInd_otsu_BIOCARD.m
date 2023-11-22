function [RepTime4eTE, RefMask] = GetRefInd_otsu_BIOCARD(varargin)
    % Function to calculate best reference image among label images using different methods
    % based on maxmizing total similarity between reference images
    % Version 2.0: seq as input to specify acquisition sequence (eTE & dyn)

    I = varargin{1};
    seq = varargin{5};

    if strcmp(varargin{2}, 'whole')
        bApplyMask = 0;
    elseif strcmp(varargin{2}, 'mask')
        bApplyMask = 1;
        fprintf('Applying eTE0 mask to all eTE to find optimal reference images\n');
    else
        error('Input incorrect');
    end
    
    if strcmp(varargin{3}, 'MI')    
        fprintf('Using MI to caluclate similarity\n');
    elseif strcmp(varargin{3}, 'Normed MI')
        fprintf('Using normalized MI to caluclate similarity\n');
    elseif strcmp(varargin{3}, 'normxcorr2')
        fprintf('Using norm cross correlation to caluclate similarity\n');
    end
    
    if nargin < 4
        LCflag = 2; % Label: LCflag=1, Control: LCflag =2
    else
        LCflag = varargin{4};
    end
    
    eteseq = sort(unique(seq));
    seq1 = find(seq==eteseq(1)); seq2 = find(seq==eteseq(2));
    seq3 = find(seq==eteseq(3)); seq4 = find(seq==eteseq(4));

    TotalSimi = zeros(size(I,3),1);
    Comb = zeros(size(I,3),4);
    OtsuMask = logical(zeros([size(I(:,:,1)) 81]));
    i = 1;
    for ieTE1 = 1:3
        for ieTE2 = 1:3
            for ieTE3 = 1:3
                for ieTE4 = 1:3       

                    if LCflag == 1
                        Ind1 = seq1(ieTE1*2-1); 
                        Ind2 = seq2(ieTE2*2-1); 
                        Ind3 = seq3(ieTE3*2-1); 
                        Ind4 = seq4(ieTE4*2-1); 
                    elseif LCflag ==2
                        Ind1 = seq1(ieTE1*2); 
                        Ind2 = seq2(ieTE2*2); 
                        Ind3 = seq3(ieTE3*2); 
                        Ind4 = seq4(ieTE4*2); 
                    end

                    I1 = I(:,:,Ind1); I2 = I(:,:,Ind2);
                    I3 = I(:,:,Ind3); I4 = I(:,:,Ind4);

                    
                    grayLevel = graythresh(I1);
                    NormI1 = (I1-min(I1(:)))./(max(I1(:)-min(I1(:))));
                    OtsuMask(:,:,i) = NormI1>grayLevel;
                    OtsuMask(:,:,i) = imfill(OtsuMask(:,:,i), 'holes');
                    OtsuMask(:,:,i) = logical(OtsuMask(:,:,i));
                    if bApplyMask
                        if ~strcmp(varargin{3}, 'normxcorr2')
                            I1 = I1.*OtsuMask(:,:,i); I2 = I2.*OtsuMask(:,:,i); I3 = I3.*OtsuMask(:,:,i); I4 = I4.*OtsuMask(:,:,i);
                        end
                    end

                    if strcmp(varargin{3}, 'MI')
                        MI12 = mi(I1,I2);
                        MI13 = mi(I1,I3);
                        MI14 = mi(I1,I4);
                        TotalSimi(i) = MI12 + MI13 + MI14;
                    elseif strcmp(varargin{3}, 'Normed MI')
                        NormMI12 = MI2(I1,I2,'Normalized');
                        NormMI13 = MI2(I1,I3,'Normalized');
                        NormMI14 = MI2(I1,I4,'Normalized');
                        TotalSimi(i) = NormMI12 + NormMI13 + NormMI14;
                    elseif strcmp(varargin{3}, 'normxcorr2')
                        if bApplyMask
                            I1 = I1(OtsuMask(:,:,i)); I2 = I2(OtsuMask(:,:,i)); I3 = I3(OtsuMask(:,:,i)); I4 = I4(OtsuMask(:,:,i));
                            corr12 = normxcorr2(I1,I2);
                            corr13 = normxcorr2(I1,I3);
                            corr14 = normxcorr2(I1,I4);
                            len = length(I1(:));
                            TotalSimi(i) = corr12(len) + corr13(len) + corr14(len);
                        else
                            corr12 = normxcorr2(I1,I2);
                            corr13 = normxcorr2(I1,I3);
                            corr14 = normxcorr2(I1,I4);
                            TotalSimi(i) = corr12(64,64) + corr13(64,64) + corr14(64,64);
                        end
                    end

                    Comb(i,1) = ieTE1; Comb(i,2) = ieTE2; Comb(i,3) = ieTE3; Comb(i,4) = ieTE4;
                    i = i+1;
                end
            end
        end
    end

    MaxInd = find(TotalSimi == max(TotalSimi));
    RepTime4eTE = Comb(MaxInd,:);
    RefMask = OtsuMask(:,:,MaxInd);
end
