function [T2,Yv,ci,deltar2] = CalKeyIndV2(dif,con,test_result,temp_all,hct_all,nlist)
    % Calculate key indicator with input size as [row col dyn scan]
    
    DataNum = length(hct_all);
    T2 = zeros(DataNum,1);
    Yv = zeros(DataNum,1);
    deltar2 = zeros(DataNum,1);
    % temp_all = false(64,64,27);
    dynnum = 6;
    etelen = 4;
    rownum = 64;
    columnnum = 64;
    roiflag = 0;
    seq = [0.440000000000000,0.440000000000000,40,40,80,80,160,160,0.440000000000000,0.440000000000000,40,40,80,80,160,160,0.440000000000000,0.440000000000000,40,40,80,80,160,160];
    te = [0.440000000000000;40;80;160];
    ti = 1020;
    bloodt1 = 1624;
    dyn = [6 6 6 6];
    
    for ii = 1:DataNum
        cnt_dif = dif(:,:,:,ii);
        cnt_con = con(:,:,:,ii);
        for i = 1:etelen
            dif_all{i} = cnt_dif(:,:,(3*(i-1)+1):3*i);
            con_all{i} = cnt_con(:,:,(3*(i-1)+1):3*i);
        end
        ind2exclude = find(test_result((12*(ii-1)+1):12*ii));
        % check if any excluded dynamics are eTE0
        ind2include_eTE0 = [1:dynnum/2];
        ind2exclude_eTE0 = ind2exclude(ind2exclude <= dynnum/2);
        ind2include_eTE0(ind2exclude_eTE0) = [];
        
        if length(ind2exclude) >= size(test_result,1)-1 
            % if less than 2 images are retained, set all results to 0 to
            % avoid error report
            ci=0;
            Yv(ii)=0;
            T2(ii) = 0;
            deltar2(ii) = 0;
        else
            % draw ROI --------------------------------------------------------
            for t=1:etelen
                % draw ROI
                anat=mean(cnt_dif(:,:,ind2include_eTE0), 3);
                droi=zeros(rownum, columnnum);
    %             figure;
        %         himage=imshow(anat,[min(anat(:)) max(anat(:))]);
                temp = temp_all(:,:,ii);
        %          temp=roipoly();colorbar;
        %          temp_all(:,:,ii) = temp;
        %         temp = false(64,64);
        %         temp(47:59,25:41)=1;
                droi=droi | temp;
                close;
                clear mylist;
                mycount=0;        
                for i=1:rownum
                    for j=1:columnnum
                        if droi(i,j)==1
                            mycount=mycount+1;
                            mylist(mycount,1)=i;
                            mylist(mycount,2)=j;
                            mylist(mycount,3)=anat(i,j,1);
                        end
                    end
                end
                [Y,I]=sort(mylist(:,3),'descend');
                if roiflag ==0
                    mask=zeros(rownum,columnnum);
                    for tt=1:nlist
                        mask(mylist(I(tt),1),mylist(I(tt),2))=1;
                    end
                    break;
                else
                    for tt=1:nlist
                        mask(mylist(I(tt),1),mylist(I(tt),2),t)=1;
                    end
                end
            end
        %     display_fmri(anat,mask);
        
            % fitting for T2 ----------------------------------------
            nte=etelen;
            for i=1:nte
                clear temp3 temp5;
                temp3=dif_all{i};
                temp5=con_all{i};
                if roiflag ==0
                    temp4=average_actpixel1(1,mask,temp3);
                    b_all{i}=temp4(:);
                    temp6=average_actpixel1(1,mask,temp5);
                    c_all{i}=temp6(:);
                else
                    temp4=average_actpixel1(1,mask(:,:,i),temp3);
                    b_all{i}=temp4(:);
                    temp6=average_actpixel1(1,mask(:,:,i),temp5);
                    c_all{i}=temp6(:);
                end   
            end
        
            %%% T2 fitting by all b points
            tetemp=seq(1:2:end);
            te_all=sort(tetemp,'ascend');
            b_all_vector = [];
            c_all_vector = [];
            for i=1:nte
                b_all_vector = [b_all_vector; b_all{i}];
                c_all_vector = [c_all_vector; c_all{i}];
            end
        
            % Dengrong Jiang, 5/14/2018: Compute labelling efficiency
            Ind4LabelEff = find(te_all == te(1));
            LabelEff = mean(b_all_vector(Ind4LabelEff))/mean(c_all_vector(Ind4LabelEff));
        
            dyn0=0;
            for i=1:etelen-1
                b_all_vector(dyn0+1:dyn0+dyn(i)/2)=exp(-(ti-te(end))/bloodt1)/exp(-(ti-te(i))/bloodt1)*b_all_vector(dyn0+1:dyn0+dyn(i)/2);
                c_all_vector(dyn0+1:dyn0+dyn(i)/2)=exp(-(ti-te(end))/bloodt1)/exp(-(ti-te(i))/bloodt1)*c_all_vector(dyn0+1:dyn0+dyn(i)/2);
                dyn0=dyn0+dyn(i)/2;
            end
        
            % Exclude dynamics with strong motion artifacts
            ind2fit = [1:length(te_all)];
            ind2fit(ind2exclude) = [];
        
            [temp1,resid, jacob]=nlinfit_hlu(te_all(ind2fit)/1000, b_all_vector(ind2fit),'monexp_model',[100,10]);
            t2=1000./temp1(2);
        
            conintval=nlparci(temp1,resid, jacob); %95% confidence interval for estimates
            aa=1000./conintval;
            ci=aa(2,:);
        
            [temp1_c,resid_c, jacob_c]=nlinfit_hlu(te_all(ind2fit)/1000, c_all_vector(ind2fit),'monexp_model',[100,10]);
            conintval_c=nlparci(temp1_c,resid_c, jacob_c); %95% confidence interval for estimates
            t2_c=1000./temp1_c(2);
            aa_c=1000./conintval_c;
            ci_c=aa_c(2,:);
            Yv(ii)=100*neoT2toY(t2,hct_all(ii));
            T2(ii) = t2;
            deltar2(ii) = 1000/ci(2)-1000/ci(1);
        
    %         if ii == 54
    %             exclusion_ind = [3 6 8 9 11 12];
    % %             exclusion_ind2 = [6];
    % 
    %             dot_color = [56, 92, 163]/255;
    %             dot_line_width = 1.5;
    %             fontsize_display = 25;
    %             figure, plot(te_all, b_all_vector, 'ok', 'linewidth', dot_line_width, 'markerfacecolor', dot_color, 'markeredgecolor', dot_color, 'markersize', 7);
    %             hold on; plot(te_all(exclusion_ind),b_all_vector(exclusion_ind),'ok', 'linewidth', dot_line_width, 'markerfacecolor', 'red', 'markeredgecolor', 'red', 'markersize', 7);
    % %             hold on; plot(te_all(exclusion_ind2),b_all_vector(exclusion_ind2),'ok', 'linewidth', dot_line_width, 'markerfacecolor', "#77AC30", 'markeredgecolor', "#77AC30", 'markersize', 7);
    %             set(gca, 'fontsize', fontsize_display, 'fontweight', 'bold', 'linewidth', 2, 'tickdir', 'out', 'ycolor', [0,0,0], 'xcolor', [0,0,0]);
    %             xlim([-10 180]);
    %             box off
    %         end
        end
    end
end