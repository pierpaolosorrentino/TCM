%% zscore =3
% load('G:\multiple sclerosis\avalanches\avalanches_MS');
% safety_aval=avalanches;
% avalanches=avalanches{1}{6};
% addpath('G:\Melbourne\markov_stuff\');

% nuovo
addpath('C:\Users\ppsor\Desktop\Code')
addpath(genpath('C:\Users\ppsor\Desktop\studi\TCM'));
clear all
% % % load('C:\Users\ppsor\Desktop\studi\TCM\avalanches\avalanches_allbins_Samelength_dwnspl2');
% load('C:\Users\ppsor\Desktop\studi\TCM\avalanches\avalanches_allbins_dspl2_84');
% load('C:\Users\ppsor\Desktop\studi\TCM\avalanches\avalanches_allbins_sameleng_dspl2_84');
% load('C:\Users\ppsor\Desktop\studi\TCM\data_MEG\ms_maggio_2020\avalanches_allbins_sameleng_dspl2_84_Maggio');
% %
% load('C:\Users\ppsor\Desktop\studi\TCM\avalanches\avalanches_allbins_dspl2_66');
load('C:\Users\ppsor\Desktop\studi\TCM\avalanches\avalanches_allbins_sameleng_dspl2_84_Maggiottt');
% QUESTO SOPRA RICOSTRUENDO DAI DATI BEAMFORMATI
% load('C:\Users\ppsor\Desktop\studi\TCM\avalanches\avalanches_allbins_Samelength_dwnsp4');
% load('C:\Users\ppsor\Desktop\TS_second_try_source_rec2_emah_dk_mni_coords\avalanches_secondtry_allbins_Samelength_NOdwnsp');
% load('C:\Users\ppsor\Desktop\studi\TCM\avalanches\avalanches_allbins_Samelength_NOdwnsp');
% load('C:\Users\ppsor\Desktop\studi\TCM\avalanches\avalanches_allbins_Samelength_dwnspl2');
% load('C:\Users\ppsor\Desktop\TS_third_try\third_try_source_rec2_pp_dk_mni_coords\avalanches_thirdTry_allbins_Samelength_NOdwnsp');
% load('C:\Users\ppsor\Desktop\studi\TCM\avalanches\avalanches_allbins_Samelength_dwnsp4');
% load('C:\Users\ppsor\Desktop\ms_maggio_2020_dk\avalanches_allbins_Maggio2020_NOdwnsp');
avalanches=avalanches{3}{6};
addpath('C:\Users\ppsor\Desktop\studi\TCM\codici_temp')

%% other zscore

% load('C:\Users\Pierpaolo\Desktop\Viktor\extra analyses\z-scores\aval_z25');
% avalanches=aval{1};
% clearvars aval

%%
min_size_aval=5;
cellsz = min(cell2mat(cellfun(@length,avalanches,'uni',false)));
clearvars temp
    for kk1=1:size(avalanches,2)
        idx=1;
        for kk2=1:size(avalanches{kk1},2)
%           if sum(sum(avalanches{kk1}{kk2},2))>min_size_aval
            if size(avalanches{kk1}{kk2},2)>min_size_aval
                curr_aval=avalanches{kk1}{kk2};
    %           ttt=delay_finder(curr_aval,1024/2,3);
    %           [~,ttt]=delay_from_starter(curr_aval);
                ttt=delay_finder_starter_2(curr_aval,512,3);
            if sum(ttt(:))~=0
                [temp{kk1}(:,:,idx)]=ttt;
                idx=idx+1;
            end
            clearvars ttt
            end    
        end
    end
     
%%
% structural data
% load('C:\Users\ppsor\Desktop\studi\TCM\data_MRI\tracts_DKT_25HC.mat');
% load('C:\Users\ppsor\Desktop\Code\DKT1reg_str')
% DKT1mm([1 2 9 10 11 14 15 18 19 20 29 32 33 34 37 38 39])=[];

% load('C:\Users\ppsor\Desktop\studi\TCM\data_MRI\workspace_tractographies_ms_spase.mat');
% sc=cat(3,mri_pzms_l,mri_ctrlms_l);
% W=cat(3,mri_pzms_tract,mri_ctrlms_tract);

load('C:\Users\ppsor\Desktop\studi\TCM\data_MRI\tracts_DKT_25HC.mat');
sc=cat(3,tract_DKT_L_MS, tract_DKT_L_CTRL);

%%

simmetric_cortex=[7 8 23:53 17 18 54:84];
sc=sc./1000; % to express it in meters

% sc=sc([7 8 17 18 23:53 54:84],[7 8 17 18 23:53 54:84],:);
% W=W(simmetric_cortex,simmetric_cortex,:);
sc=sc(simmetric_cortex,simmetric_cortex,:);
sc_all_ms=sc(:,:,1:18);
sc_all_ctrl=sc(:,:,19:38);
sc_ms=sum(sc(:,:,1:18),3)./sum(sc(:,:,1:18)~=0,3);
sc_ctrl=sum(sc(:,:,19:end),3)./sum(sc(:,:,19:end)~=0,3);
sc_ms(isnan(sc_ms))=0;
sc_ctrl(isnan(sc_ctrl))=0;

% delays per patient
for kk1=1:size(temp,2)
    aver_del(:,:,kk1)=sum(temp{kk1},3)./sum(temp{kk1}~=0,3);
end
aver_del(isnan(aver_del))=0;
aver_del=aver_del(simmetric_cortex,simmetric_cortex,:);

% figure
% subplot(1,2,1)
% imagesc(mean(aver_del,3));
% title('normal scale')
% subplot(1,2,2)
% imagesc(log10(mean(aver_del,3)));
% title('log')
%% run for weighted dels
% weight_del=zeros(size(W));
% weight_del(W~=0 & aver_del~=0)=W(W~=0 & aver_del~=0) .* aver_del(W~=0 & aver_del~=0);

% delays per group
dels_ms=sum(aver_del(:,:,1:18),3)./sum(aver_del(:,:,1:18)~=0,3);
dels_ctrl=sum(aver_del(:,:,19:end),3)./sum(aver_del(:,:,19:end)~=0,3);
dels_ms(isnan(dels_ms))=0;
dels_ctrl(isnan(dels_ctrl))=0;
dels_each_pz=aver_del(:,:,1:18);
dels_each_ctrl=aver_del(:,:,19:end);

% velocities
vels=zeros(size(sc));
mask= aver_del~=0 & sc~=0;
vels(mask)= sc(mask)./aver_del(mask);
vels_ms=sum(vels(:,:,1:18),3)./sum(vels(:,:,1:18)~=0,3);
vels_ctrl=sum(vels(:,:,19:end),3)./sum(vels(:,:,19:end)~=0,3);
vels_ms(isnan(vels_ms))=0;
vels_ctrl(isnan(vels_ctrl))=0;
vels_each_pz=vels(:,:,1:18);
vels_each_ctrl=vels(:,:,19:end);



% % % vels_ctrl=mean(vels(:,:,19:end));
%%
% vels_tot=sum(vels,3)./sum(vels~=0,3);
% vels_tot(isnan(vels_tot))=0;
% Cortical_Regions=[7 8 17 18 23:84];
% vels_tot=vels_tot(Cortical_Regions,:,:);
% vels_tot=vels_tot(:,Cortical_Regions,:);
% vels_tot=(vels_tot+vels_tot')/2;
% 
% vels_ctrl=vels_ctrl(Cortical_Regions,:,:);
% vels_ctrl=vels_ctrl(:,Cortical_Regions,:);
% vels_ctrl=(vels_ctrl+vels_ctrl')/2;
% 
% writematrix(round(vels_ctrl,3),'Avels_DKT_66_3.txt','Delimiter','tab');
%%
% correlazione percentili lunghezza con delay
for kk1=1:size(sc,3)
        mask_temp=mask(:,:,kk1);
        sc_temp=sc(:,:,kk1);
        dels_temp=aver_del(:,:,kk1); 
        perc_pos=zeros(size(sc,1),size(sc,2));
        discr=discretize(sc_temp(mask_temp),prctile(sc_temp(mask_temp),0:1:100));
        perc_pos(mask_temp)=discr;
        for kk2=1:max(discr)
            perc_len(kk1,kk2)=sum(sc_temp(perc_pos==kk2))/sum(sc_temp(perc_pos==kk2)~=0);  

             perc_del(kk1,kk2)=sum(dels_temp(perc_pos==kk2))/sum(dels_temp(perc_pos==kk2)~=0);  
%             perc_del(kk1,kk2)=mean(dels_temp(perc_pos==kk2));      
        end
end    

% e con le velocita'
for kk1=1:size(sc,3)
        mask_temp=mask(:,:,kk1);
        sc_temp=sc(:,:,kk1);
        vels_temp=vels(:,:,kk1);
        perc_pos=zeros(size(sc,1),size(sc,2));
        discr=discretize(sc_temp(mask_temp),prctile(sc_temp(mask_temp),0:1:100));
        perc_pos(mask_temp)=discr;
        for kk2=1:max(discr)
            perc_vel(kk1,kk2)=sum(vels_temp(perc_pos==kk2))/sum(vels_temp(perc_pos==kk2)~=0);      
        end
end   


%% velocita' omogenee
Omo_vels=[2 5 7 10];
for kk_vels=1:size(Omo_vels,2)   
    for kk1=1:size(sc,3)
            mask_temp=mask(:,:,kk1);
            sc_temp=sc(:,:,kk1);
            vels_temp=sc_temp./Omo_vels(kk_vels); % normalmente qua lo facevo con le velocita'
            perc_pos=zeros(size(sc,1),size(sc,2));
            discr=discretize(sc_temp(mask_temp),prctile(sc_temp(mask_temp),0:1:100));
            perc_pos(mask_temp)=discr;
            for kk2=1:max(discr)
                perc_delOmog(kk1,kk2,kk_vels)=sum(vels_temp(perc_pos==kk2),'all')/sum(vels_temp(perc_pos==kk2)~=0,'all');      
            end
    end    
end

%% surrogati valanghe

% computation of randomized matrices
group_generator = @(x)x(randperm(numel(x)));

num_perm=100;
sum_per_goup=zeros(84,84);
for kk1=1:size(avalanches,2)
    sum_per_pat=zeros(84,84);
    c=cellfun('size',avalanches{kk1},2);
    long_aval=sum(c>min_size_aval);
    if long_aval>0
        big_avals=1;
        for kk2=1:size(avalanches{kk1},2)
            if size(avalanches{kk1}{kk2},2)>min_size_aval
                curr_aval=avalanches{kk1}{kk2};
                temp=zeros(size(curr_aval,1),size(curr_aval,1),num_perm);
                for kk3=1:num_perm
                    rand_temp=group_generator(1:size(curr_aval,2));
                    curr_aval_rand=curr_aval(:,rand_temp);
                    temp(:,:,kk3)=delay_finder_starter_2(curr_aval_rand,512,3);
                    clearvars curr_diff_rand
                end
                tm_per_aval_t(:,:,big_avals)=reshape(permute(temp,[1 3 2]),[],size(temp,2),1);
                big_avals=big_avals+1;
            end
        end
        tm_per_pat_t(:,:,kk1)=mean(tm_per_aval_t,3);
        clearvars  hh_struct tm_per_aval_t   
    end
end

% surrogate velocities
idx=1;
for kk_perm=1:84:size(tm_per_pat_t,1)
    probe=tm_per_pat_t(kk_perm:kk_perm+size(tm_per_pat_t,2)-1,:,19:end);
    for kk1=1:size(probe,3)
        sc_temp=sc(:,:,kk1);
        probe_temp=probe(:,:,kk1);
        dels_temp=zeros(size(sc_temp));
        vels_temp=zeros(size(sc_temp));
        mask_temp=sc_temp~=0 & probe_temp~=0;
        dels_temp(mask_temp)=probe_temp(mask_temp); % delays
        vels_temp(mask_temp)=sc_temp(mask_temp)./probe_temp(mask_temp); % velocities
        perc_pos=zeros(size(sc,1),size(sc,2));
        discr=discretize(sc_temp(mask_temp),prctile(sc_temp(mask_temp),0:1:100));
        perc_pos(mask_temp)=discr;
        for kk2=1:max(discr)
            perc_del_perm(kk1,kk2,idx)=sum(dels_temp(perc_pos==kk2),'all')/sum(dels_temp(perc_pos==kk2)~=0,'all');      
            perc_vel_perm(kk1,kk2,idx)=sum(vels_temp(perc_pos==kk2),'all')/sum(vels_temp(perc_pos==kk2)~=0,'all');      
        end
    end   
    idx=idx+1;
end

perc_vel_perm_Mean=squeeze(mean(perc_vel_perm,1));
perc_vel_perm_Std=squeeze(std(perc_vel_perm,1));

% figure
% plot(perc_vel_perm_Mean,'b*')
perc_del_perm_Mean=squeeze(mean(perc_del_perm,1));
perc_del_perm_Std=squeeze(std(perc_del_perm,1));

% figure
% plot(perc_del_perm_Mean,'b*')
% velocities calculated with mean delays

mean_del=zeros(size(aver_del));
    for kk1=1:size(aver_del,3)
        mean_del(:,:,kk1)=sum(aver_del(:,:,kk1),'all')/sum(aver_del(:,:,kk1)~=0,'all');
    end
vels_cost_del=zeros(size(aver_del));

   for kk1=1:size(aver_del,3)
       del_pz=unique(mean_del(:,:,kk1));
       t=sc(:,:,kk1);
       mask_t=t~=0;
       yy=zeros(84,84);
       yy(mask_t)=t(mask_t)./mean_del(mask_t);
       vels_cost_del(:,:,kk1)=yy;
    end

for kk1=1:size(vels_cost_del,3)
    sc_temp=sc(:,:,kk1);
    probe_temp=vels_cost_del(:,:,kk1);
    perc_pos=zeros(size(sc,1),size(sc,2));
    discr=discretize(sc_temp(sc_temp~=0 & probe_temp~=0),prctile(sc_temp(sc_temp~=0 & probe_temp~=0),0:1:100));
    perc_pos(sc_temp~=0 & probe_temp~=0)=discr;
    for kk2=1:max(discr)
        perc_vels_cost_del(kk1,kk2)=sum(probe_temp(perc_pos==kk2),'all')/sum(probe_temp(perc_pos==kk2)~=0,'all');      
    end
end   

perc_vels_cost_del_Mean=mean(perc_vels_cost_del,1);

% figure
% plot(perc_vels_cost_del_Mean)
%% FIGURE 2
% vels_ctrl=sc_ctrl./dels_ctrl;
simmetric_cortex=[7 8 23:53 17 18 54:84];
sc_ctrl=sc_ctrl(simmetric_cortex,simmetric_cortex);
dels_ctrl=dels_ctrl(simmetric_cortex,simmetric_cortex);
vels_ctrl=vels_ctrl(simmetric_cortex,simmetric_cortex);
%
[r,p]=corr(sc_ctrl(sc_ctrl~=0 & dels_ctrl~=0),dels_ctrl(sc_ctrl~=0 & dels_ctrl~=0))
%%
figure
set(gcf, 'InvertHardcopy', 'on')

% sgtitle('THE HUMAN DELAYOME')
subplot(4,2,1)
sc_ctrl(sc_ctrl==0)=NaN;
him1=imagesc(sc_ctrl);
set(him1,'alphadata',~isnan(sc_ctrl))
sc_ctr(isnan(sc_ctrl))=0;
colorbar
title('TRACK LENGTHS (meters)','FontWeight','bold');
xlabel('BRAIN REGIONS','FontWeight','bold');
ylabel('BRAIN REGIONS','FontWeight','bold');
xticks([])
yticks([])
axis square

subplot(4,2,2)
hist(sc_ctrl(sc_ctrl~=0),80)
title('DISTRIBUTION OF LENGTHS','FontWeight','bold');
xlabel('TRACK LENGTHS (meters)','FontWeight','bold');
yticklabels([])
xlim([0.03 0.3])
ax=gca;
ax.FontWeight='bold';
axis square

subplot(4,2,3)
imagesc(log10(dels_ctrl))
colorbar
title('TIME DELAYS (sec - log scale)')
xlabel('BRAIN REGIONS','FontWeight','bold');
ylabel('BRAIN REGIONS','FontWeight','bold');
xticks([])
yticks([])
axis square

subplot(4,2,4)
hist(log10(dels_ctrl(dels_ctrl~=0)),40)
title('TIME DELAYS (sec - log scale)','FontWeight','bold');
xlabel('DELAYS (sec)','FontWeight','bold');
yticklabels([])
ax=gca;
ax.FontWeight='bold';
axis square

subplot(4,2,5)
vels_ctrl(vels_ctrl==0)=NaN;
him2=imagesc(vels_ctrl);
set(him2,'alphadata',~isnan(vels_ctrl))
vels_ctrl(isnan(vels_ctrl))=0;
colorbar
title('VELOCITIES  (m/sec)')
xlabel('BRAIN REGIONS','FontWeight','bold');
ylabel('BRAIN REGIONS','FontWeight','bold');
xticks([])
yticks([])
axis square

subplot(4,2,6)
hist(vels_ctrl(vels_ctrl~=0),80)
title('VELOCITIES','FontWeight','bold');
xlabel('VELOCITIES  (m/sec)','FontWeight','bold');
yticklabels([])
ax=gca;
ax.FontWeight='bold';
axis square


%% FIGURE 3

num_perm=100;
sym_whole_pos=[1:10 21:53 11:20 54:84];
% sym_whole_pos=[1:66];
sc_ctrl_LR=sc_ctrl(sym_whole_pos,sym_whole_pos);
dels_ctrl_LR=dels_ctrl(sym_whole_pos,sym_whole_pos);
vels_ctrl_LR=vels_ctrl(sym_whole_pos,sym_whole_pos);
maskISO=zeros(size(sc_ctrl));
maskISO(1:42,1:42)=1;
maskISO(43:84,43:84)=1;
% maskISO(1:33,1:33)=1;
% maskISO(34:66,34:66)=1;
maskCROSS=zeros(size(sc_ctrl));
maskCROSS(1:42,43:84)=1;
maskCROSS(43:84,1:42)=1;
A=mean(sc_ctrl_LR,3);
B=mean(dels_ctrl_LR,3);
C=sc_ctrl_LR./B;
D=mean(C,3);
D(D==Inf)=0;
D(isnan(D))=0;
g=mean(perc_len(19:end,:),1);
gg=mean(perc_del(19:end,:),1);
ggg=mean(perc_vel(19:end,:),1);

% maskCROSS(1:33,34:66)=1;
% maskCROSS(34:66,1:33)=1;

len_del_cors=zeros(size(dels_each_ctrl,3),1);
len_del_cors_perm=zeros(size(dels_each_ctrl,3),num_perm);
% mask = tril(true(size(sc_ctrl)),-1);
group_generator = @(x)x(randperm(numel(x)));
for kk1=1:size(dels_each_ctrl,3)
    sc_temp=sc_all_ctrl(:,:,kk1);
    del_temp=dels_each_ctrl(:,:,kk1);
    len_del_cors(kk1)=corr(sc_temp(sc_temp~=0 & del_temp~=0),del_temp(sc_temp~=0 & del_temp~=0),'type','S');
    for kk_perm=1:num_perm
        temp_pos=group_generator(1:size(sc_temp(sc_temp~=0 & del_temp~=0),1)); 
        len_del_cors_perm(kk1,kk_perm)=corr(sc_temp(sc_temp~=0 & del_temp~=0),del_temp(temp_pos)','type','S');
    end
end
%%
% only do this if you want renormalized lines
sc_inorm=zeros(84);
sc_inorm(sc_ctrl_LR~=0)=inormal(sc_ctrl_LR(sc_ctrl_LR~=0));
dels_inorm=zeros(84);
dels_inorm(dels_ctrl_LR~=0)=inormal(dels_ctrl(dels_ctrl_LR~=0));

%%

figure
subplot(2,3,1)
scatter(sc_ctrl_LR(sc_ctrl_LR~=0 & dels_ctrl_LR~=0 & maskISO),dels_ctrl_LR(sc_ctrl_LR~=0 & dels_ctrl_LR~=0 & maskISO),5,'o','b')
% scatter(sc_inorm(sc_inorm~=0 & dels_inorm~=0 & maskISO),dels_inorm(sc_inorm~=0 & dels_inorm~=0 & maskISO),5,'o','b')
hold on
scatter(sc_ctrl_LR(sc_ctrl_LR~=0 & dels_ctrl_LR~=0 & maskCROSS),dels_ctrl_LR(sc_ctrl_LR~=0 & dels_ctrl_LR~=0 & maskCROSS),5,'o','r')
% scatter(sc_inorm(sc_inorm~=0 & dels_inorm~=0 & maskCROSS),dels_inorm(sc_inorm~=0 & dels_inorm~=0 &maskCROSS),5,'o','r')
title({'CORRELATION LENGTHS/DELAYS' , 'GROUP LEVEL'},'FontSize',18,'FontWeight','bold')
xlabel('LENGTHS (m)','FontSize',14,'FontWeight','bold')
ylabel('DELAYS (sec)','FontSize',14,'FontWeight','bold')
xlim([0.03,0.3]);

ylim([0.006,0.04]);
% legend('AutoUpdate','off')
% legend('ISO','CROSS','WHOLE')
ax=gca;
ax.FontSize=16;
ax.FontWeight='bold';
% coeffs = polyfit(sc_inorm(sc_inorm~=0 & dels_inorm~=0 & maskISO),dels_inorm(sc_inorm~=0 & dels_inorm~=0 & maskISO),1);
% coeffs = polyfit(sc_ctrl_LR(sc_ctrl_LR~=0 & dels_ctrl_LR~=0),log10(dels_ctrl(sc_ctrl_LR~=0 & dels_ctrl_LR~=0)),1);
coeffs = polyfit(sc_ctrl_LR(sc_ctrl_LR~=0 & dels_ctrl_LR~=0),dels_ctrl(sc_ctrl_LR~=0 & dels_ctrl_LR~=0),1);
hline = refline(coeffs);
hline.Color = 'k';
hline.LineWidth=2;
legend('ISO','CROSS','WHOLE')
hold off

subplot(2,3,[2 3])
hold on
h1=plot(len_del_cors_perm,'b.','MarkerSize',10);
h2=plot(len_del_cors,'r.','MarkerSize',10);
title({'CORRELATION LENGTH/DELAYS', 'SUBJECT-SPECIFIC'},'FontSize',18,'FontWeight','bold');
xlabel('SUBJECT','FontSize',14,'FontWeight','bold');
ylabel('CORRELATION COEFFICIENT','FontSize',14,'FontWeight','bold');
xticks([]);
xlim([0 21]);
ax=gca;
ax.FontSize=16;
ax.FontWeight='bold';
legend([h1(1), h2(1)],'RANDOMIZED','OBSERVED')
hold off

subplot(2,3,4)
hold on
plot(mean(perc_del(19:end,:),1),'bo','MarkerSize',8,'MarkerFaceColor','b')
plot(squeeze(mean(perc_delOmog(19:end,:,:),1)),'LineWidth',2)
title({'EMPIRICAL DELAYS','vs BY CONSTANT VELOCITIES'}','FontSize',18,'FontWeight','bold');
xlabel('%ILE OF TRAIT LENGTH','FontSize',14,'FontWeight','bold');
ylabel('MEAN DELAY (sec)','FontSize',14,'FontWeight','bold');
legend({'OBSERVED DELAYS','10 m/sec','7 m/sec','5 m/sec','2 m/sec'},'location','northwest')
ax=gca;
ax.FontSize=16;
ax.FontWeight='bold';
hold off


subplot(2,3,5)
hold on
plot(mean(perc_vel(19:end,:),1),'LineWidth',4,'Color','r','LineStyle','-')
title({'CORRELATION LENGHTS/VELOCITIES', 'NULL-VALIDATION'},'FontSize',18,'FontWeight','bold');
xlabel('%ile of trait length','FontSize',14,'FontWeight','bold');
ylabel('mean velocity (m/sec)','FontSize',14,'FontWeight','bold');
yticklabels([]);
xx = 1:numel(mean(perc_vel(19:end,:),1));
plot(mean(perc_vel_perm_Mean,2))
curve1 = mean(perc_vel_perm_Mean,2)' + mean(perc_vel_perm_Std,2)';
curve2 = mean(perc_vel_perm_Mean,2)' - mean(perc_vel_perm_Std,2)';
x2 = [xx, fliplr(xx)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'b','facealpha',.1);
axis square
legend({'OBSERVED VELOCITIES','VELOCITIES - NULL DELAYS','INTERVALS OF CONFIDENCE'},'location','northwest','FontSize',12,'FontWeight','bold')
ax=gca;
ax.FontSize=16;
ax.FontWeight='bold';
hold off


subplot(2,3,6)
hold on
scatter(sc_ctrl_LR(maskISO==1 & D~=0),D(maskISO==1 & D~=0),5,'o','b')
scatter(sc_ctrl_LR(maskCROSS==1 & D~=0),D(maskCROSS==1 & D~=0),5,'o','r')
title({'CORRELATION LENGTHS/VELOCITIES', 'GROUP LEVEL'},'FontSize',18,'FontWeight','bold')
xlabel('LENGTHS (m)','FontSize',14,'FontWeight','bold')
ylabel('VELOCITY (m/sec)','FontSize',14,'FontWeight','bold')
xlim([0.03,0.3]);
legend('ISO','CROSS','location','northwest')
ax=gca;
ax.FontSize=16;
ax.FontWeight='bold';
coeffs = polyfit(sc_ctrl_LR(sc_ctrl_LR~=0 & D~=0),D(sc_ctrl_LR~=0 & D~=0),1);
hline = refline(coeffs);
hline.Color = 'k';
hline.LineWidth=2;
legend('ISO','CROSS','WHOLE')
hold off
hold off

%% confronti con MS
addpath('C:\Users\ppsor\Desktop\Code\violinplotfigure')
load('C:\Users\ppsor\Desktop\studi\TCM\data_MRI\lesions2');
MEDIA= sum(dels_each_ctrl,3)./sum(dels_each_ctrl~=0,3);
dels_each_ctrl(dels_each_ctrl==0)=NaN;
STD= nanstd(dels_each_ctrl,[],3);
dels_each_ctrl(isnan(dels_each_ctrl))=0;
reliab_dels=zeros(size(dels_each_ctrl));
% Get a logical vector of the locations where the value is more than 3 sd away from the mean.
withinMean = abs(dels_each_ctrl - MEDIA) < repmat(3.*STD,1,1,20);
% Extract only those elements
reliab_dels(withinMean) = dels_each_ctrl(withinMean);

% mask_reliab=repmat(sum(dels_each_ctrl~=0,3)>10,1,1,size(dels_each_pz,3));
dels_ctrlRel=sum(reliab_dels,3)./sum(reliab_dels~=0,3);
dels_ctrlRel(isnan(dels_ctrlRel))=0;

% lesions=lesions(simmetric_cortex,simmetric_cortex,:);
mask_les=lesions~=0 & dels_each_pz~=0;
diff_dels=zeros(size(dels_each_pz,3));
temp1=repmat(dels_ctrlRel,1,1,size(dels_each_pz,3));

%% randomizzazione per pannello 4
diff_dels=dels_each_pz - repmat(dels_ctrlRel,1,1,size(dels_each_pz,3));
% health_diffs=diff_dels(~mask_les & dels_each_pz~=0);
% lesion_diffs=diff_dels(mask_les);
% num_edge_les=numel(lesion_diffs);
% obs_diff=mean(diff_dels(mask_les));
% for kk1=1:1000
%     rand_mean_diffs(kk1)=mean(health_diffs(randperm(numel(health_diffs),num_edge_les)));
% end

%%

load('C:\Users\ppsor\Desktop\studi\TCM\Figure\Figure 4\Figure4_data.mat')
statistics=length(find(rand_mean_diffs>obs_diff))/1000;
addpath(genpath('C:\Users\ppsor\Desktop\Code\raincloud_plots\'));
[cb] = cbrewer('qual', 'Set3', 12, 'pchip');
cl(1, :) = cb(4, :);
cl(2, :) = cb(1, :);
fig_position = [200 200 600 400]; % coordinates for figures    

figure
set(gcf,'color','w');
subplot(1,4,1)
hold on
plot(mean(perc_del(1:18,:),1),'r.','MarkerSize',12)
plot(mean(perc_del(19:end,:),1),'b.','MarkerSize',12)
title('DELAYS','FontSize',20,'FontWeight','bold');
xlabel('%ile of trait length','FontSize',20,'FontWeight','bold');
ylabel('MEAN DELAY (sec)');
legend('MS','CTRL','location','Northwest','FontSize',14,'FontWeight','bold');
axis square
ax=gca;
ax.FontSize=16;
ax.FontWeight='bold';
hold off

subplot(1,4,2)
hold on 
[f1,x1]=ecdf(mean(perc_del(1:18,:),1));
plot(x1,f1,'r','LineWidth',2)
[f2,x2]=ecdf(mean(perc_del(19:end,:),1));
plot(x2,f2,'b','LineWidth',2)
title('DELAYS - eCDF','FontSize',20,'FontWeight','bold');
xlabel('DELAYS (sec)','FontSize',20,'FontWeight','bold');
% legend([ 'KStest=' num2str(kstest2(mean(perc_del(1:18,:),1),mean(perc_del(19:end,:),1)))],'location','Northwest');
axis square
legend('MS','CTRLS','location','Northwest','FontSize',16,'FontWeight','bold');
axis square
ax=gca;
ax.FontSize=16;
ax.FontWeight='bold';
hold off

subplot(1,4,3)
h_1=violinplot(diff_dels(mask_les));
h_1(1,1).ViolinPlot.FaceColor=[1 0 0];
xticklabels([]);
axis square
title({'EDGE-WISE DELAY DIFFERENCES','patients - controls'},'FontSize',24,'FontWeight','bold')
ylabel('DELAY DIFFERENCE','FontSize',20,'FontWeight','bold');
fig = gcf;
fig.InvertHardcopy = 'off';
ax=gca;
ax.FontSize=16;
ax.FontWeight='bold';
axis square

subplot(1,4,4)
h2 = raincloud_plot(rand_mean_diffs, 'box_on', 0, 'box_dodge', 0, 'box_dodge_amount',...
0, 'dot_dodge_amount', .3, 'color', cb(1,:), 'cloud_edge_col', cb(1,:));
title('MEAN DELAY DIFFERENCE','FontSize',20);
xline(obs_diff,'r',{'OBSERVED DIFFERENCE'},'FontSize',16)
yticks([]);
ax=gca;
ax.FontSize=16;
ax.FontWeight='bold';
axis square
box off

%%

g=mean(perc_len(19:end,:),1);
gg=mean(perc_del(19:end,:),1);
[r,p]=corr(g',gg');

figure
scatter(g,gg,'.')
lsline
xlim([0.03,0.25])


[r,p]=corr(sc_ctrl_LR(sc_ctrl_LR~=0 & dels_ctrl_LR~=0),dels_ctrl_LR(sc_ctrl_LR~=0 & dels_ctrl_LR~=0));



