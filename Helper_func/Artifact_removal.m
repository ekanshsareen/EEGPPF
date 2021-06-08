% Usage:
%   >> [art, horiz, vert, blink, disc,...
%         soglia_DV, diff_var, soglia_K, med2_K, meanK, soglia_SED, med2_SED, SED, soglia_SAD, med2_SAD, SAD, ...
%         soglia_GDSF, med2_GDSF, GDSF, soglia_V, med2_V, nuovaV, soglia_D, maxdin]=ADJUST (EEG,out)
% Inputs:
%   EEG        - current dataset structure or structure array (has to be epoched)
%   out        - (string) report file name 
% Outputs:
%   art        - List of artifacted ICs
%   horiz      - List of HEM ICs 
%   vert       - List of VEM ICs   
%   blink      - List of EB ICs     
%   disc       - List of GD ICs     
%   soglia_DV  - SVD threshold      
%   diff_var   - SVD feature values
%   soglia_K   - TK threshold      
%   meanK      - TK feature values
%   soglia_SED - SED threshold      
%   SED        - SED feature values
%   soglia_SAD - SAD threshold      
%   SAD        - SAD feature values
%   soglia_GDSF- GDSF threshold      
%   GDSF       - GDSF feature values
%   soglia_V   - MEV threshold      
%   nuovaV     - MEV feature values
%
% function [art, horiz, vert, blink, disc,...
%         soglia_DV, diff_var, soglia_K, meanK, soglia_SED, SED, soglia_SAD, SAD, ...
%         soglia_GDSF, GDSF, soglia_V, nuovaV, soglia_D, maxdin]=ADJUST (EEG,out)
function [art, horiz, vert, blink, disc,...
        soglia_DV, diff_var, soglia_K, med2_K, meanK, soglia_SED, med2_SED, SED, soglia_SAD, med2_SAD, SAD, ...
        soglia_GDSF, med2_GDSF, GDSF, soglia_V, med2_V, nuovaV, soglia_D, maxdin]=Artifact_removal (EEG,out)  
%% Settings
if length(size(EEG.data))==3   
    num_epoch=size(EEG.data,3);      
else   
    num_epoch=1;
end
topografie=EEG.icawinv'; %computes IC topographies
for i=1:size(EEG.icawinv,2) % number of ICs  
    ScalingFactor=norm(topografie(i,:));
    topografie(i,:)=topografie(i,:)/ScalingFactor;
    if length(size(EEG.data))==3
        EEG.icaact(i,:,:)=ScalingFactor*EEG.icaact(i,:,:);
    else
        EEG.icaact(i,:)=ScalingFactor*EEG.icaact(i,:);
    end  
end
blink=[]; horiz=[]; vert=[]; disc=[];
%% Check EEG channel position information
nopos_channels=[];
pos_channels=setdiff(1:length(EEG.chanlocs),nopos_channels);

%% Feature extraction
disp(' ')
disp('Features Extraction:')
% GDSF - General Discontinuity Spatial Feature
disp('GDSF - General Discontinuity Spatial Feature...')
GDSF = compute_GD_feat(topografie,EEG.chanlocs(1,pos_channels),size(EEG.icawinv,2));

% SED - Spatial Eye Difference
disp('SED - Spatial Eye Difference...')
[SED,medie_left,medie_right]=computeSED_NOnorm(topografie,EEG.chanlocs(1,pos_channels),size(EEG.icawinv,2)); 

% SAD - Spatial Average Difference
disp('SAD - Spatial Average Difference...')
[SAD,var_front,var_back,mean_front,mean_back]=computeSAD(topografie,EEG.chanlocs(1,pos_channels),size(EEG.icawinv,2));

% SVD - Spatial Variance Difference between front zone and back zone
diff_var=var_front-var_back;

% epoch dynamic range, variance and kurtosis
K=zeros(num_epoch,size(EEG.icawinv,2)); %kurtosis
Kloc=K;
Vmax=zeros(num_epoch,size(EEG.icawinv,2)); %variance
for i=1:size(EEG.icawinv,2) % number of ICs
    
    for j=1:num_epoch              
        Vmax(j,i)=var(EEG.icaact(i,:,j));        
        K(j,i)=kurt(EEG.icaact(i,:,j));
    end  
end

disp('Temporal Kurtosis...')
meanK=zeros(1,size(EEG.icawinv,2));
for i=1:size(EEG.icawinv,2)
    if num_epoch>100
    meanK(1,i)=trim_and_mean(K(:,i)); 
    else meanK(1,i)=mean(K(:,i));
    end
end

% MEV - Maximum Epoch Variance
disp('Maximum epoch variance...')
maxvar=zeros(1,size(EEG.icawinv,2));
meanvar=zeros(1,size(EEG.icawinv,2));
for i=1:size(EEG.icawinv,2)
    if num_epoch>100
     maxvar(1,i)=trim_and_max(Vmax(:,i)');
     meanvar(1,i)=trim_and_mean(Vmax(:,i)');
    else 
     maxvar(1,i)=max(Vmax(:,i));
     meanvar(1,i)=mean(Vmax(:,i));
    end
end
% MEV in reviewed formulation:
nuovaV=maxvar./meanvar;

%% Thresholds computation
disp('Computing thresholds...')
[soglia_K,med1_K,med2_K]=EM(meanK);
[soglia_SED,med1_SED,med2_SED]=EM(SED);
[soglia_SAD,med1_SAD,med2_SAD]=EM(SAD);
[soglia_GDSF,med1_GDSF,med2_GDSF]=EM(GDSF);
[soglia_V,med1_V,med2_V]=EM(nuovaV); 
%% Horizontal eye movements (HEM)
horiz=intersect(intersect(find(SED>=soglia_SED),find(medie_left.*medie_right<0)),...
    (find(nuovaV>=soglia_V)));

%% Vertical eye movements (VEM)
vert=intersect(intersect(find(SAD>=soglia_SAD),find(medie_left.*medie_right>0)),...
    intersect(find(diff_var>0),find(nuovaV>=soglia_V)));

%% Eye Blink (EB)
blink=intersect ( intersect( find(SAD>=soglia_SAD),find(medie_left.*medie_right>0) ) ,...
    intersect ( find(meanK>=soglia_K),find(diff_var>0) ));
%% Generic Discontinuities (GD)
disc=intersect(find(GDSF>=soglia_GDSF),find(nuovaV>=soglia_V));
aic=unique([blink disc horiz vert]);
art = nonzeros( union (union(blink,horiz) , union(vert,disc)) )'; %artifact ICs
soglia_D=0;
soglia_DV=0;
maxdin=zeros(1,size(EEG.icawinv,2));
return