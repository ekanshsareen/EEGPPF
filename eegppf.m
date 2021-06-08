% eegppf() - Matlab code for filtering and pre-processing of EEG data built
%   by incorporating ADJUST toolbox (Marco Buiatti et al.) developed at
%   Cognitive Neuroimaging Unit, INSERM - as eegppf() at Signal processing
%   Bio-medical Imaging Lab (SBILab), IIIT-Delhi.
%   User feedback welcome: email ekansh15139@iiitd.ac.in
%
% Note: The pipleine is linked to Data in Brief article titled "EEG dataset of individuals with 
%       intellectual and developmental disorder and healthy controls under rest and music stimuli"
%
% Authors:
%   Ekansh Sareen and Lakshya Singh
%
% Requirements:
%   MATLAB, EEGLAB toolbox, and CleanLine plugin.
%
% Description:
%   eegppf() can be directly used to filter and pre-process raw EEG data.
%   It follows five-step procedure:
%   1) Load data
%   2) Filter data
%   3) Channel location file
%   4) Independent Component Analysis
%   5) Removal of components
%   The code automatically detects the spatial and temporal features of
%   each of the ICs. Amongst the determined ICs, omly the worst two are
%   removed which is also automatically detected in this code. This is kept
%   to control the loss of data while working with lesser channel data.
%   The additional codes required for eegppf() to run are provided under
%   'Helper_func' folder.
%   Usage:
%     (1) Add the pipeline and main EEGLAB folder in the working directory
%     path. 
%     (2) Required arguments
%     File_name='sample_dataset.set'; (Sample file is provided)
%     File_path='../pipeline/Helper_func/sample_data/';
%     chan_loc='../pipeline/Helper_func/emotiv_electrode.ced';
%     Fs=128;
%     lpf=1;
%     hpf=60;
%     linenoise = 50;
%     electrode_select={'AF3','F7','F3','FC5','T7','P7','O1','O2','P8','T8','FC6','F4','F8','AF4'};
%     rem_num=2;
%     >>clean_data=eegppf(File_name,File_path,chan_loc,Fs,lpf, hpf,linenoise,electrode_select,rem_num)
%     
%     Note: If the function does not work try executing eeglab command to initialize all the required variables, 
%           close the GUI, and then run eegppf() with required arguments as mentioned in step-(2).
% 
%
% If you use this pipeline for your research work, kindly cite these related articles:
% (1) Sareen, E., Singh, L., Varkey, B., Achary, K., Gupta, A. (2020). 
%     EEG dataset of individuals with intellectual and developmental disorder
%     and healthy controls under rest and music stimuli. Data in Brief, 105488, ISSN 2352-3409,
%     DOI:https://doi.org/10.1016/j.dib.2020.105488.
% (2) Sareen, E., Gupta, A., Verma, R., Achary, G. K., Varkey, B (2019), 
%     Studying functional brain networks from dry electrode EEG set during music and resting states
%     in neurodevelopment disorder, bioRxiv 759738 [Preprint]. 
%     Available from: https://www.biorxiv.org/content/10.1101/759738v1
%
% Copyright (C) 2020, Ekansh Sareen (1) and Lakshya Singh (1), 
% (1)SBILab, Indraprastha Institite of Information and Technology - Delhi
% (IIIT-D). Email: ekansh15139@iiitd.ac.in, sbilab@iiitd.ac.in

function clean_data=eegppf(File_name,File_path,chan_loc,Fs,lpf,hpf,linenoise,electrode_select,rem_num)   
    if nargin <9
        rem_num=2;
    end   
    if nargin <7
        linenoise=-1;   % Default value
    end
    if nargin <6
        lpf=1;          % Default value
        hpf=0.4*Fs;
    end
    disp('The pipeline is now in action!')
    disp('Sit back and let it clean the data for you..')
% Initializing variables
    ALLEEG=[];
    EEG=[];
    CURRENTSET=0;
    ALLCOM={};
    
% Loading data
    File_details=strsplit(File_name,'.');
    File_extension=File_details{length(File_details)};
    % For EDF file
    if strcmp(File_extension,'edf')
        EEG = pop_biosig(full_filename);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off');
    % For Mat file
    elseif strcmp(File_extension,'mat')
        % nbchan: No. of channels in the EEG data (automatically selects from data, default is '0')
        % srate: Sampling rate of the EEG data 
        EEG = pop_importdata('dataformat','matlab','nbchan',0,'data',full_filename,'srate',Fs);
    % For set file
    elseif strcmp(File_extension,'set')
        EEG = pop_loadset('filename',File_name,'filepath',File_path);
    else
        disp('extension not supported');
    end
    
% Channel locations
    % Add the channel location file for the EEG device
    % Custom channel location file for EMOTIV EPOC has been added (sample)
    EEG=pop_chanedit(EEG,'load',{chan_loc , 'filetype' ,'autodetect'});
    
% Select channels     
    EEG = pop_select( EEG,'channel',electrode_select);

% Filtering
    % Specify the filter band for the bandpass filter usinhg lpf and hpf
    EEG = pop_eegfiltnew(EEG, 'locutoff',lpf,'hicutoff',hpf);
    % Notch filter
    % This part executes only if hpf>50Hz
    if all(hpf>linenoise) && all(linenoise~=-1)
        % Line noise of 50Hz
        if linenoise==50
            EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:14] ,'computepower',1,...
                'linefreqs',linenoise,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,...
                'scanforlines',1,'sigtype','Channels','tau',100,'verb',1,'winsize',4,'winstep',1);
        % Line noise of 60Hz    
        elseif linenoise==60
            EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:14] ,'computepower',1,...
                'linefreqs',linenoise,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,...
                'scanforlines',1,'sigtype','Channels','tau',100,'verb',1,'winsize',4,'winstep',1);
        end
    end 
    % If you do not have CleanLine plugin, you can use the command
    % mentioned below to implement a notch filter centered at the
    % linenoise. Make sure to comment the above section for linenoise
    % removal using CleanLine before using this alternate option. 
    % EEG = pop_eegfiltnew(EEG, 'locutoff',linenoise-2,'hicutoff',linenoise+2, 'revfilt', 1);
    
% Independent Component Analysis (ICA)
    EEG = pop_runica(EEG,'extended',1);
    report_name='report.txt';
    lag=(EEG.xmax-EEG.xmin);
    disp(['Continuous dataset epoched in ' num2str(lag) ' sec long epochs to compute feature statistics over time']);
    % If epoch is not present: ifepoc=0
    % Change ifepoch if epochs are present in the data
    ifepoc=0;
    index=1;
    if ifepoc==0
        nevents=length(EEG.event);
        EEG.event(index+nevents).type='Full';
        EEG.event(index+nevents).latency=1+(index-1)*lag*EEG.srate; %EEG.srate is the sampling frequency
        EEG=eeg_checkset(EEG,'eventconsistency');
    end
    EEGep = pop_epoch( EEG, {  'Full'  }, [0 lag], 'newname', [EEG.setname '_ep5'] , 'epochinfo', 'yes');
    EEGep = eeg_checkset(EEGep);
    if isempty(EEGep.icaact)
        disp('Warning: EEG.icaact missing! Recomputed from EEG.icaweights, EEG.icasphere and EEG.data');
        EEGep.icaact = reshape(EEGep.icaweights*EEGep.icasphere*reshape(EEGep.data(1:size(EEGep.icaweights,1),:,:),[size(EEGep.icaweights,1) size(EEGep.data,2)*size(EEGep.data,3)]),[size(EEGep.icaweights,1) size(EEGep.data,2) size(EEGep.data,3)]);
    end
    
 % Artifact removal  
    [~, ~, ~, ~, ~, ~, ~, soglia_K, ~, meanK, soglia_SED, ~, SED, soglia_SAD, ~, SAD, ...
        soglia_GDSF, ~, GDSF, soglia_V, ~, nuovaV, ~, ~]=Artifact_removal(EEGep,report_name);
  
 % Description of output arguments of TK, SED, SAD, GDSF and MEV
    %   soglia_K   - TK threshold, meanK - TK feature values
    %   soglia_SED - SED threshold, SED - SED feature values
    %   soglia_SAD - SAD threshold, SAD - SAD feature values
    %   soglia_GDSF - GDSF threshold, GDSF - GDSF feature values
    %   soglia_V - MEV threshold, nuovaV - MEV feature values
    
    %   Removing components by sorting
    disp('Computing the artifacts from spatial and temporal features...')
    num_chan=length(electrode_select);
    Sorting_mat_ICA=zeros(num_chan,1);
    comp_mat=zeros(num_chan,5);
    comp_mat(:,1)=(-soglia_K+meanK)./soglia_K;
    comp_mat(:,2)=(-soglia_SED+SED)./soglia_SED;
    comp_mat(:,3)=(-soglia_SAD+SAD)./soglia_SAD;
    comp_mat(:,4)=(-soglia_GDSF+GDSF)./soglia_GDSF;
    comp_mat(:,5)=(-soglia_V+nuovaV)./soglia_V;
    for i =1:num_chan
        for k=1:5
            if comp_mat(i,k)>0
                Sorting_mat_ICA(i)=Sorting_mat_ICA(i)+comp_mat(i,k);
            end
        end
    end   
    
% Computing the worst rem_num ICs 
    var = strcat('Computing the worst ', int2str(rem_num), ' components...');
    disp(var)
    W=EEG.icaweights*EEG.icasphere;
    ics = W*EEG.data;
    [~,II]=sort(Sorting_mat_ICA);
    for temp_rem=1:rem_num
        if Sorting_mat_ICA(II((num_chan+1)-temp_rem))>0
            ics(II((num_chan+1)-temp_rem),:)=0;
        end
    end
    clean_data=inv(W)*ics;
end