% compute_GD_feat() - Computes Generic Discontinuity spatial feature 
%
% Usage:
%   >> res = compute_GD_feat(topografie,canali,num_componenti);
%
% Inputs:
%   topografie - topographies vector
%   canali     - EEG.chanlocs struct
%   num_componenti  - number of components
%
% Outputs:
%   res       - GDSF values



function res = compute_GD_feat(topografie,canali,num_componenti)

% Computes GDSF, discontinuity spatial feature
% topografie is the topography weights matrix
% canali is the structure EEG.chanlocs
% num_componenti is the number of ICs
% res is GDSF values

xpos=[canali.X];ypos=[canali.Y];zpos=[canali.Z];
pos=[xpos',ypos',zpos'];

res=zeros(1,num_componenti);   

for ic=1:num_componenti
    
    % consider the vector topografie(ic,:)
    
    aux=[];
    
    for el=1:length(canali)-1
        
        P=pos(el,:); %position of current electrode
        d=pos-repmat(P,length(canali),1);
        %d=pos-repmat(P,62,1);
        dist=sqrt(sum((d.*d),2));
        
        [y,I]=sort(dist);
        repchas=I(2:11); % list of 10 nearest channels to el
        weightchas=exp(-y(2:11)); % respective weights, computed wrt distance
        
        aux=[aux abs(topografie(ic,el)-mean(weightchas.*topografie(ic,repchas)'))];
                % difference between el and the average of 10 neighbors
                % weighted according to weightchas
    end
    
    res(ic)=max(aux);
    
end

