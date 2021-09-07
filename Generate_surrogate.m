function surrogate_fmri = Generate_surrogate(fmri,num_randomisations)
%Generate_surrogate CREATES SURROGATE DATASETS FOR fMRI DATA
%
% Create surrogate datasets to use as null models for fMRI dynamic network
% analyses
%
% Input:
% fmri = original fmri 2D matrix of n*t where n is the brain regions and t
% the timepoints
%
% num_randomisations = number of surrogate datasets to be created
%
% Output:
% surrogate_fmri = cell array containing surrogate datasets in each cell

disp("Fourier transforming original dataset...")
fmritrans = fft(fmri,[],2); %transform timeseries
magnitude = abs(fmritrans); %save magnitude
surrogate_fmri = cell(num_randomisations,1);

for qq = 1:num_randomisations %Perform following steps for all required surrogates

    disp(strcat("  Creating surrogate dataset ",num2str(qq),"/",num2str(num_randomisations),"..."))
    
    %Randomise angle for each frequency; if even datapoints, leave freq bin
    %N/2 alone
    randangle = zeros(size(fmritrans));
    randangle(:,1) = angle(fmritrans(:,1));
    if mod(size(fmritrans), 2) == 0
        randangle(:,round(size(randangle,2)/2)+1)=angle(fmritrans(:,round(size(randangle,2)/2)+1)); 
    end
    for ii=2:round(size(randangle,2)/2)
        randnum = -pi + rand*2*pi;
        randangle(:,ii)=angle(fmritrans(:,ii))+randnum;
        randangle(:,end+2-ii) = angle(fmritrans(:,end+2-ii))-randnum;
    end

    %Inverse tranformation of randomized data and save in cell
    surrogate_fmri{qq} = ifft(magnitude.*exp(1i*randangle),[],2, 'symmetric');
    clearvars randangle
    
end
end