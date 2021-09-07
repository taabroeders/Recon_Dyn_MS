function [netwassign_final,varargout] = CommunityDetection(timeseries,initial_assignment,window_size,step_size,varargin)
%INPUT:
% timeseries= t*n matrix containing the subjects' timeseries. t is the
% number of timepoints and n the number of regions.
% initial_assignment= n vector containing the initial assignment (e.g.
% literature-based)
% window_size= integer indicating the window-size
% step_size= integer indicating the step-size
%optional:
% gausswin= use gaussian windows (true/false). Need weightedcorrs.m in
% path (www.mathworks.com/matlabcentral/fileexchange/20846-weighted-correlation-matrix)
% maxit= maximum number of iterations
%
%OUTPUT:
% netwassign_final= n*w matrix containing the windowed network assignments.
% n is the number of regions and w is the number of windows.
%optional:
% itw= the number of iterations until convergence for each window

%Variable input arguments
if nargin == 4
    gausswin=false;
    maxit=1000;
elseif nargin == 5
    gausswin=varargin{1};
    maxit=1000;
elseif nargin == 6
    gausswin=varargin{1};
    maxit=varargin{2};
else
    error('incorrect number of input variables')
end

%Check variable sizes
if length(initial_assignment)~=size(timeseries,2)
    error('Initial assignment variable does not have the same number of regions as the timeseries variable')
elseif length(window_size)~=1; error('Window size: need single integer value')
elseif length(step_size)~=1; error('Step size: need single integer value')
elseif length(gausswin)~=1; error('Gaussian window: need single logical value')
elseif length(maxit)~=1; error('Maximum number of iterations: need single integer value')    
elseif size(timeseries,1)<size(timeseries,2)
    warning('Warning: more regions in timeseries than timepoints, could indicate incorrect orientation (t*n)')
end

%Check variable datatypes
if ~isnumeric(timeseries); error('Timeseries: wrong datatype entered')
elseif ~isnumeric(initial_assignment) || ~all(mod(initial_assignment,1)==0); error('Initial assignment: wrong datatype entered')
elseif ~isnumeric(window_size) || ~mod(window_size,1)==0; error('Window size: need single integer value')
elseif ~isnumeric(step_size) || ~mod(step_size,1)==0; error('Step size: need single integer value')
elseif ~islogical(gausswin); error('Gaussian window: need single logical value');
elseif ~isnumeric(maxit) || ~mod(maxit,1)==0; error('Maximum number of iterations: need single integer value')
end

%Set gaussian window parameters if needed
if gausswin
    window_shape=gausswin(window_size,(window_size-1)/6);window_shape=window_shape/mean(window_shape); % stdev = (N-1)/(2*alpha), so stdev=3
end

%Count number of networks/regions/windows
num_netw=length(unique(initial_assignment));
numreg=length(initial_assignment);
numwin=length(1:step_size:(size(timeseries,1)-(window_size-1)));
    
qq=0;itw=zeros(numwin,1);netwassign_final=zeros(numreg,numwin);
for jj = 1:step_size:(size(timeseries,1)-(window_size-1))    %Create windows
    qq=qq+1;
    fmri_epoch=timeseries(jj:jj+(window_size-1),:);
    
    %Calculate connectivity
    if gausswin
        connect = abs(atanh(weightedcorrs(fmri_epoch,window_shape)));
    else
        connect = abs(atanh(corr(fmri_epoch)));
    end
    connect(logical(eye(size(connect)))) = NaN;
    
    costmat = zeros(numreg,1);
    netwassign=initial_assignment;    % first set assignment to initial assignment
    tmporderprev=0;    %keep track of worst-assigned region for each previous iteration
    
    % iteratively update the assignment until convergence
    for it=1:maxit
        
        % calculate the quality of each current assignment for all regions
        for ii = 1:numreg
            winode = mean(connect(ii,netwassign == netwassign(ii)),'omitnan');
            btnode = mean(connect(ii,netwassign ~= netwassign(ii)),'omitnan');
            costmat(ii) = (winode-btnode)/(winode+btnode);
        end
        
        % select the region that shows the worst assignment quality
        [~,tmporder]=sort(costmat);
        
        % stop the algorithm if the assignment is the same as the previous iteration (i.e. convergence)
        if tmporderprev~=tmporder(1)
            q=tmporder(1);
        else
            break
        end
        tmporderprev=q; % save worst assigned region for the next iteration
        
        % calculate the network with the best match for the worst region and reassign
        newcost = zeros(num_netw,1);
        for ii=1:num_netw
            newcost(ii) = mean(connect(q,netwassign==ii),'omitnan');
        end
        [~,z] = max(newcost); % determine best match
        netwassign(q)=z; % update assignment
        
    end
    itw(qq)=it;
    netwassign_final(:,qq)=netwassign; %save final assignment for each window
end
    
varargout{1}=itw;

end

