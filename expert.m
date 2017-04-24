function [ pred ] = expert( APCoord, topology, numC, selOpt )
%EXPERT1:Expert prediction for the first expert
% Input:
%        APCoord: Co-ordinate of the access  point
%        topology: Topology of the system
%        numC : Number of clusters
%        selOpt : Selection option
% Output:
%        pred: Expert prediction, a Map

pred = containers.Map ();

switch (selOpt)
    case 'closestAP' 
        % Output distribution based on closeness to AP
        for indx=1:numC
            key = char([99 48+indx]);
            nodePos = topology(key);
            
            len = size( nodePos , 1); 
            
            dist = sum((repmat(APCoord,len,1) - nodePos).^2,2);
            
            pred(key) = (dist./sum(dist))';
        end
        
    case 'closestClus'
        % Output distribution based on closeness to every other node
        for indx=1:numC
            key = char([99 48+indx]);
            nodePos = topology(key);
            nnode = size( nodePos , 1); 
            pmf = zeros(1,nnode);
            
            for cnt=1:nnode
                nodeCoord = nodePos(cnt,:);
                
                dist = sqrt(sum((repmat(nodeCoord,nnode,1) - nodePos).^2,2));
                pmf(cnt) = mean(dist);
            end
            
            pred(key) = pmf./sum(pmf);
            
        end
        
    case 'uniform'
        % Output distribution is uniform, i.e., all nodes are equally
        % likely
        for indx=1:numC
            key = char([99 48+indx]);
            nodePos = topology(key);
            len = size( nodePos , 1); 
          
            pred(key) = (1/len)*ones(1,len);
        end
        
        
    otherwise
        error('Unknown prediction option specified');       
end



end

