function [ pred ] = expert( APCoord, topology, numC, selOpt,opt)
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
            nodePos = nodePos + randn(size(nodePos));
            
            len = size( nodePos , 1); 
            
            dist = sum((repmat(APCoord,len,1) - nodePos).^2,2);
            
            prob = 1./dist;
            
            if ( strcmp(opt, 'mad') == 1)
                pmf = rand(1,len);
                pred(key) = pmf./sum(pmf);
            else
                pred(key) = (prob./sum(prob))';
            end
        end
        
    case 'closestClus'
        % Output distribution based on closeness to every other node
        for indx=1:numC
            key = char([99 48+indx]);
            nodePos = topology(key); nodePos = nodePos + randn(size(nodePos));
            nnode = size( nodePos , 1); 
            pmf = zeros(1,nnode);
            
            for cnt=1:nnode
                nodeCoord = nodePos(cnt,:);
                
                dist = sqrt(sum((repmat(nodeCoord,nnode,1) - nodePos).^2,2));
                pmf(cnt) = 1/mean(dist);
            end
            
            if ( strcmp(opt, 'mad') == 1)
                pmf = abs(randn(1,nnode));
                pred(key) = pmf./sum(pmf);
            else    
                pred(key) = pmf./sum(pmf);
            end
            
        end
        
    case 'uniform'
        % Output distribution is uniform, i.e., all nodes are equally
        % likely
        for indx=1:numC
            key = char([99 48+indx]);
            nodePos = topology(key); nodePos = nodePos + randn(size(nodePos));
            len = size( nodePos , 1); 
          
            pred(key) = (1/len)*ones(1,len);
        end
        
        
    otherwise
        error('Unknown prediction option specified');       
end



end

