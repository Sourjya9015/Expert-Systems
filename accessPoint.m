classdef accessPoint < hgsetget
    %ACCESSPOINT: Defines the accesspoint or the base station
    % Needs to select the co-ordinators based on expert advice
    
    properties
        location;
        numClusters;
        topology; % stores the network topology = coordinates of every node
                  %  declare this as a Map with the key as the cluster index c1, c2,
                  %  c3 etc.
        
        cqiFeedback; % channel quality feedback; again a Map set from the main simulation
        
        % Expert setting
        numExperts = 3; % Number of experts. Let's keep this fixed for now.
        expertWt;   % Weights on each expert
        expertShare; % 'static','fixed', or 'variable'
        
        eta = 0.1; % can set eta using set
        alpha = 0.4; % set a value for alpha
        
    end
    
    methods
        
        function Initialize (obj)
            % NclusX1 vector with weights on each experts
            obj.expertWt = (1/obj.numExperts)*ones(obj.numExperts,1);
            %obj.expertWt = [1; zeros(obj.numExperts,1)];
        end
        
        function xyLeaders = selectCoordinators(obj)
            % I do not think we should do it like this. Where is the
            % part on expert advice?? This is just centroid calculation.
            
            % Get expert predictions. These are all maps with key c1, c2
            % etc.
            
            xyLeaders = zeros(obj.numClusters, 2);
            
            e1 = expert( obj.location, obj.topology, obj.numClusters, 'closestAP' );
            e2 = expert( obj.location, obj.topology, obj.numClusters, 'closestClus' );
            e3 = expert( obj.location, obj.topology, obj.numClusters, 'uniform' );
            
            Losses = zeros(obj.numExperts,1);
            
            for ind = 1:obj.numClusters
                
                key = char([99 48+ind]);
                pos = obj.topology(key);
                
                nnodes = size(pos,1);
                
                e1x = e1(key);
                e2x = e2(key);
                e3x = e3(key);
                
                cqi = obj.cqiFeedback(key);
                
                % Prediction
                pmf = sum(repmat(obj.expertWt,1,nnodes).*([e1x; e2x; e3x]));
                pmf = pmf.*cqi'; % weigh the prediction with the channel quality
                pmf = pmf./sum(pmf);
                [~,index] = max(pmf);

                xyLeaders(ind,:) = pos(index,:);
                
                % Compute losses
                Losses(1) = Losses(1) + (1/obj.numClusters)*loss( e1x, xyLeaders(ind,:), pos );
                Losses(2) = Losses(2) + (1/obj.numClusters)*loss( e2x, xyLeaders(ind,:), pos );
                Losses(3) = Losses(3) + (1/obj.numClusters)*loss( e3x, xyLeaders(ind,:), pos );
                              
            end
            
            
            obj.expertWt = obj.expertWt.* exp(-obj.eta*Losses);
            
            switch(obj.expertShare)
                case 'no'
                    % no changes required
                case 'fixed'
                    pool = sum(obj.alpha*obj.expertWt);
                    obj.expertWt = (1-obj.alpha)*obj.expertWt + ...
                        1/(obj.numClusters-1)*(pool - obj.alpha*obj.expertWt);
                case 'variable'
                    pool = sum((1 - (1-obj.alpha).^Losses).*obj.expertWt);
                    obj.expertWt = ((1-obj.alpha).^Losses).* obj.expertWt + ...
                        1/(obj.numClusters-1)*(pool - (1 - (1-obj.alpha).^Losses).*obj.expertWt);
                    
                otherwise
                    error('Invalid update method');
            end

        end
    end
    
end

