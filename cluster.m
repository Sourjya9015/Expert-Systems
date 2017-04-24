classdef cluster < hgsetget
    %CLUSTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %clusterId;
        numNodes;
        nodesPos;
        nodeEnergyUsage;
        channelLoss2Coord;
        channelLoss2AP;
        
        clusterCenter; % cluster center; default at origin
        clusterRadius; % cluster radius; default 10
        
        fc = 2.4*1e9; % operating frequency 2 GHz
        
        bw = 10*1e6; % 10 MHz of total BW
        Nf = 4; % noise factor in dB
        N0 = -174; % Noise psd in dB
        Ptx = 10; % Tx power in dBm
        nBitsTx = 800; % each node sends 100 bytes of data : fixed for now.
        gainRx = 12; % 12 dB of Rx gain at the AP
        gainTx = 4; % 4 dB of Tx gain at the nodes
        
    end
    
    methods
        function obj = cluster ( num, cntr, radius)
            
            %obj.clusterId = cId;
            
            obj.numNodes = num;
            obj.clusterCenter = cntr;
            obj.clusterRadius = radius;
            
            p1 = [rand(obj.numNodes,1) rand(obj.numNodes,1)];
            obj.nodesPos = p1*obj.clusterRadius + repmat(cntr, obj.numNodes,1);
            
            obj.nodeEnergyUsage = zeros(obj.numNodes,1);
        end
        
        % computes the channel loss to the specified accessPoint
        % works as a CQI feed back
        function computeChannelLoss (obj,apCoord)
            
            pmf = [0.1 0.55 0.35]; % LOS, NLOS and Blocking probability
            loss = [0 20 100]; % loss in dB due to channel state
            
            pdist = [0, cumsum(pmf)];
            u = rand(obj.numNodes,1);
            [~,indx] = histc(u,pdist);
            chLoss = loss(indx); % in dB
            
            dist = sqrt(sum((obj.nodesPos - repmat(apCoord, obj.numNodes,1)).^2,2));            
            
            lam = (3*1e8)/obj.fc;
            PL = 20.*log10(((4*pi*dist)./lam));
            
            obj.channelLoss2AP = PL + chLoss';
        end
        
        function computeLoss2Coordinator (obj, LeaderCoord)
            pmf = [0.4059 0.594 0.0001]; % very less blocking probabilty
            
            loss = [0 10 100]; % loss in dB due to channel state
            
            pdist = [0, cumsum(pmf)];
            u = rand(obj.numNodes,1);
            [~,indx] = histc(u,pdist);
            chLoss = loss(indx); % in dB
            
            dist = sqrt(sum((obj.nodesPos - repmat(LeaderCoord, obj.numNodes,1)).^2,2));
            
            lam = (3*1e8)/obj.fc;
            PL = 20.*log10(((4*pi*dist)./lam));
            
            PL(dist == 0) = 0;
            
            
            obj.channelLoss2Coord = PL + (dist ~= 0).*(chLoss');
        end
        
        % secondary nodes TX to coordinator. Coordinator Tx to AP/BS
        function transmit (obj)
            bwPerNode = obj.bw/(obj.numNodes - 1);
            SNR = obj.Ptx + obj.gainTx + obj.gainRx ...
                - obj.Nf - obj.N0 - 10*log10(obj.bw) - obj.channelLoss2Coord;
            SNR2AP = obj.Ptx + obj.gainTx + obj.gainRx ...
                - obj.Nf - obj.N0 - 10*log10(obj.bw) - obj.channelLoss2AP;

            % Obtaine rate from SNR using Shannon's law
            rate = (obj.channelLoss2Coord ~=0).*(bwPerNode.*log2(1+10.^(0.1*SNR))) + ...
                (obj.channelLoss2Coord == 0).*(obj.bw*log2(1+ 10.^(0.1*SNR2AP)));
            
            numBits = (obj.channelLoss2Coord ~=0).*(obj.nBitsTx*ones(obj.numNodes,1)) + ...
                (obj.channelLoss2Coord == 0).*(obj.numNodes*obj.nBitsTx*ones(obj.numNodes,1));
                    
            
            % TX time is num of bits/ rate
            txTime = numBits ./ rate;
            
            obj.nodeEnergyUsage = obj.nodeEnergyUsage + 10.^(0.1*obj.Ptx).*txTime*1e-3;
        end
        
        function flush (obj)
            obj.nodeEnergyUsage = zeros(obj.numNodes,1);
        end
    end
    
end

