%% Runs the simulation

clc; clear; close all;

numCluster = 8; % number clusters of machine nodes
xyAP = [0 0]; % coordinates of the access point

%clusterCenters = [50 50; 0 -50; -50 50];
radii = 100;
ang = (2*pi)*rand(numCluster,1);
p = radii*exp(-1i*ang);
clusterCenters = [real(p) imag(p)];



nNodes = 100;
radius = 20; % meters


for indx = 1:numCluster
    networkCluster(indx) = cluster (nNodes,clusterCenters(indx,:), radius);
end

nEpochs = 100;

AP = accessPoint ();

AP.set('location',xyAP, 'numClusters', numCluster);

AP.Initialize ();

options = {'no','fixed','variable'};

for opt = 1:length(options)
    
    AP.set('expertShare',cell2mat(options(opt)));
    for i=1:nEpochs

        cqiReport = containers.Map ();
        topology = containers.Map ();

        % Report the path loss to the AP/BS
        for indx = 1:numCluster
            key = char([99 48+indx]);
            networkCluster(indx).computeChannelLoss(xyAP);
            cqiReport(key) = networkCluster(indx).channelLoss2AP;
            topology(key) = networkCluster(indx).nodesPos;
        end


        AP.set('topology',topology, 'cqiFeedback', cqiReport);

        % BS has all different experts residing. Job of each expert is to find
        % leader coordinates with confidences using any alogirhtm
        xyLeaders = AP.selectCoordinators();  % Has 'average' distance measures for all choices of leaders for each cluster
        %[dis,leaderIndx]=min(Leaders);

        % BS returns the "leader coordinates"
        % Populate the leader co-ordinates in xyLeaders
        % this can also be a map of sorts.

        for indx = 1:numCluster        
            networkCluster(indx).computeLoss2Coordinator(xyLeaders(indx,:));
            networkCluster(indx).transmit();
        end

    end

    %% Gathering the energy used data
    energy = zeros(numCluster,3);

    for indx = 1:numCluster   
        arrEn = networkCluster(indx).nodeEnergyUsage;
        energy(indx,1) = mean(arrEn);
        energy(indx,2) = max(arrEn);
        energy(indx,3) = sqrt(mean((arrEn-  energy(indx,1)).^2));
        
        networkCluster(indx).flush();
    end

    energy = 10*log10(energy);

    figure(opt);
    str = sprintf('%s share update',cell2mat(options(opt)));
    bar(1:numCluster, energy);
    set(gca,'Fontsize',16);
    xlabel('Cluster'); ylabel('Power consumed (dB)');
    legend('Mean','Maximum','Std. Div','Location','SE');
    title(str);
    grid on;
    %ylim([-5 35]);
    
end

