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


options = {'no','fixed','variable'};
%options = {'no'};

for opt = 1:length(options)
    
    weightsExprts = zeros(AP.numExperts, nEpochs+1);
    
    AP.Initialize ();
    
    weightsExprts(:,1) = AP.expertWt;
    
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
        
        weightsExprts(:,i+1) = AP.expertWt;

        for indx = 1:numCluster        
            networkCluster(indx).computeLoss2Coordinator(xyLeaders(indx,:));
            networkCluster(indx).transmit();
        end

    end

    %% Gathering the energy used data
    energy = zeros(numCluster,2);
    stdDiv = zeros(numCluster,1);
    for indx = 1:numCluster   
        arrEn = networkCluster(indx).nodeEnergyUsage;
        energy(indx,1) = mean(arrEn);
        energy(indx,2) = max(arrEn);
        stdDiv(indx) = sqrt(mean((arrEn-  energy(indx,1)).^2));
        
        networkCluster(indx).flush();
    end

    energy = 10*log10(energy) + 30 ;
    stdDiv = 10*log10(stdDiv) + 30 ;

    figure(opt);
    str = sprintf('%s share update',cell2mat(options(opt)));
    subplot(2,1,1);bar(1:numCluster, energy);
    set(gca,'Fontsize',12);
    xlabel('Cluster'); ylabel('Power consumed (dBm)');
    legend('Mean','Maximum','Location','SE');
    title(str);
    subplot(2,1,2);errorbar(1:numCluster, energy(:,1),stdDiv);
    set(gca,'Fontsize',12);
    xlabel('Cluster'); ylabel('Power consumed (dBm)');
    title('Mean energy with standard div.');
    grid on;
    %ylim([-5 35]);
    
    figure(4+opt);
    plot(0:nEpochs, weightsExprts(1,:),'-', 0:nEpochs, weightsExprts(2,:),'-', ...
        0:nEpochs, weightsExprts(3,:),'-', 'Linewidth',2);
    set(gca,'Fontsize',16);
    xlabel('Num of epochs'); ylabel('Weights');
    legend('Expert 1','Expert 2', 'Expert 3');
    title(str);
    grid on;
end

