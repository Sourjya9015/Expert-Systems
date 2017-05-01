%% Runs the simulation

clc; clear; close all;

numCluster = 8; % number clusters of machine nodes
xyAP = [0 0]; % coordinates of the access point

%clusterCenters = [50 50; 0 -50; -50 50];
radii = 100*ones(numCluster,1) + 20*randn(numCluster,1);
ang = linspace(0,2*pi-0.5,numCluster);%(2*pi)*rand(numCluster,1);
p = radii.*exp(-1i*ang');
clusterCenters = [real(p) imag(p)];
%clusterCenters = [100*randn(numCluster,1) 100*randn(numCluster,1)];


nNodes = 20;
radius = 40; % meters
learnRate = 0.1;
share = 0.1;


for indx = 1:numCluster
    networkCluster(indx) = cluster (nNodes,clusterCenters(indx,:), radius);
end

figure(11);
plot(networkCluster(1).nodesPos(:,1), networkCluster(1).nodesPos(:,2), 'xr', ...
    networkCluster(2).nodesPos(:,1), networkCluster(2).nodesPos(:,2), 'og',...
    networkCluster(3).nodesPos(:,1), networkCluster(3).nodesPos(:,2), 'sb', ...
    networkCluster(4).nodesPos(:,1), networkCluster(4).nodesPos(:,2), 'xk',...
    networkCluster(5).nodesPos(:,1), networkCluster(5).nodesPos(:,2), 'oc',...
    networkCluster(6).nodesPos(:,1), networkCluster(6).nodesPos(:,2), 'sm',...
    networkCluster(7).nodesPos(:,1), networkCluster(7).nodesPos(:,2), 'xb',...
    networkCluster(8).nodesPos(:,1), networkCluster(8).nodesPos(:,2), 'ok');
grid on;

nEpochs = 200;

AP = accessPoint ();

AP.set('location',xyAP, 'numClusters', numCluster,'eta',learnRate,'alpha', share, ...
    'expertType','sane');


options = {'no','fixed','variable'};
%options = {'no'};

losses = zeros(length(options),nEpochs);

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
        [xyLeaders, lossIter] = AP.selectCoordinators();  % Has 'average' distance measures for all choices of leaders for each cluster
        
        losses(opt,i) = sum(lossIter);
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

     energy = log(energy*1e6) ;
     stdDiv = log(stdDiv*1e6);

    figure(opt);
    str = sprintf('%s share update',cell2mat(options(opt)));
    subplot(2,1,1);bar(1:numCluster, energy);
    set(gca,'Fontsize',12);
    xlabel('Cluster'); ylabel('Energy (log micro J)');
    legend('Mean','Maximum','Location','SE');
    title(str);
    subplot(2,1,2);errorbar(1:numCluster, energy(:,1),stdDiv);
    set(gca,'Fontsize',12);
    xlabel('Cluster'); ylabel('Energy (log micro J)');
    title('Mean energy with standard div.');
    grid on;

%     figure(opt)
%     str = sprintf('%s share update',cell2mat(options(opt)));
%     semilogy(1:numCluster, energy(:,1), '-s', 1:numCluster, energy (:,2), '-s','Linewidth',2);
%     set(gca,'Fontsize',12);
%     xlabel('Cluster'); ylabel('Energy (Joules)');
%     legend('Mean','Maximum','Location','SE'); title(str); grid on;
    
    figure(4+opt);
    plot(0:nEpochs, weightsExprts(1,:),'-', 0:nEpochs, weightsExprts(2,:),'-', ...
        0:nEpochs, weightsExprts(3,:),'-', 'Linewidth',2);
    set(gca,'Fontsize',16);
    xlabel('Num of epochs'); ylabel('Weights');
    legend('Expert 1','Expert 2', 'Expert 3');
    title(str);
    grid on;
end

% figure(10);
% plot(1:nEpochs,losses(1,:),'-', 1:nEpochs,losses(2,:),'-',1:nEpochs,losses(3,:),'-', 'Linewidth', 2);
% set(gca,'Fontsize',16);
% xlabel('Iterations'); ylabel('\sum L(t,i)');
% grid on;
