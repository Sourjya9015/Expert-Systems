function [ L ] = loss( pmf, leaderPos, nodePos )
% LOSS: Computes the loss on one expert for one cluster
% Input:
%       pmf : the distribution predicted by the expert for cluster 1
%       leaderPos : Predicted position of the leader
%       nodePos : network topology for a given cluster

len = size(nodePos,1);

dist = sum((repmat(leaderPos,len,1) - nodePos).^2,2);

L = sum(pmf.* dist');

end

