%% 1: Centrality in IO Network of Goods
clear; clc; close all; % Good practice
load IOdownload.mat; % Loading the data
n = 47; % The number of sectors in both economies
%% (a) In-Degree & Out-Degree Centrality
disp('In-Degree & Out-Degree Centrality');
disp(' ');

% Sweden
SE = io.swe2000; % Incidence matrix
GSE = digraph(SE); % Digraph

% In-degree centrality
inSE = centrality(GSE, 'indegree', 'Importance', GSE.Edges.Weight); % In-degree matrix
[~, tempIndex] = maxk(inSE, n);
disp('The three most central sectors in Sweden the year 2000 based on in-degree centrality:');
disp(['#1', name(tempIndex(1)), inSE(tempIndex(1))]);
disp(['#2', name(tempIndex(2)), inSE(tempIndex(2))]);
disp(['#3', name(tempIndex(3)), inSE(tempIndex(3))]);

% Out-degree centrality
outSE = centrality(GSE, 'outdegree', 'Importance', GSE.Edges.Weight); % Out-degree matrix
[~, tempIndex] = maxk(outSE, n);
disp('The three most central sectors in Sweden the year 2000 based on out-degree centrality:');
disp(['#1', name(tempIndex(1)), outSE(tempIndex(1))]);
disp(['#2', name(tempIndex(2)), outSE(tempIndex(2))]);
disp(['#3', name(tempIndex(3)), outSE(tempIndex(3))]);

% Indonesia
ID = io.idn2000; % Incidence matrix
GID = digraph(ID); % Digraph

% In-degree centrality
inID = centrality(GID, 'indegree', 'Importance', GID.Edges.Weight); % In-degree matrix
[~, tempIndex] = maxk(inID, n);
disp('The three most central sectors in Indonesia the year 2000 based on in-degree centrality:');
disp(['#1', name(tempIndex(1)), inID(tempIndex(1))]);
disp(['#2', name(tempIndex(2)), inID(tempIndex(2))]);
disp(['#3', name(tempIndex(3)), inID(tempIndex(3))]);

% Out-degree centrality
outID = centrality(GID, 'outdegree', 'Importance', GID.Edges.Weight); % Out-degree matrix
[~, tempIndex] = maxk(outID, n);
disp('The three most central sectors in Indonesia the year 2000 based on in-degree centrality:');
disp(['#1', name(tempIndex(1)), outID(tempIndex(1))]);
disp(['#2', name(tempIndex(2)), outID(tempIndex(2))]);
disp(['#3', name(tempIndex(3)), outID(tempIndex(3))]);
%% (b) Eigenvector Centrality on the Largest Connected Component
disp('------------------------------------------------------------------');
disp('Eigenvector Centrality on the Largest Connected Component');
disp(' ');

% Sweden
binsSE = conncomp(GSE); % 7 bins (6 isolated nodes)
[~, indexBinsSE] = find(binsSE == max(binsSE));
connGSE = subgraph(GSE, indexBinsSE);
connWSE = adjacency(connGSE, 'weighted');
[VSE, DSE] = eigs(connWSE);
[~, indSE] = max(diag(DSE));
xEigSE = VSE(:, indSE);
xEigSEMod = abs(xEigSE);
[~, tempIndex] = maxk(xEigSEMod, n);
disp('The three most central sectors in Sweden the year 2000 based on eigenvector centrality:');
disp(['#1', name(indexBinsSE(tempIndex(1))), xEigSEMod(tempIndex(1))]);
disp(['#2', name(indexBinsSE(tempIndex(2))), xEigSEMod(tempIndex(2))]);
disp(['#3', name(indexBinsSE(tempIndex(3))), xEigSEMod(tempIndex(3))]);

% Indonesia
binsID = conncomp(GID); % All belong to bin 1
[~, indexBinsID] = find(binsID == max(binsID)); % Size 47 -> already LCC
connGID = subgraph(GID, indexBinsID); % Same graph as GID
connWID = adjacency(connGID, 'weighted');
[VID, DID] = eigs(connWID);
[~, indID] = max(diag(DID));
xEigID = VID(:, indID);
xEigIDMod = abs(xEigID);
[~, tempIndex] = maxk(xEigIDMod, n);
disp('The three most central sectors in Indonesia the year 2000 based on eigenvector centrality:');
disp(['#1', name(indexBinsID(tempIndex(1))), xEigIDMod(tempIndex(1))]);
disp(['#2', name(indexBinsID(tempIndex(2))), xEigIDMod(tempIndex(2))]);
disp(['#3', name(indexBinsID(tempIndex(3))), xEigIDMod(tempIndex(3))]);
%% (c) Katz Centrality w/ Beta = .15 & 2 Different Values of Mu
disp('------------------------------------------------------------------');
disp('Katz Centrality w/ Beta = .15 & 2 Different Values of Mu');
disp(' ');

beta = .15;

% Case I: Mu = 1 for all nodes
muOnes = ones(length(name), 1);

% Case II: Mu = 1 for (31); 0 otherwise
muOne = zeros(length(name), 1);
muOne(31) = 1;

unitMatrix = eye(length(name));

% Sweden
WSE = adjacency(GSE, 'weighted'); % Weighted adjacency matrix
WSETransposed = WSE';
[VSEKatz, DSEKatz] = eigs(WSETransposed);
[~, indSEKatz] = max(diag(DSEKatz));
xEigSEKatz = VSEKatz(:, indSEKatz);
domEigSEKatz = max(abs(xEigSEKatz)); % The dominant eigen value
inverseSE = inv(unitMatrix - ((1 - beta) / domEigSEKatz) * WSETransposed);

% Case 1
katz1SE = inverseSE \ muOnes;
katzAbs1SE = abs(katz1SE);
[~, indKatzSE1] = maxk(katzAbs1SE, n);
disp('The three most central sectors in Sweden the year 2000 based on Katz centrality (case I):');
disp(['#1', name(indKatzSE1(1)), katzAbs1SE(indKatzSE1(1))]);
disp(['#2', name(indKatzSE1(2)), katzAbs1SE(indKatzSE1(2))]);
disp(['#3', name(indKatzSE1(3)), katzAbs1SE(indKatzSE1(3))]);

% Case 2
katz2SE = inverseSE \ muOne;
katzAbs2SE = abs(katz2SE);
[~, indKatzSE2] = maxk(katzAbs2SE, n);
disp('The three most central sectors in Sweden the year 2000 based on Katz centrality (case II):');
disp(['#1', name(indKatzSE2(1)), katzAbs2SE(indKatzSE2(1))]);
disp(['#2', name(indKatzSE2(2)), katzAbs2SE(indKatzSE2(2))]);
disp(['#3', name(indKatzSE2(3)), katzAbs2SE(indKatzSE2(3))]);

% Indonesia
WID = adjacency(GID, 'weighted'); % Weighted adjacency matrix
WIDTransposed = WID';
[VIDKatz, DIDKatz] = eigs(WIDTransposed);
[~, indIDKatz] = max(diag(DIDKatz));
xEigIDKatz = VIDKatz(:, indIDKatz);
domEigIDKatz = max(abs(xEigIDKatz)); % The dominant eigen value
inverseID = inv(unitMatrix - ((1 - beta) / domEigIDKatz) * WIDTransposed);

% Case 1
katz1ID = inverseID \ muOnes;
katzAbs1ID = abs(katz1ID);
[~, indKatzID1] = maxk(katzAbs1ID, n);
disp('The three most central sectors in Indonesia the year 2000 based on Katz centrality (case I):');
disp(['#1', name(indKatzID1(1)), katzAbs1ID(indKatzID1(1))]);
disp(['#2', name(indKatzID1(2)), katzAbs1ID(indKatzID1(2))]);
disp(['#3', name(indKatzID1(3)), katzAbs1ID(indKatzID1(3))]);

% Case 2
katz2ID = inverseID \ muOne;
katzAbs2ID = abs(katz2ID);
[~, indKatzID2] = maxk(katzAbs2ID, n);
disp('The three most central sectors in Indonesia the year 2000 based on Katz centrality (case II):');
disp(['#1', name(indKatzID2(1)), katzAbs2ID(indKatzID2(1))]);
disp(['#2', name(indKatzID2(2)), katzAbs2ID(indKatzID2(2))]);
disp(['#3', name(indKatzID2(3)), katzAbs2ID(indKatzID2(3))]);
%% 2: Influence on Twitter
load -ascii users.mat; % Loading the users

% Loading the adjacency matrix
load -ascii twitter.mat;
G = digraph(twitter(:, 1), twitter(:, 2), twitter(:, 3));
W = adjacency(G);
%% (a) PageRank w/ Beta = .15 & Mu = 1 for All Nodes
disp('------------------------------------------------------------------');
disp('PageRank w/ Beta = .5 & Mu = 1 for All Nodes');
disp(' ');

beta = .15;
mu = ones(length(W), 1);

% Calculate normalized adjacency matrix
w = sum(W, 2); % Out-degree
W(diag(w == 0)) = 1; % Self-loops
w(~w) = 1; % Adjusting the out-degree
P = diag(w) \ W; % Final calculation

pi = rand(length(W), 1); % Initial value (assigned randomly)

% The iterative approach
for i = 1:1e2
    pi = (1 - beta) * P' * pi + beta * mu;
end

% The five most central nodes (in order)
[~, IDs] = sort(pi, 'descend');
topFiveIDs = IDs(1:5);
disp('The top five most central user IDs according to PageRank: ');
tweeters = users(topFiveIDs);
disp(['1: ', num2str(tweeters(1))]);
disp(['2: ', num2str(tweeters(2))]);
disp(['3: ', num2str(tweeters(3))]);
disp(['4: ', num2str(tweeters(4))]);
disp(['5: ', num2str(tweeters(5))]);
%% (b) Change in Opinion w/ 2 Stubborn Nodes
% Selection of stubborn nodes
stubbornOne = 1999;
stubbornTwo = numnodes(G) - stubbornOne;

% Initialize actual stubbornness
stubbornG = rmedge(G, stubbornOne, 1:numnodes(G));
stubbornG = rmedge(stubbornG, stubbornTwo, 1:numnodes(G));
stubbornW = adjacency(stubbornG);

% Observing a handful of nodes
obsNodes(1) = stubbornOne;
obsNodes(2) = stubbornTwo;
obsNodes(3) = 1000;
obsNodes(4) = 3000;
obsNodes(5) = 6000;
sort(obsNodes);

% Calculate normalized adjacency matrix
wStubborn = sum(stubbornW, 2); % Out-degree
stubbornW(diag(wStubborn == 0)) = 1; % Self-loops
wStubborn(~wStubborn) = 1; % Adjusting the out-degree
stubbornP = diag(wStubborn) \ stubbornW; % Final calculation

% Initial opinions
z = rand(length(stubbornP), 1); % All values = .50 -> no change (for these nodes)
z(stubbornOne) = 1;
z(stubbornTwo) = 0;

iterations = 3e2;
savedObs = zeros(length(obsNodes), iterations);

for i = 1:iterations
    if i ~= 1
        z = stubbornP * z;
    end
    
    for j = 1:length(obsNodes)
        savedObs(j, i) = z(obsNodes(j));
    end
end

figure(); % Figure 1

for j = 1:length(obsNodes)
    plot(savedObs(j, :));

    if j ~= length(obsNodes)
        hold on;
    end
end

tweetersNew = users(obsNodes);
legend(num2str(tweetersNew(1)), num2str(tweetersNew(2)), ...
    num2str(tweetersNew(3)), num2str(tweetersNew(4)), num2str(tweetersNew(5)));
title('How the Opinions Change Over Time for Different Users');
xlabel('Iteration');
ylabel('Opinion');
hold off;
%% (c1) Opinion Distribution w/ Stubbornness Taken from 1st & 2nd PageRank
% Selection of stubborn nodes
stubbornOne = topFiveIDs(1);
stubbornTwo = topFiveIDs(2);

% Initialize actual stubbornness
stubbornG = rmedge(G, stubbornOne, 1:numnodes(G));
stubbornG = rmedge(stubbornG, stubbornTwo, 1:numnodes(G));
stubbornW = adjacency(stubbornG);

% Calculate normalized adjacency matrix
wStubborn = sum(stubbornW, 2); % Out-degree
stubbornW(diag(wStubborn == 0)) = 1; % Self-loops
wStubborn(~wStubborn) = 1; % Adjusting the out-degree
stubbornP = diag(wStubborn) \ stubbornW; % Final calculation

% Initial opinions
z = ones(length(stubbornP), 1) * .50;
z(stubbornOne) = 1;
z(stubbornTwo) = 0;

iterations = 1e4;

for i = 1:iterations
    z = stubbornP * z;
end

figure(); % Figure 2
histogram(z);
title('Opinion Distribution Based on PageRank Critera (1st = 1, 2nd = 0)');
xlabel('Opinion');
ylabel('Frequency');
%% (cII) Opinion Distribution w/ Stubbornness Taken from 1st & 5th PageRank
% Selection of stubborn nodes
stubbornOne = topFiveIDs(1);
stubbornTwo = topFiveIDs(5);

% Initialize actual stubbornness
stubbornG = rmedge(G, stubbornOne, 1:numnodes(G));
stubbornG = rmedge(stubbornG, stubbornTwo, 1:numnodes(G));
stubbornW = adjacency(stubbornG);

% Calculate normalized adjacency matrix
wStubborn = sum(stubbornW, 2); % Out-degree
stubbornW(diag(wStubborn == 0)) = 1; % Self-loops
wStubborn(~wStubborn) = 1; % Adjusting the out-degree
stubbornP = diag(wStubborn) \ stubbornW; % Final calculation

% Initial opinions
z = ones(length(stubbornP), 1) * .50;
z(stubbornOne) = 1;
z(stubbornTwo) = 0;

for i = 1:iterations
    z = stubbornP * z;
end

figure();
histogram(z);
title('Opinion Distribution Based on PageRank Critera (1st = 1, 5th = 0)');
xlabel('Opinion');
ylabel('Frequency');
%% (cIII) Opinion Distribution w/ Stubbornness Taken from 3rd & 4th PageRank
% Selection of stubborn nodes
stubbornOne = topFiveIDs(3);
stubbornTwo = topFiveIDs(4);

% Initialize actual stubbornness
stubbornG = rmedge(G, stubbornOne, 1:numnodes(G));
stubbornG = rmedge(stubbornG, stubbornTwo, 1:numnodes(G));
stubbornW = adjacency(stubbornG);

% Observing a handful of nodes
obsNodes = randi(length(stubbornW), 5);
obsNodes(1) = stubbornOne;
obsNodes(2) = stubbornTwo;
sort(obsNodes);

% Calculate normalized adjacency matrix
wStubborn = sum(stubbornW, 2); % Out-degree
stubbornW(diag(wStubborn == 0)) = 1; % Self-loops
wStubborn(~wStubborn) = 1; % Adjusting the out-degree
stubbornP = diag(wStubborn) \ stubbornW; % Final calculation

% Initial opinions
z = ones(length(stubbornP), 1) * .50;
z(stubbornOne) = 1;
z(stubbornTwo) = 0;

for i = 1:iterations
    z = stubbornP * z;
end

figure(); % Figure 4
histogram(z);
title('Opinion Distribution Based on PageRank Critera (3rd = 1, 4th = 0)');
xlabel('Opinion');
ylabel('Frequency');