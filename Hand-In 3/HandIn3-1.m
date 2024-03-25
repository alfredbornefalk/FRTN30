%% 1: Single-Particle Random Walk
clear; clc; close all; % Good practice
nodes = ['o', 'a', 'b', 'c', 'd']; % Defining the nodes

% Transition matrix (o, a, b, c, d)
W = [0, 2/5, 1/5, 0, 0;
     0, 0, 3/4, 1/4, 0;
     1/2, 0, 0, 1/2, 0;
     0, 0, 1/3, 0, 2/3;
     0, 1/3, 0, 1/3, 0];

w = sum(W, 2); % Out-degree
D = diag(1 - w); % Diagonal matrix
Q = D + W; % Normalized transition matrix
%% (a) The Simulated Average Return Time
disp('The Return Time');
disp(' ');

tickMax = 1e6; % 1,000,000 steps
nbrOfParticles = 1;
node = find(nodes == 'a'); % Node index

% Simulation process
ticks = cumsum(exprnd(1 / nbrOfParticles, tickMax, 1)); % Generate times
particlesPosition = node * ones(1, nbrOfParticles);
particlesToMove = randi(nbrOfParticles, 1, tickMax);
movements = zeros(1, tickMax); % Preallocation

% One particle is moved at a time
for i = 1:tickMax
    currentParticle = particlesToMove(i);
    currentNode = particlesPosition(currentParticle);
    nextNode = find(cumsum(Q(currentNode, :)) > rand, 1); % Randomize
    particlesPosition(currentParticle) = nextNode; % Update
    movements(i) = nextNode; % Save
end

moves = [particlesToMove; movements];

% Calculating the simulated average return time
simulatedAverageReturnTime = mean(diff(ticks(moves(2, :) == node)));

disp(['The simulated average return time of node a is ', ...
    num2str(simulatedAverageReturnTime), '.']);
%% (b) The Theoretical Average Return Time
% Calculating the normalized transition matrix
QTemp = Q;
qTemp = sum(Q, 2);
QTemp(diag(qTemp == 0)) = 1; % Self-loops
qTemp(~qTemp) = 1; % Adjusting the out-degree
P = diag(qTemp) \ QTemp;

[pi, ~] = eigs(P', 1); % Eigenvector corresponding to eigenvalue 1

if sum(pi) < 0
    pi = -pi;
end

pi = pi / norm(pi, 1);
trueAverageReturnTime = 1 / pi(node);

disp(['The theoretical average return time of node a is ', ...
    num2str(trueAverageReturnTime), '.']);
%% (c) The Simulated Average Hitting Time
disp('------------------------------------------------------------------');
disp('The Hitting Time');
disp(' ');

startNode = find(nodes == 'o'); % Start node index
endNode = find(nodes == 'd'); % End node index

% Simulation process
ticks = cumsum(exprnd(1 / nbrOfParticles, tickMax, 1)); % Generate times
particlesPosition = startNode * ones(1, nbrOfParticles);
particlesToMove = randi(nbrOfParticles, 1, tickMax);
movements = zeros(1, tickMax); % Preallocation

% One particle is moved at a time
for i = 1:tickMax
    currentParticle = particlesToMove(i);
    currentNode = particlesPosition(currentParticle);
    nextNode = find(cumsum(Q(currentNode, :)) > rand, 1); % Randomize
    particlesPosition(currentParticle) = nextNode; % Update
    movements(i) = nextNode; % Save
end

moves = [particlesToMove; movements];

% Find every hit and save the clock times
start = false;
result = []; % Size unknown -> no preallocation

for j = 1:length(movements)
    if movements(j) == startNode && ~start
        start = true;
        count = j;
    end

    if movements(j) == endNode && start
        result = [result, ticks(j) - ticks(count)];
        start = false;
    end
end

simulatedAverageHittingTime = mean(result); % Simulated average
disp(['The simulated average hit time from node o to d is ', ...
    num2str(simulatedAverageHittingTime), '.']);
%% (d) The Theoretical Average Hitting Time
subQ = Q;
subQ(endNode, :) = [];
subQ(:, endNode) = [];
trueAverageHittingTime = (eye(length(subQ)) - subQ) \ ones(length(subQ), 1);
disp(['The theoretical average hit time from node o to d is ', ...
    num2str(trueAverageHittingTime(startNode)), '.']);
%% 2: Graph Coloring & Network Games
clear; % Clearing the workspace
%% (a) Line Graph With Ten Nodes
disp('------------------------------------------------------------------');
disp('Line Graph With Ten Nodes');
disp(' ');

G = graph(1:9, 2:10); % The line graph
redAndGreen = zeros(1, 10); % 0 for red, 1 for green
maximumIterations = 5e2; % The (maximum) # iterations (changed from 1e2)
speed = .01; % Given in seconds

figure();

p = plot(G);
currentRed = find(redAndGreen == 0);
currentGreen = find(redAndGreen == 1);
highlight(p, currentGreen, 'NodeColor', 'green');
highlight(p, currentRed, 'NodeColor', 'red');

U = zeros(1, maximumIterations + 1);
U(1) = 9;

% Iterate and update the states
for time = 1:maximumIterations
    % Animation: setup
    p = plot(G);
    currentRed = find(redAndGreen == 0);
    currentGreen = find(redAndGreen == 1);
    highlight(p, currentGreen, 'NodeColor', 'green');
    highlight(p, currentRed, 'NodeColor', 'red');

    % Animation: pausing
    pause(speed);

    % Node selection (uniformly drawn)
    I = randi(10);

    % Getting the neighbours
    neighborNodes = neighbors(G, I);

    % Pre-calculations
    lSum = sum(redAndGreen(neighborNodes) == 0); % Sum of red
    fSum = sum(redAndGreen(neighborNodes) == 1); % Sum of green
    eta = time * .01;

    % Probability of getting the green state (1)
    P = exp(-eta * fSum) / (exp(-eta * fSum) + exp(-eta * lSum));

    % Random selection of new state, and updating the state vector
    redAndGreen(I) = rand <= P;

    % Estimation of the potential
    iSum = zeros(10, 1);
    
    for i = 1:10
        iSum(i) = sum(redAndGreen(i) == redAndGreen(neighbors(G, i)));
    end
    
    U(time + 1) = .5 * sum(iSum);
end

% Final animation
p = plot(G);
currentRed = find(redAndGreen == 0);
currentGreen = find(redAndGreen == 1);
highlight(p, currentGreen, 'NodeColor', 'green');
highlight(p, currentRed, 'NodeColor', 'red');

% Plot the potential function
plot(linspace(0, maximumIterations, length(U)), U);

disp(['The potential is ', num2str(U(maximumIterations + 1))]);
disp(' ');
%% (b1) Assigning WiFi-Channels to Routers w/ eta_3
disp('------------------------------------------------------------------');
disp('Assigning WiFi-Channels to Routers w/ eta_3');
disp(' ');

clear; % Again, clearing the workspace is necessary

load -ASCII wifi.mat; % Adjacency matrix
load -ASCII coord.mat; % The routers' coordinates

G = graph(wifi); % The Wi-Fi graph
routers = ones(100, 1); % The state vector
maximumIterations = 1e3; % The (maximum) number of iterations
eta = (1:maximumIterations) * 1e-2; % Dividing w/ higher (multiplying w/ lower) -> higher potential

U = zeros(1, maximumIterations + 1);

% Calculation of the total costs of the network
for i = 1:100
    diff = abs(routers(i) - routers(neighbors(G, i)));
    v = zeros(size(diff));
    v(diff == 0) = 2;
    v(diff == 1) = 1;
    iSum(i) = sum(v);
end

% Estimation of the potential
U(1) = .5 * sum(iSum);

% Iteration process
for time = 1:maximumIterations
    I = randi(100); % Randomly select a node
    neigh = neighbors(G, I); % Extract its neighbors
    lSum = zeros(1, 8); % Pre-calculate sums

    % The cost function
    for j = 1:8
        diff = abs(routers(neigh) - j);
        v = zeros(size(diff));
        v(diff == 0) = 2;
        v(diff == 1) = 1;
        lSum(j) = sum(v);
    end

    P = exp(-eta(time) * lSum) / sum(exp(-eta(time) * lSum)); % Probability calculation
    routers(I) = find(cumsum(P) > rand, 1); % Randomly select a new state for node I

    iSum = zeros(100, 1);

    % Calculation of the total costs of the network
    for i = 1:100
        diff = abs(routers(i) - routers(neighbors(G, i)));
        v = zeros(size(diff));
        v(diff == 0) = 2;
        v(diff == 1) = 1;
        iSum(i) = sum(v);
    end
    
    % Estimation of the potential
    U(time + 1) = .5 * sum(iSum);
end

% Plot the assigned Wi-Fi channels
figure();
p = plot(G);
highlight(p, find(routers == 1), 'NodeColor', 'red');
highlight(p, find(routers == 2), 'NodeColor', 'green');
highlight(p, find(routers == 3), 'NodeColor', 'blue');
highlight(p, find(routers == 4), 'NodeColor', 'yellow');
highlight(p, find(routers == 5), 'NodeColor', 'magenta');
highlight(p, find(routers == 6), 'NodeColor', 'cyan');
highlight(p, find(routers == 7), 'NodeColor', 'white');
highlight(p, find(routers == 8), 'NodeColor', 'black');

% Plot the potential function
plot(linspace(0, maximumIterations, length(U)), U);
hold on;

disp(['The potential is ', num2str(U(maximumIterations + 1))]);
disp(' ');
%% (b2) Assigning WiFi-Channels to Routers w/ all etas
disp('------------------------------------------------------------------');
disp('Assigning WiFi-Channels to Routers w/ all etas');
disp(' ');

eta = [(1:maximumIterations) * 1e-4, (1:maximumIterations) * 1e-3, (1:maximumIterations) * 1e-2, (1:maximumIterations) * 1e-1]; 
figure();

for e = 0:3
    routers = ones(100, 1); % The initial state vector
    U = zeros(1, maximumIterations + 1);

    % Calculation of the total costs of the network
    for i = 1:100
        diff = abs(routers(i) - routers(neighbors(G, i)));
        v = zeros(size(diff));
        v(diff == 0) = 2;
        v(diff == 1) = 1;
        iSum(i) = sum(v);
    end
    
    % Estimation of the potential
    U(1) = .5 * sum(iSum);
    
    % Iteration process
    for time = 1:maximumIterations
        I = randi(100); % Randomly select a node
        neigh = neighbors(G, I); % Extract its neighbors
        lSum = zeros(1, 8); % Pre-calculate sums
    
        % The cost function
        for j = 1:8
            diff = abs(routers(neigh) - j);
            v = zeros(size(diff));
            v(diff == 0) = 2;
            v(diff == 1) = 1;
            lSum(j) = sum(v);
        end
    
        P = exp(-eta(time + e * maximumIterations) * lSum) / sum(exp(-eta(time + e * maximumIterations) * lSum)); % Probability calculation
        routers(I) = find(cumsum(P) > rand, 1); % Randomly select a new state for node I
    
        iSum = zeros(100, 1);
    
        % Calculation of the total costs of the network
        for i = 1:100
            diff = abs(routers(i) - routers(neighbors(G, i)));
            v = zeros(size(diff));
            v(diff == 0) = 2;
            v(diff == 1) = 1;
            iSum(i) = sum(v);
        end
        
        % Estimation of the potential
        U(time + 1) = .5 * sum(iSum);
    end
    
    % Plot the potential function for the current eta
    plot(linspace(0, maximumIterations, length(U)), U);
    hold on;
    
    disp(['The potential is ', num2str(U(maximumIterations + 1))]);
    disp(' ');
end

hold off;