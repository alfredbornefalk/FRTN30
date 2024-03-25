%% 1: Preliminary Parts
clear; clc; close all; % Good practice
%% 1.1: Epidemic on a Known Graph
% Setting up the parameters
n = 5e2;
beta = .3; 
rho = .7;
weeks = 15;
infectedAtStart = 10;
N = 1e2;

% Generating the sparce adjacency matrix
W = zeros(n);
W = W + diag(ones(n - 1, 1), 1);
W = W + diag(ones(n - 1, 1), -1);
W = W + diag(ones(n - 2, 1), 2);
W = W + diag(ones(n - 2, 1), -2);
W = W + diag(ones(1, 1), n - 1);
W = W + diag(ones(1, 1), 1 - n);
W = W + diag(ones(2, 1), n - 2);
W = W + diag(ones(2, 1), 2 - n);
W = sparse(W);
G = graph(W);

S = zeros(N, weeks + 1);
R = zeros(N, weeks + 1);
I = zeros(N, weeks + 1);
IPW = zeros(N, weeks + 1);

% Simulation
for j = 1:N
    [S(j, :), I(j, :), R(j, :), IPW(j, :)] = SimulateEpidemic(G, infectedAtStart, beta, rho, weeks);
end

% Calculating the means
meanS = mean(S);
meanI = mean(I);
meanR = mean(R);
meanIPW = mean(IPW);

PlotSimulation(meanS, meanI, meanR, meanIPW);
%% 1.2: Generating a Random Graph
clear; % Necessary

k1 = 2;
k2 = 3;
k3 = 4.6;
k4 = 6;
k5 = 10;
n = 1e3; % # nodes in the network

G1 = RandomGraph(k1, n);
G2 = RandomGraph(k2, n);
G3 = RandomGraph(k3, n);
G4 = RandomGraph(k4, n);
G5 = RandomGraph(k5, n);

averageDegree1 = mean(degree(G1));
averageDegree2 = mean(degree(G2));
averageDegree3 = mean(degree(G3));
averageDegree4 = mean(degree(G4));
averageDegree5 = mean(degree(G5));

disp(['k1 is ', num2str(k1), ' and the average degree of G1 is ', num2str(averageDegree1)]);
disp(['k2 is ', num2str(k2), ' and the average degree of G2 is ', num2str(averageDegree2)]);
disp(['k3 is ', num2str(k3), ' and the average degree of G3 is ', num2str(averageDegree3)]);
disp(['k4 is ', num2str(k4), ' and the average degree of G4 is ', num2str(averageDegree4)]);
disp(['k5 is ', num2str(k5), ' and the average degree of G5 is ', num2str(averageDegree5)]);
%% 2: Simulating a Pandemic Without Vaccination
clear; % Necessary

% Setting up the parameters
n = 5e2; % # nodes in the network
k = 6;
beta = .3;
rho = .7;
weeks = 15;
infectedAtStart = 10;
N = 1e2;

G = RandomGraph(k, n);

S = zeros(N, weeks + 1);
R = zeros(N, weeks + 1);
I = zeros(N, weeks + 1);
IPW = zeros(N, weeks + 1);

% Simulation
for j = 1:N
    [S(j, :), I(j, :), R(j, :), IPW(j, :)] = SimulateEpidemic(G, infectedAtStart, beta, rho, weeks);
end

% Calculating the means
meanS = mean(S);
meanI = mean(I);
meanR = mean(R);
meanIPW = mean(IPW);

PlotSimulation(meanS, meanI, meanR, meanIPW);
%% 3: Simulating a Pandemic With Vaccination
clear; % Necessary

n = 5e2; % # nodes in the network
k = 6;
beta = .3;
rho = .7;
weeks = 15;
infectedAtStart = 10;
N = 1e2;

% Total fraction of vaccinated people
vaccinated = [0, 5, 15, 25, 35, 45, 55, 60, 60, 60, 60, 60, 60, 60, 60];

G = RandomGraph(k, n);

S = zeros(N, weeks + 1);
R = zeros(N, weeks + 1);
I = zeros(N, weeks + 1);
IPW = zeros(N, weeks + 1);

% Simulation
for j = 1:N
    [S(j, :), I(j, :), R(j, :), IPW(j, :)] = SimulateEpidemic(G, infectedAtStart, beta, rho, weeks, vaccinated);
end

% Calculating the means
meanS = mean(S);
meanI = mean(I);
meanR = mean(R);
meanIPW = mean(IPW);

PlotSimulation(meanS, meanI, meanR, meanIPW);
%% 4: The H1N1 Pandemic in Sweden 2009
clear; % Necessary

% Setting up the parameters
n = 934;
weeks = 15;
infectedAtStart = 1; 
N = 10;

% Vaccination vector
vaccinated = [5, 9, 16, 24, 32, 40, 47, 54, 59, 60, 60, 60, 60, 60, 60];

% True vector for newly infected each week (scaled down)
I0 = [1, 3, 5, 9, 17, 32, 32, 17, 5, 2, 1, 0, 0, 0, 0];

minRMSE = realmax; % Initial RMSE (as high as possible)

%Initial guesses as suggested by the problem description
kDelta = 1;
betaDelta = .1;
rhoDelta = .1;

% Best parameter estimates in the simulations
k0 = 6.7;
beta0 = .3;
rho0 = .7;

IBest = zeros(16, 1); % Used for figure 8 in the report

figure();

% The algorithm
while true
    for k = k0 - kDelta : kDelta : k0 + kDelta
        
        G = RandomGraph(k, n);
        
        for beta = beta0 - betaDelta : betaDelta : beta0 + betaDelta
            for rho = rho0 - rhoDelta : rhoDelta : rho0 + rhoDelta
                IPWT = zeros(weeks + 1, 1);
                
                for I = 1:N
                    [~, ~, ~, IPW] = SimulateEpidemic(G, infectedAtStart, beta, rho, weeks, vaccinated);
                    IPWT = IPWT + IPW;
                end
               
                I = IPWT / N;
                RMSE = sqrt(mean((I(2:end) - I0').^2)); % Calculating RMSE
                parameters = [k, beta, rho];
                
                % Checking if a new minimum RMSE is found
                if RMSE < minRMSE
                    minRMSE = RMSE; % Update
                    minParameters = parameters; % Save corresponding parameters
                    
                    % Plot new minimum (figure 7)
                    t = 0:15;
                    plot(t, I);
                    hold on;
                    plot(t, [1 I0]);
                    drawnow;

                    IBest = I; % Saving the best infected vector
                end
            end
        end
    end
    
    % If the same parameters are found, the algorithm should be stopped
    if isequal(minParameters, [k0, beta0, rho0])
        break
    end
    
    % Update newly best found parameters
    k0 = minParameters(1);
    beta0 = minParameters(2);
    rho0 = minParameters(3);
end

% Plot the minimum (figure 8)
figure()
t = 0:15;
plot(t, IBest, 'blue');
hold on;
plot(t, [1 I0], 'red');
drawnow;
title('Newly Infected Individuals');
xlabel('Week');
ylabel('Individuals');

disp(['The algorithm is done; the best parameters are k = ', num2str(k0), ', beta = ', num2str(beta0), ', rho = ', num2str(rho0)]);
disp(['RMSE = ', num2str(minRMSE)]);

% Plotting figure 9
N = 100;

G = RandomGraph(k0, n);

S = zeros(N, weeks + 1);
R = zeros(N, weeks + 1);
I = zeros(N, weeks + 1);
IPW = zeros(N, weeks + 1);

% Simulation
for j = 1 : N
    [S(j, :), I(j, :), R(j, :), IPW(j, :)] = SimulateEpidemic(G, infectedAtStart, beta0, rho0, weeks, vaccinated);
end

% Calculating the means
meanS = mean(S);
meanI = mean(I);
meanR = mean(R);
meanIPW = mean(IPW);

% Plotting the average simulation result
PlotSimulation(meanS, meanI, meanR, meanIPW)