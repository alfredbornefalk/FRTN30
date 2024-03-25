%% Good practice
clear; clc; close all;
%% Loading the data
load -ASCII traffic.mat;
load -ASCII capacities.mat;
load -ASCII traveltime.mat;
load -ASCII flow.mat;
%% Part (a): The Shortest Path
% Source nodes
s = [1, 2, 3, 4, 1, 6, 7, 8, 9, 2, 3, 3, 4, 5, 6, 10, 10, 7, 8, 9, 11, 12, 13, 11, 13, 14, 15, 16];

% Target nodes
t = [2, 3, 4, 5, 6, 7, 8, 9, 13, 7, 8, 9, 9, 14, 10, 11, 15, 10, 11, 12, 12, 13, 14, 15, 17, 17, 16, 17];

G1 = digraph(s, t, traveltime); % Travel time graph
figure(); % Figure 1
plot(G1, 'Layout', 'force', 'EdgeLabel', G1.Edges.Weight);

sp = shortestpath(G1, 1, 17); % Shortest path from 1 to 17

disp(['The shortest path from 1 to 17: ', num2str(sp), '.']);
disp(' ');
%% Part (b): The Maximum Flow
G2 = digraph(s, t, capacities); % Capacity graph
figure(); % Figure 2
plot(G2, 'Layout', 'force', 'EdgeLabel', G2.Edges.Weight);

mf = maxflow(G2, 1, 17); % Maximum flow from 1 to 17

disp(['The maximum flow from 1 to 17: ', num2str(mf), '.']);
disp(' ');
%% Part (c): The Exogenous Inflow or Outflow at Each Node
exflows = traffic * flow; % sum(exflows) = 0 -> exogenous net flow

G3 = digraph(s, t, flow); % Flow graph
figure(); % Figure 3
plot(G3, 'Layout', 'force', 'EdgeLabel', G3.Edges.Weight);

disp(['The external inflow (or outflow) from 1 to 17: ', num2str(exflows') '.']);
disp(' ');
%% Preparations for (d)–(g)
M = 28;
lambda = zeros(17, 1); % All net inflows set to zero (except on node 1)
mu = lambda;
lambda(1) = exflows(1); % The net inflow on node 1
mu(end) = lambda(1);
%% Part (d): The Social Optimum
% Algorithm
cvx_begin
    variable f(M)
    minimize sum(traveltime .* capacities .* inv_pos(1 - f ./ capacities) - traveltime .* capacities)
    subject to
        traffic * f == lambda - mu;
        0 <= f <= capacities;
cvx_end

G4 = digraph(s, t, f);
figure(); % Figure 4
plot(G4, 'Layout', 'force', 'EdgeLabel', round(G4.Edges.Weight));
%% Part (e): The Wardrop Equilibrium
% Algorithm
cvx_begin
    variable g(M)
    minimize sum(-traveltime .* capacities .* log(1 - g ./ capacities))
    subject to
        traffic * g == lambda - mu;
        0 <= g <= capacities;
cvx_end

G5 = digraph(s, t, g);
figure(); % Figure 5
plot(G5,'Layout', 'force', 'EdgeLabel', round(G5.Edges.Weight))
%% Part (f): The New Wardrop Equilibrium
omega1 = f .* (traveltime .* capacities ./ ((capacities - f) .^ 2)); % Tolls

cvx_begin
    variable h(M)
    minimize sum(-traveltime .* capacities .* log(1 - h ./ capacities) + omega1 .* h)
    subject to
        traffic * h == lambda - mu;
        0 <= h <= capacities;
cvx_end

G6 = digraph(s,t,h);
figure(); % Figure 6
plot(G6,'Layout','force','EdgeLabel', round(G6.Edges.Weight));
%% Part (g): The New System Optimum
% The new system optimum
cvx_begin
	variable p(M)
	p1 = 0;

	for k = 1:M
		p1 = p1 + traveltime(k) * quad_over_lin(p(k), capacities(k) - p(k));
    end

	minimize p1
	subject to
        traffic * p == lambda - mu;
        0 <= p <= capacities;
cvx_end	

G7 = digraph(s, t, p);
figure(); % Figure 7
plot(G7,'Layout','force','EdgeLabel', round(G7.Edges.Weight))

% Verification w/ new Wardrop equilibrium
omega2 = p .* (traveltime .* capacities ./ ((capacities - p) .^ 2)) - traveltime;

cvx_begin
    variable q(M)
    minimize sum(-traveltime .* capacities .* log(1 - q ./ capacities) + omega2 .* q)
    subject to
        traffic * q == lambda - mu;
        0 <= q <= capacities;
cvx_end

G8 = digraph(s, t, q);
figure(); % Figure 8
plot(G8, 'Layout', 'force', 'EdgeLabel', round(G8.Edges.Weight));