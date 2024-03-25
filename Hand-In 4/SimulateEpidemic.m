function [S, I, R, IPW] = SimulateEpidemic(G, infectedAtStart, beta, rho, weeks, vaccinated)
    if nargin < 6
        vaccinated = zeros(1, weeks + 1);
    else
        if length(vaccinated) ~= weeks
            error('The vaccination vector must have the same length as the corresponding simulation time');
        else
            vaccinated(end + 1) = vaccinated(end);
        end
    end

    n = numnodes(G);
    X = zeros(n, 1);
    V = zeros(n, 1);

    X(randperm(n, infectedAtStart), 1) = 1;
    V(randsample(n, round(vaccinated(1) * n / 1e2))) = 1;

    susceptible = zeros(weeks, 1);
    infected = zeros(weeks, 1);
    recovered = zeros(weeks, 1);
    infected_per_week = zeros(weeks, 1);

    susceptible(1) = n - infectedAtStart;
    infected(1) = infectedAtStart;
    recovered(1) = 0;
    infected_per_week(1) = infectedAtStart;

    for t = 2:weeks + 1
        frac_vacc = (vaccinated(t) - vaccinated(t - 1)) / 1e2;
        vacc_ind = round(frac_vacc * n);
        V(randsample(find(V == 0), vacc_ind)) = 1;

        Xt_next = zeros(n, 1);
        infected_t = 0;

        for I = 1:n
            if X(I) == 0 && V(I) == 0
                neighbor = neighbors(G, I);
                m = sum(X(neighbor) == 1 & V(neighbor) == 0);
                Xt_next(I) = rand < (1 - (1 - beta)^m);
                infected_t = infected_t + Xt_next(I);
            elseif X(I) == 1
                Xt_next(I) = (rand < rho) + 1;
            elseif X(I) == 2
                Xt_next(I) = 2;
            else
                % Do nothing
            end
        end

        X = Xt_next;

        susceptible(t) = sum(X == 0);
        infected(t) = sum(X == 1);
        recovered(t) = sum(X == 2);
        infected_per_week(t) = infected_t;
    end

    S = susceptible;
    I = infected;
    R = recovered;
    IPW = infected_per_week;
end