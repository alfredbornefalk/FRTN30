function G = RandomGraph(k, n)
    nodesAtStart = floor(k) + 1;
    W = ~diag(ones(1, nodesAtStart));
    c = k / 2;

    for t = 1:n - nodesAtStart
        k0 = length(W);
        w = sum(W, 2);
        P = w ./ sum(w);

        links = floor(c) + uint32(rand < (c - floor(c)));

        for j = 1:links
            neighbor = randsample(k0, 1, true, P);
            P(neighbor) = 0;
            W(k0 + 1, neighbor) = 1;
            W(neighbor, k0 + 1) = 1;
        end
    end

    G = graph(W);
end