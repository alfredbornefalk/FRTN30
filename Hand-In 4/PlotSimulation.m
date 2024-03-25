function PlotSimulation(S, I, R, IPW)
    t = 0:length(S) - 1;
    
    figure();

    subplot(311);
    plot(t, S, 'blue');
    title('Susceptible Individuals');
    xlabel('Week');
    ylabel('Individuals');

    subplot(312);
    plot(t, I, 'red');
    title('Infected Individuals');
    xlabel('Week');
    ylabel('Individuals');

    subplot(313);
    plot(t, R, 'green');
    title('Recovered Individuals');
    xlabel('Week');
    ylabel('Individuals');

    figure();
    plot(t, IPW, 'blue');
    title('Newly Infected Individuals');
    xlabel('Week');
    ylabel('Individuals');
end