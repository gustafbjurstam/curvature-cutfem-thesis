function plot_cond(taus, A, h, rates, filename, refs, export)
    % Plots all
    figure(Position="default")
    markers = ["-o", "-+", "-x", "-square", "-^", "-*"];
    taus_str = [];
    for k = 1:length(taus)
        if taus(k) == 0
            taus_str = [taus_str;"0"];
        else
            taus_str = [taus_str;"10^{" + num2str(log10(taus(k))) + "}"];
        end
    end

    for k = 1:length(taus)
        loglog(h, A(:,k), markers(k))
        hold on
    end
    legger = "$\tau = " + taus_str + "$, rate $\approx" + num2str(round(rates(:),2)) + "$";
    legend(legger, Location="northeast", Interpreter="latex", AutoUpdate="off")
    for ref = refs
        loglog([h(end),h(1)], geomean((A(1,:)))*[(h(end)/h(1))^ref, 1], "--black")
    end
    xlabel("$h$")
    ylabel("Approximate Condition Number $\kappa_1(A)$")
    axis tight
    hold off
    if export
        exportgraphics(gcf, filename, 'ContentType', 'vector');
    end
end