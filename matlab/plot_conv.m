function plot_conv(taus, A, h, rates, filename, refs, export, ylims)
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
    legger = "$\tau = " + taus_str + "$, rate $\approx" +num2str(round(rates,2)) + "$";
    legend(legger, Location="southeast", Interpreter="latex", AutoUpdate="off")
    for ref = refs
        loglog([h(end),h(1)], geomean(max(A,[], 1))*[(h(end)/h(1))^ref, 1], "--black")
    end
    xlabel("$h$")
    ylabel("$| \kern-0.7pt | H_h - H^e | \kern-0.7pt |_{L^2(\Gamma_h)}$")
    axis tight
    hold off
    ylim(ylims)
    if export
        exportgraphics(gcf, filename, 'ContentType', 'vector');
    end
end