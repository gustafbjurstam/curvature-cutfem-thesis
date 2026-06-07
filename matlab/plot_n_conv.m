function plot_n_conv(n, rate, h, exps)
    figure(Position="default")
    loglog(h, n, "-o")
    hold on
    xlabel("$h$")
    for expp = exps
        loglog([h(end),h(1)], n(1)*[(h(end)/h(1))^expp, 1], "--black")
    end
    axis tight
end