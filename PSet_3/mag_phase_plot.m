function mag_phase_plot(filtername, stop_lb, passbounds, stop_ub, db, phase, wout, xmax, ymin, ymax)
    figure('Name', filtername)
    
    % Magnitude Response plot
    subplot(2,1,1);
    hold on;
    plot(wout, db);
    plot(stop_lb, -30*ones(length(passbounds)), 'k--');
    plot(passbounds, -2*ones(length(passbounds)), 'k--');
    plot(passbounds, 0*ones(length(passbounds)), 'k--');
    plot(stop_ub, -30*ones(length(stop_ub)), 'k--');
    hold off
    title(['Magnitude Response']);
    ylabel("Magnitude (dB)");
    xlabel("Frequency (MHz)");
    xlim([0 xmax]);
    ylim([ymin ymax]);
    
    % Phase Response plot
    subplot(2,1,2);
    plot(wout, phase);
    title(['Phase Response'])
    ylabel("Phase (degrees)");
    xlabel("Frequency (MHz)");
    xlim([0 xmax]);
end