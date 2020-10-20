function 1fplot(filtername, h, stop_lb, passbounds, stop_ub, db, phase, wout)
    figure;
    subplot(2,1,1);
    hold on;
    plot(wout/(2e3*pi), db);
    plot(stop_lb, -30*ones(length(passbounds)), 'k--');
    plot(passbounds, -2*ones(length(passbounds)), 'k--');
    plot(passbounds, 0*ones(length(passbounds)), 'k--');
    plot(stop_ub, -30*ones(length(stop_ub)), 'k--');
    hold off
    title(['Magnitude Response of ', filtername, 'Filter']);
    ylabel("Magnitude (dB)");
    xlabel("Frequency (kHz)");
    xlim([0 3000]);
    ylim([-40 5]);
    subplot(2,1,2);
    plot(wout/(2e3*pi), unwrap(angle(h)*180/pi));
    title(['Phase Response of ', filtername, 'Filter']);
    ylabel("Phase (degrees)");
    xlabel("Frequency (kHz)");
end