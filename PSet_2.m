%% PSet #2
% Omar Thenmalai

%% Question 1 Part D
ap = isallpass([16 84 120 25 0], [25 120 84 16]);
mp = isminphase([8 -2 -1 0], [225 105 -24 -12]);

[phiA, wA] = phasez([16 84 120 25 0], [25 120 84 16], 100);
[phiH, wH] = phasez([2 7 -4], [36 168 165 -75], 100);
[phiM, wM] = phasez([8 -2 -1 0], [225 105 -24 -12], 100);

phiA = phiA*180/(pi);
phiH = phiH*180/(pi);
phiM = phiM*180/(pi);

figure;
hold on
plot(wA, phiA);
plot(wH, phiH);
plot(wM, phiM);
legend('All Pass', 'H(z)', 'Minimum Phase');
legend('Location', 'southwest');
xlabel('Frequency (w)');
title('Filter Phase Responses');
ylabel('Phase (degrees)');
hold off

    