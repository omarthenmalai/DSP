close all

[z,p,k] = ellip(3,2,30,[0.2,0.4]);
[b,a] = zp2tf(z,p,k);
[d0,d1] = tf2ca(b,a);
p0 = fliplr(d0);
p1 = fliplr(d1);


% Part A
[b1,a1] = comp_num_dem(p0,d0,p1,d1);

% Part B
[h,w] = freqz(b,a,1000000);
[h1,w1] = freqz(b1,a1,1000000);
abs_err = max(abs(h-h1)); % Maximum absolute error = 1.4673e-12

% Part C
mag_h = 20*log10(abs(h));
figure;
subplot(2,1,1);
plot(w./pi, mag_h);
title('Magnitude Response of H');
ylabel('Magnitude (dB)');
xlabel('Frequency (normalized radian)');
ylim([-50, 2]);
subplot(2,1,2);
plot(w./pi, unwrap(angle(h))*180/pi)
title('Phase Response of H')
ylabel('Phase (degrees)');
xlabel('Frequency (normalized radian)');

% Part D
bscale = 8;
b_16 = round((b*bscale)*16)/(bscale*16);
b_4 = round((b*bscale)*4)/(bscale*4);
a_16 = round((a*16))/(bscale*16);
a_4 = round((a*4))/(bscale*4);

[h_16, w_16] = freqz(b_16,a_16,1000000);
[h_4,w_4] = freqz(b_4,a_4,1000000);
mag_h_4 = 20*log10(abs(h_4));
mag_h_16 = 20*log10(abs(h_16));

figure;
hold on
plot(w./pi, mag_h);
plot(w_16./pi, mag_h_16);
plot(w_4./pi, mag_h_4);
title('Magnitude Responses of H and Quantized Forms');
ylabel('Magnitude (dB)');
xlabel('Frequency (normalized radian)');
legend('H', 'H rounded to the nearest 16th', 'H rounded to the nearest 4th');
ylim([-50, 10]);
hold off

% Part E
[z_16, p_16, k_16] = tf2zp(b_16, a_16);
[z_4, p_4, k_4] = tf2zp(b_4, a_4);
figure;
subplot(3,1,1);
zplane(z, p);
title('Pole-Zero Plot of Original Filter');
subplot(3,1,2);
zplane(z_16, p_16);
title('Pole-Zero Plot of Original Filter Rounded to Nearest 16th');
subplot(3,1,3);
zplane(z_4, p_4);
title('Pole-Zero Plot of Original Filter Rounded to Nearest 4th');
% The poles for the the original H and when b and a are rounded off to the
% nearest 16th are nearly identical.
% When a and b are rounded to the nearest 4th, the poles have all moved all
% over the place depending on if they were rounded up or rounded down.


% Part F
p0_16 = round((p0*16))/(16);
d0_16 = round((d0*16))/(16);
p1_16 = round((p1)*16)/(16);
d1_16 = round((d1)*16)/(16);

p0_4 = round((p0)*4)/(4);
d0_4 = round((d0)*4)/(4);
p1_4 = round((p1)*4)/(4);
d1_4 = round((d1)*4)/(4);

[b3,a3] = comp_num_dem(p0_16, d0_16, p1_16, d1_16);
[b4, a4] = comp_num_dem(p0_4, d0_4, p1_4, d1_4);

[h1_16, w1_16] = freqz(b3,a3,1000000);
[h1_4, w1_4] = freqz(b4,a4,1000000);
figure;
hold on
plot(w1./pi, 20*log10(abs(h1)));
plot(w1_16./pi, 20*log10(abs(h1_16)));
plot(w1_4./pi, 20*log10(abs(h1_4)));
hold off
ylim([-50, 2]);
title('Magnitude Responses of Original and Rounded off Parallel Allpass Realizations');
ylabel('Magnitude (dB)');
xlabel('Frequency (normalized radian)');
legend('Original', 'Rounded to nearest 16th', 'Rounded to nearest 4th');
[z1, p1, k1] = tf2zp(b1, a1);
[z3, p3, k3] = tf2zp(b3, a3);
[z4, p4, k4] = tf2zp(b4, a4);
figure;
subplot(3,1,1);
zplane(z1, p1);
title('Pole-Zero Plot of Original Parallel Allpass Realization');
subplot(3,1,2);
zplane(z3, p3);
title('Pole-Zero Plot of Parallel Allpass Realization Rounded to Nearest 16th')
subplot(3,1,3);
zplane(z4, p4);
title('Pole-Zero Plot of Parallel Allpass Realization Rounded to Nearest 4th')
% The poles have stayed closer to their positions in the original parallel
% allpass realization. Unlike with the original filter, when rounded to the
% nearest 4th, the poles of the parallel allpass realization stay close to
% the unit circle and do not shift closer to the imaginary-axis


