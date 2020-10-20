%% Question 1
% Part A
close all

A_pass = 2;
A_stop = 30;

f_pass = [1.2e6, 1.5e6];
f_stop = [1.0e6, 1.6e6];

w_pass = 2*pi*f_pass;
w_stop = 2*pi*f_stop;

B = w_pass(2) - w_pass(1);
w0 = sqrt(w_pass(2)*w_pass(1));

w_proto_10 = abs(((w_stop(1)).^2 - w0^2)./(B*w_stop(1)));
w_proto_16 = abs(((w_stop(2)).^2 - w0^2)./(B*w_stop(2)));
% Choose w_proto_16, as magnitude of w_proto_16 < w_proto_10

% Part B
% explitit formula for Butterworth filter order
nButter = ceil(0.5*log10((10^(A_stop/10) - 1)/(10^(A_pass/10)-1))/log10(w_proto_16/1));
% explicit formula for Chebychev filter order
nCheby = ceil(acosh(sqrt((10^(A_stop/10) - 1)/(10^(A_pass/10)-1)))/acosh(w_proto_16/1));
% Bandpass filter orders are 2*n. Bandpass Butter Order: 18. Bandpass Cheby
% Order: 10

% Part C
[n_b, Wn_b] = buttord(w_pass, w_stop, A_pass, A_stop, 's'); % n = 9 -> 2n = 18
[z_b, p_b, k_b] = butter(n_b, Wn_b, 's');

[n_c1, Wn_c1] = cheb1ord(w_pass, w_stop, A_pass, A_stop, 's');
[z_c1, p_c1, k_c1] = cheby1(n_c1, A_pass, Wn_c1, 's'); % n = 5 -> 2n = 10

[n_c2, Wn_c2] = cheb2ord(w_pass, w_stop, A_pass, A_stop, 's');
[z_c2, p_c2, k_c2] = cheby2(n_c2, A_stop, Wn_c2, 's'); % n = 5 -> 2n = 10

[n_e, Wp_e] = ellipord(w_pass, w_stop, A_pass, A_stop, 's');
[z_e, p_e, k_e] = ellip(n_e, A_pass, A_stop, Wp_e, 's');

% Part D
figure;
subplot(2,2,1);
zplane(z_b, p_b);
grid on
title("Pole-Zero Plot for Butterworth Bandpass Filter");
subplot(2,2,2);
zplane(z_c1, p_c1);
grid on
title("Pole-Zero Plot for Chebychev I Bandpass Filter");
subplot(2,2,3);
zplane(z_c2, p_c2);
grid on
title("Pole-Zero Plot for Chebychev II Bandpass Filter");
subplot(2,2,4);
zplane(z_e, p_e);
grid on
title("Pole-Zero Plot for Elliptic Bandpass Filter");

% Part E ??? WHY IS IT 1 for WprotoP
[n_e_lp, Wp_e_lp] = ellipord(1, w_proto_16, A_pass, A_stop, 's');
[z_e_lp, p_e_lp, k_e_lp] = ellip(n_e_lp, A_pass, A_stop, Wp_e_lp, 's');
figure;
zplane(z_e_lp, p_e_lp);
grid on;
title("Pole-Zero Plot for Elliptic Lowpass Prototype Filter");

% Part F
[b_b, a_b] = zp2tf(z_b, p_b, k_b);
[b_c1, a_c1] = zp2tf(z_c1, p_c1, k_c1);
[b_c2, a_c2] = zp2tf(z_c2, p_c2, k_c2);
[b_e, a_e] = zp2tf(z_e, p_e, k_e);

[h_b, wout_b] = freqs(b_b, a_b, 100000);
[h_c1, wout_c1] = freqs(b_c1, a_c1, 100000);
[h_c2, wout_c2] = freqs(b_c2, a_c2, 100000);
[h_e, wout_e] = freqs(b_e, a_e, 100000);

stop_lb = linspace(0, 1000e3);
stop_ub = linspace(1600e3, 3000e3);
passbounds = linspace(1200e3, 1500e3);


butterdb = 20*log10(abs(h_b));
butterphase = unwrap(angle(h_b)*180/pi);
cheby1db = 20*log10(abs(h_c1));
cheby1phase = unwrap(angle(h_c1)*180/pi);
cheby2db = 20*log10(abs(h_c2));
cheby2phase = unwrap(angle(h_c2)*180/pi);
ellipdb = 20*log10(abs(h_e));
ellipphase = unwrap(angle(h_e)*180/pi);
xmax = 3000e3;
ymin = -40;
ymax = 1;
mag_phase_plot("Analog Butterworth", stop_lb, passbounds, stop_ub, butterdb, butterphase, wout_b/(2*pi), xmax, ymin, ymax);
mag_phase_plot("Analog Chebychev I", stop_lb, passbounds, stop_ub, cheby1db, cheby1phase, wout_c1/(2*pi),xmax, ymin, ymax);
mag_phase_plot("Analog Chebychev II", stop_lb, passbounds, stop_ub, 20*log10(abs(h_c2)), cheby2phase, wout_c2/(2*pi),xmax, ymin, ymax);
mag_phase_plot("Analog Elliptic", stop_lb, passbounds, stop_ub, ellipdb, ellipphase, wout_e/(2*pi),xmax, ymin, ymax);

% Part G
figure
hold on
plot(wout_b/(2*pi), butterdb);
plot(wout_c1/(2*pi), 20*log10(abs(h_c1)));
plot(wout_c2/(2*pi), 20*log10(abs(h_c2)));
plot(wout_e/(2*pi), 20*log10(abs(h_e)));
plot(linspace(1.05e6, 1.55e6), -2*ones(100), 'k--');
plot(linspace(1.05e6, 1.55e6), zeros(100), 'k--');
hold off
title('Zoomed In Magnitude Response of the 4 Filters in the Passband');
legend({"Butterworth", "Chebychev I", "Chebychev II", "Elliptic"}, 'Location', 'South'); 
ylabel('Magnitude (dB)');
xlabel('Frequency (MHz)');
ylim([-3, 1]);
xlim([1.05e6, 1.55e6]);

% Part H
butter_lower = interp1(wout_b/(2e3*pi), butterdb, 1000);
butter_upper = interp1(wout_b/(2e3*pi), butterdb, 1600);

cheby1_lower = interp1(wout_c1/(2e3*pi), cheby1db, 1000);
cheby1_upper = interp1(wout_c1/(2e3*pi), cheby1db, 1600);

cheby2_lower = interp1(wout_c2/(2e3*pi), cheby2db, 1000);
cheby2_upper = interp1(wout_c2/(2e3*pi), cheby2db, 1600);

ellip_lower = interp1(wout_e/(2e3*pi), ellipdb, 1000);
ellip_upper = interp1(wout_e/(2e3*pi), ellipdb, 1600);

attenTable = [butter_lower, butter_upper; cheby1_lower, cheby1_upper; cheby2_lower, cheby2_upper; ellip_lower, ellip_upper]


%% Question 2
close all;
clc;
% Part B
sampling_rate = 6000e3;
% f_pass = [1.2e6, 1.5e6];
% f_stop = [1.0e6, 1.6e6];
% A_pass = 2;
% A_stop = 30;

f_pass_bilinear = f_pass/(sampling_rate/2);
f_stop_bilinear = f_stop/(sampling_rate/2);

[n1, Wn1] = buttord(f_pass_bilinear, f_stop_bilinear, A_pass, A_stop); % n = 9 -> 2n = 18
[z1 p1, k1] = butter(n1, Wn1);

[n2, Wn2] = cheb1ord(f_pass_bilinear, f_stop_bilinear, A_pass, A_stop);
[z2, p2, k2] = cheby1(n2, A_pass, Wn2); % n = 5 -> 2n = 10

[n3, Wn3] = cheb2ord(f_pass_bilinear, f_stop_bilinear, A_pass, A_stop);
[z3, p3, k3] = cheby2(n3, A_stop, Wn3); % n = 5 -> 2n = 10

[n4, Wn4] = ellipord(f_pass_bilinear, f_stop_bilinear, A_pass, A_stop);
[z4, p4, k4] = ellip(n4, A_pass, A_stop, Wn4);

% Orders for Digital IIR filters
butter_order_digital = n1*2;
cheby1_order_digital = n2*2;
cheby2_order_digital = n3*2;
ellip_order_digital = n4*2;


%Part C
[b1,a1] = zp2tf(z1,p1,k1);
[b2,a2] = zp2tf(z2,p2,k2);
[b3,a3] = zp2tf(z3,p3,k3);
[b4,a4] = zp2tf(z4,p4,k4);

[h1,w1] = freqz(b1,a1,100000, sampling_rate);
[h2,w2] = freqz(b2,a2,100000, sampling_rate);
[h3,w3] = freqz(b3,a3,100000, sampling_rate);
[h4,w4] = freqz(b4,a4,100000, sampling_rate);

stop_lb = linspace(0, 1000e3);
stop_ub = linspace(1600e3, 3000e3);
passbounds = linspace(1200e3, 1500e3);

mag_phase_plot("Digital Butterworth Filter via Bilinear Transform", stop_lb, passbounds, stop_ub, 20*log10(abs(h1)), unwrap(angle(h1)*180/pi), w1, sampling_rate/2, -40, 1);
mag_phase_plot("Digital Chebychev I Filter via Bilinear Transform", stop_lb, passbounds, stop_ub, 20*log10(abs(h2)), unwrap(angle(h2)*180/pi), w2, sampling_rate/2, -40, 1);
mag_phase_plot("Digital Chebychev II Filter via Bilinear Transform", stop_lb, passbounds, stop_ub, 20*log10(abs(h3)), unwrap(angle(h3)*180/pi), w3, sampling_rate/2, -40, 1);
mag_phase_plot("Digital Elliptic Filter via Bilinear Transform", stop_lb, passbounds, stop_ub, 20*log10(abs(h4)), unwrap(angle(h4)*180/pi), w4, sampling_rate/2, -40, 1);

% Part D
figure('Name', 'Pole-Zero Plots for Digital Filters via Bilinear Transfom');
subplot(2,2,1);
zplane(z1,p1);
title('Butterworth');
subplot(2,2,2);
zplane(z2,p2);
title('Chebychev I');
subplot(2,2,3);
zplane(z3,p3);
title('Chebychev II');
subplot(2,2,4);
zplane(z4,p4);
title('Elliptic')


%% Question 3
clc
close all
% Part A
[b5, a5] = impinvar(b_b, a_b, sampling_rate);
[b6, a6] = impinvar(b_c1, a_c1, sampling_rate);
[b7, a7] = impinvar(b_c2, a_c2, sampling_rate);
[b8, a8] = impinvar(b_e, a_e, sampling_rate);

o1 = filtord(b5, a5);
o2 = filtord(b6, a6);
o3 = filtord(b7, a7);
o4 = filtord(b8, a8);

[h5,w5] = freqz(b5,a5,100000, sampling_rate);
[h6,w6] = freqz(b6,a6,100000, sampling_rate);
[h7,w7] = freqz(b7,a7,100000, sampling_rate);
[h8,w8] = freqz(b8,a8,100000, sampling_rate);

stop_lb = linspace(0, 1000e3);
stop_ub = linspace(1600e3, 3000e3);
passbounds = linspace(1200e3, 1500e3);
mag_phase_plot("Digital Butterworth Filter via Impulse Invariance", stop_lb, passbounds, stop_ub, 20*log10(abs(h5)), unwrap(angle(h5)*180/pi), w1, sampling_rate/2, -40, 1);
mag_phase_plot("Digital Chebychev I Filter via Impulse Invariance", stop_lb, passbounds, stop_ub, 20*log10(abs(h6)), unwrap(angle(h6)*180/pi), w2, sampling_rate/2, -40, 1);
mag_phase_plot("Digital Chebychev II Filter via Impulse Invariance", stop_lb, passbounds, stop_ub, 20*log10(abs(h7)), unwrap(angle(h7)*180/pi), w3, sampling_rate/2, -40, 1);
mag_phase_plot("Digital Elliptic Filter via Impulse Invariance", stop_lb, passbounds, stop_ub, 20*log10(abs(h8)), unwrap(angle(h8)*180/pi), w4, sampling_rate/2, -40, 1);

[z5, p5, k5] = tf2zp(b5,a5);
[z6, p6, k6] = tf2zp(b6,a6);
[z7, p7, k7] = tf2zp(b7,a7);
[z8, p8, k8] = tf2zp(b8,a8);

figure('Name', 'Pole-Zero Plots for Digital Filters via Impulse Invariance');
subplot(2,2,1);
zplane(z5, p5);
title('Butterworth');
subplot(2,2,2);
zplane(z6, p6);
title('Chebychev 1');
subplot(2,2,3);
zplane(z7, p7);
title('Chebychev 2');
subplot(2,2,4);
zplane(z8, z8);
title('Elliptic');

% Part D
figure('Name', 'Butterworth');
subplot(2,1,1);
impulse(zpk(z_b, p_b, k_b));
title('Impulse Response of Analog Filter')
subplot(2,1,2);
impz(b5, a5, 211, sampling_rate);
title('Impulse Response of Digital Filter');

figure('Name', 'Chebychev I');
subplot(2,1,1);
impulse(zpk(z_c1, p_c1, k_c1));
title('Impulse Response of Analog Filter')
subplot(2,1,2);
impz(b6, a6, 481, sampling_rate);
title('Impulse Response of Digital Filter');

figure('Name', 'Chebychev II');
subplot(2,1,1);
impulse(zpk(z_c2, p_c2, k_c2));
title('Impulse Response of Analog Filter')
subplot(2,1,2);
impz(b7, a7, 181, sampling_rate);
title('Impulse Response of Digital Filter');

figure('Name', 'Elliptic');
subplot(2,1,1);
impulse(zpk(z_e, p_e, k_e));
title('Impulse Response of Analog Filter')
subplot(2,1,2);
impz(b8, a8, 241, sampling_rate);
title('Impulse Response of Digital Filter');

%% Question 4
close all
%a
cheb_wind = chebwin(31, 30);
[h_c, w_c] = freqz(cheb_wind, 1000);
h_c = h_c / h_c(1);

%b
[z_c, p_c, k_c] = tf2zpk(cheb_wind);
figure;
zplane(z_c, p_c);
title('Zeros of the Chebychev Window');
c_mlobe_w = abs((2*z_c(1))/(4*pi/31));
%mainlobe width= 4.9338

%c
figure;
plot(w_c/pi, 20*log10(abs(h_c)));
title('Magnitude Response of the Chebychev and Kaiser Windows');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');

%e & f
beta = 3.02;
hold on;
kaiser_wind = kaiser(31, beta);
[h_k, f_k] = freqz(kaiser_wind, 1000); 
h_k = h_k / h_k(1);
[z_k, p_k, k_k] = tf2zpk(kaiser_wind);
plot(f_k/pi, 20*log10(abs(h_k)));
legend('Chebychev Window', 'Kaiser Window');

figure;
zplane(z_c, z_k);
title('Zeros of the Chebychev and Kaiser Windows');
legend('Chebychev', 'Kaiser');

%g    
%beta = 3.02

%h
wvtool(kaiser_wind);
%peak sidelobe level = 24.8db

%i
cheb_etotal = sum((h_c.^2));
w0_cheb = islocalmin(abs(h_c));
w0_cheb = find(w0_cheb, 1, 'first');
cheb_elobe = sum(abs(h_c(w0_cheb:end).^2));
frac_cheb_elob = cheb_elobe / cheb_etotal;

kaiser_etotal = sum((h_k.^2));
w0_kaiser = islocalmin(abs(h_k));
w0_kaiser = find(w0_kaiser, 1, 'first');
kaiser_elobe = sum(abs(h_k(w0_kaiser:end).^2));
frac_kaiser_elob = kaiser_elobe / kaiser_etotal;


%% Question 5
close all
%a
deltaPass = abs((10^(2/20)-1)/(10^(2/20)+1));
deltaStop = abs(10^(-30/20));

fsamp = 6e6;
fcuts = [1e6 1.2e6 1.5e6 1.6e6];
mags = [0 1 0];
devs = [deltaStop deltaPass deltaStop];
[n_k, Wn_k, beta_k, ftype_k] = kaiserord(fcuts, mags, devs, fsamp);
%Kaiser window order = 93
f_k = kaiser(n_k+2, beta_k);
b_k = fir1(n_k+1, Wn_k, ftype_k, f_k, 'noscale');

[n_f, fo_f, ao_f, w_f] = firpmord(fcuts, mags, devs, fsamp);
%Park-McClellan window order = 56
b_f = firpm(n_f+2, fo_f, ao_f, w_f);

%b
figure;
subplot(2,1,1)
stem(1:length(b_k),b_k);
title('Stem Plot for Kaiser Filter');
xlabel('n');
ylabel('h[n]');

subplot(2,1,2)
stem(1:length(b_f), b_f);
title('Stem Plot for Parks-McClellan Filter');
xlabel('n');
ylabel('h[n]');

stop_lb = linspace(0, 1);
stop_ub = linspace(1.6, 3);
passbounds = linspace(1.2, 1.5);

[h_k, f_k] = freqz(b_k, 1, 1024, fsamp);
[h_f, fo_f] = freqz(b_f, 1, 1024, fsamp);

figure;
subplot(2,1,1);
hold on
plot(f_k/(1e6), 20*log10(abs(h_k)));
plot(stop_lb, -30*ones(length(stop_lb)), 'k--');
plot(stop_ub, -30*ones(length(stop_ub)), 'k--');
plot(passbounds, -2*ones(size(passbounds)), 'k--')
plot(passbounds, zeros(size(passbounds)), 'k--');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
title('Magnitude Response of Kaiser Filter');
ylim([-50 2]);

subplot(2,1,2);
hold on
plot(fo_f/(1e6), 20*log10(abs(h_f)));
plot(stop_lb, -30*ones(length(stop_lb)), 'k--');
plot(stop_ub, -30*ones(length(stop_ub)), 'k--');
plot(passbounds, -2*ones(size(passbounds)), 'k--')
plot(passbounds, zeros(size(passbounds)), 'k--');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
title('Magnitude Response of Parks-McClellan Filter');
ylim([-50 2]);

%c
mag_pass_lb_k = h_k(find(((f_k/(1e6) >= 1.2) & (f_k/(1e6) <= 1.5) )));
mag_stop_lb_k = h_k(find(((f_k/(1e6) >= 0) & (f_k/(1e6) <= 1) )));
mag_stop_ub_k = h_k(find(((f_k/(1e6) >= 1.6) & (f_k/(1e6) <= 3) )));
precise_var_k = 20*log10(abs(max(mag_pass_lb_k)-min(mag_pass_lb_k)))
min_atten_lb_k = 20*log10(abs(max(mag_stop_lb_k)))
min_atten_ub_k = 20*log10(abs(max(mag_stop_ub_k)))
%Kaiser
%precise variation = 2.3773dB
%minimum attenuation of lower stopband = -34.1127dB
%minimum attenuation of upper stopband = -28.9072dB

mag_pass_lb_f = h_f(find(((f_f/(1e6) >= 1.2) & (f_f/(1e6) <= 1.5) )));
mag_stop_lb_f = h_f(find(((f_f/(1e6) >= 0) & (f_f/(1e6) <= 1) )));
mag_stop_ub_f = h_f(find(((f_f/(1e6) >= 1.6) & (f_f/(1e6) <= 3) )));
precise_var_f = 20*log10(abs(max(mag_pass_lb_f)-min(mag_pass_lb_f)))
min_atten_lb_f = 20*log10(abs(max(mag_stop_lb_f)))
min_atten_ub_f = 20*log10(abs(max(mag_stop_ub_f)))
%Park-McClellan
%precise variation = 2.3773dB
%minimum attenuation of lower stopband = -34.1127dB
%minimum attenuation of upper stopband = -28.9072dB

%d
weight = deltaPass/deltaStop;
checkweight = (weight == w_f(find(w_f, 1, 'first')));
%checkweight = 1

%e
%Increasing the order by two decreases the strange ripple in the lower
%transition band.

