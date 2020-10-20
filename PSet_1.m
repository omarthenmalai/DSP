%% Omar Thenmalai PSet #1

%% Question 4
% Part C
T = 2;
N = 10;
t = -1000:1000;
syms n;
x = symsum(1-abs((t-n*T-T/2)/(T/2)),n,0,N-1);
x = double(x);
figure;
plot(t, x);
ylabel('x(t)');
xlabel('time');
f = linspace(-2*pi,2*pi, 301);
X = abs(N*diric(2*pi*f*T, N).*exp(-1i*pi*f*T*N).*(sinc(f*T/2).^2));
figure;
plot(f, X);
ylabel('|X(f)|');
xlabel('Frequency');
%% Question 6
f = 20000;
amp = 2;
l = 1000;
ts = 1/100000; 
t = (0:l-1)*ts;
x = amp*sin(2*pi*f*t); % Generate 20kHz signal sampled at 100kHz
x_with_noise = x + sqrt(0.2).*randn(1,size(t,2)); % Add noise with variance 0.2
X = fft(x_with_noise, 1024); % take the 1024 point DFT
mag = abs(X);
mag = fftshift(mag);
window = chebwin(1000, 30); % 1000 with 30dB chebyshev window

windowed = mag.*window;
n=length(X);
freq = (-n/2:n/2-1)*100/n;
figure;
plot(freq,windowed); % fft with noise
ylabel('Power');
xlabel('Frequency');
%% Question 7
% Part A
H = tf([-0.72, 0.1, 1], [0.9025, 0.95, 1]);
H = tf([1 0.1 -0.72], [1 0.95 0.9025]);
[z,p,k] = tf2zp([1 0.1 -0.72],[1 0.95 0.9025]);
figure;
zplane(z, p);

v = sqrt(4).*randn(1,100000);
x = filter([1 0.1 -0.72], [1 0.95 0.9025], v);
[s_est, w] = pwelch(x, hamming(512), 256, 512);
s_est = s_est(:)/s_est(1);
S = [];
for j=1:size(w) 
    S = [S; 4*abs((1+0.1*exp(-1i*w(j)) + -0.72*exp(-1i*2*w(j)))/(1+0.95*exp(-1i*w(j))+0.9025*exp(-1i*2*w(j)))).^2];
end
S = S(:)/S(1);
figure; 
hold on;  
plot(w, s_est);
plot(w, S);
hold off;
legend('Periodogram Estimate', 'PSD');
xlabel('Radian Frequency');
ylabel('Power');
