% Part A
h = [0.15774243, 0.69950381, 1.06226376, 0.44583132, -0.31998660, ...
    -0.18351806, 0.13788809, 0.03892321, -0.04466375, ...
    -7.83251152E-4, 6.75606236E-3, -1.52353381E-3];
N = length(h);
H0 = h;
H1 = zeros(1,N);
F0 = zeros(1,N);
F1 = zeros(1,N);
for k=0:(N-1)
    H1(k+1) = (-1)^k*h(N-k);
    F0(k+1) = h(N-k);
    F1(k+1) = (-1)^(k+1)*h(k+1);
end
   
% Part B
M = 2;
t = (1/M)*(conv(F0, H0) + conv(F1,H1));
t_error = max([abs(t(1:11)) abs(t(13:length(t)))]);
c = t(N); % c = 2, N = 11
H1_alternate = zeros(1,N);
H0_alternate = zeros(1,N);
for k=1:N
    H1_alternate(k) = H1(k)*(-1)^(k-1);
    H0_alternate(k) = H0(k)*(-1)^(k-1);
end
a = (1/M)*(conv(F0,H0_alternate) + conv(F1, H1_alternate));
a_error = max(abs(a))

% Part C
[H0_w, w] = freqz(H0, 1, 1000);
[H1_w, w] = freqz(H1, 1, 1000);
[F0_w, w] = freqz(F0, 1, 1000);
[F1_w, w] = freqz(F1, 1, 1000);

figure('Name', 'Magnitude Responses of |H0(w)| and |H1(w)|');
plot(w, abs(H0_w));
hold on;
plot(w, abs(H1_w));
xlim([0,pi]);
xlabel('Frequency');
ylabel('Magnitude');
hold off;
legend('|H0(w)|', '|H1(w)|');


% Part D
const = abs(H0_w).^2 + abs(H1_w).^2; % The constant is 4!
err = max(abs(H0_w) - flipud(abs(H1_w))); % Max error is approximately 0.0069

% Part E
E00 = H0(1:2:end);
E01 = H0(2:2:end);
E10 = H1(1:2:end);
E11 = H1(2:2:end);

p00 = conv(flip(E00), E00) + conv(flip(E10), E10);
plength = length(p00);
mid = ceil(plength/2)
p00_error = max([abs(p00(1:mid-1)) abs(p00(mid+1:plength))])% 0.0011
p01 = conv(flip(E00), E01) + conv(flip(E10),E11); 
p01_error = max([abs(p01(1:mid-1)) abs(p01(mid+1:plength))]) % 2.7756e-17
p10 = conv(flip(E01), E00) + conv(flip(E11), E10);
p10_error = max([abs(p10(1:mid-1)) abs(p10(mid+1:plength))])% 2.7756e-17
p11 = conv(flip(E01), E01) + conv(flip(E11),E11); 
p11_error = max([abs(p11(1:mid-1)) abs(p11(mid+1:plength))]) % 0.0011
d = p00(mid); % d = 2
d_error = abs(p00(mid) - p11(mid)); % 2.2204e-16