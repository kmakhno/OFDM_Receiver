clear all;
A = load('rfeDump8k_AWGN20db_TU6_50kmh.mat');

rd = A.rfeDump;
rd_conj = conj(rd);

Nfft = 8192;%number of fft samples for one symbol ofdm
GI = Nfft/4; %guard interval
window_size = GI;
s = [];
c = ones(1,window_size);

for i = 1:(length(rd) - Nfft)
    s(i) = rd_conj(i)*rd(i + Nfft);
end

y = conv(s,c);

% figure(1)
% plot(abs(y));
% grid;

tr = zeros(1, GI+Nfft);
for t = 0:length(y)/(GI+Nfft)-1
    tr = tr + y(t*(window_size+Nfft)+1:(t+1)*(window_size+Nfft));
end

% figure(2)
% plot(abs(tr));
% grid;

[val, ind] = max(tr);
phi = atan2(imag(tr(ind)), real(tr(ind)));

S = [];

for k = 1:length(rd)
    S(k) = rd(k)*exp(-1i*phi*k/(Nfft));
end
%% Check

S_conj = conj(S);
i = 1;
for i = 1:(length(S) - Nfft)
    s_check(i) = S_conj(i)*S(i + Nfft);
end
y_check = conv(s_check,c);
tr_check = zeros(1, GI+Nfft);
for t = 0:length(y_check)/(GI+Nfft)-1
    tr_check = tr_check + y_check(t*(window_size+Nfft)+1:(t+1)*(window_size+Nfft));
end

% figure(3)
% plot(abs(tr_check));
% grid;

[val_check, ind_check] = max(tr_check);
phi_check = atan2(imag(tr_check(ind_check)), real(tr_check(ind_check)));
%% Post-FFT synchronization
symbols_vec = S(ind_check:end-5709);
symbols = [];
for i = 0:25
        symbols(i+1,:) = symbols_vec(i*(Nfft+GI)+1:(i+1)*Nfft+i*GI);
end

% 
% first_symbol = symbols_vec(1:Nfft);
% 
fft_symbol = fft(symbols.', Nfft);
%fft_symbol = fftshift(fft_symbol);
% figure(4);
% plot(abs(fft_symbol(:,1)));
% grid;

H = comm.PNSequence('Polynomial', [11 2 0], ...
                        'SamplesPerFrame', 5616, ...
                        'InitialConditions', [1 1 1 1 1 1 1 1 1 1 1]);

ff = step(H);                  
                    
h = -step(H)*2+1;
h = h';
mask = [1,0,0,0,0,0,0,0,0,0,0,0];
mask = repmat(mask,1,468);
mask1 = mask.*h;

m2 = [0,0,0,1,0,0,0,0,0,0,0,0];
m2 = repmat(m2,1,468);
mask2 = m2.*h;

m3 = [0,0,0,0,0,0,1,0,0,0,0,0];
m3 = repmat(m3,1,468);
mask3 = m3.*h;

m4 = [0,0,0,0,0,0,0,0,0,1,0,0];
m4 = repmat(m4,1,468);
mask4 = m4.*h;

mask1 = [zeros(1,1288) mask1 zeros(1,1288)];
mask2 = [zeros(1,1288) mask2 zeros(1,1288)];
mask3 = [zeros(1,1288) mask3 zeros(1,1288)];
mask4 = [zeros(1,1288) mask4 zeros(1,1288)];
r = zeros(8192,1);
M1 = (fft(mask1'));
M2 = (fft(mask2'));
M3 = (fft(mask3'));
M4 = (fft(mask4'));
for i = 4:4:26
    r = r + abs((ifft((fft(fft_symbol(:,i))).*conj(M3))));
end

plot(r);
