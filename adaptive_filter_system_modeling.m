%% signal parameters - sP: signal, viP: input measurement noise, voP: output measurement noise, N: length of signal vector
sP = 1;
viP = 0.01;
voP = 0.01;
N = 1e5;
%% filter parameters - size: filter size, h: input coloring filter, w0: system to be learned
size = 3;
h = [1 1 1];
w0 = [1 2 -1];
%% Generation of signals - s: uncolored signal, vi: input measurement noise, vo: output measurement noise, u: system input, x: adaptive filter input, d: reference signal
s = sqrt(sP)*randn(1,N);
vi = sqrt(viP)*randn(1,N);
vo = sqrt(voP)*randn(1,N);
u = filter(h,1,s);
x = vi+u;
d = filter(w0,1,u)+vo;
%% estimation of Rx: autocorrelation matrix and p: cross correlation of reference signal and input signal vector
Rx = getRx(x, size);
p = getp(x,d,size);
%% get solution
w = inv(Rx)*p'
function [Rx] = getRx(x, size)
    for ii = 1:size
        xc(ii) = circshift(x,ii-1)*x'/length(x);
    end
    Rx = toeplitz(xc);
end
function p = getp(x,d,size)
    for ii = 1:size
        p(ii) = circshift(x,ii-1)*d'/length(x);
    end
end













