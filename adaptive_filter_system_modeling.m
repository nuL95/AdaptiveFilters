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
p = getp(x,d,size)';
%% get solution
w = inv(Rx)*p;
%% find solution using steepest descent
steps = 100;
ws = zeros(size,steps);
Rxe = eig(Rx);
mu = 1/(min(Rxe)+max(Rxe));
for k = 1:steps
    ws(:,k+1)=(eye(size)-2*mu*Rx)*ws(:,k)-2*mu*p;
    pf(k)= ws(:,k)'*Rx*ws(:,k)-2*ws(:,k)'*p;
end
%% plot - pf: performance surface showing iterative approach to solution, w: theoretical ideal filter, ws: ideal filter found using steepest descent
plot(pf)
w(:,end)
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













