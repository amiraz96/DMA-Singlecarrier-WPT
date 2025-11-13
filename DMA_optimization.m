clear all;
clc;

%% System Parameters
f_0 = 5e9; % frequency
c0 = 3e8; %% Light Speed
lambda = c0/f_0; %% WaveLength
antenna_length = 0.2; % antenna length in meters
P_max = 1e-3;
D = sqrt(2)*antenna_length;
frensel = ((D^4)/(8*lambda))^(1/3);
fraunhofer = (2*(D^2))/lambda;
boresight_gain = 2; 
x_min = 0; x_max = 10; y_min = 0; y_max = 10; z_min = 0; z_max = 3; % coordinates of your area

%% Define your user locations and number
M = 1; % number of users
user_loc = zeros(M, 3);
RF_Thr = ones(M, 1).*1e-6; % your power thresholds
%% DMA
DMA_loc = [(x_min + x_max)/2 (y_min + y_max)/2 z_max];
[DMA_element_loc, H, N_d, N_e] = DMA_deploy(f_0, antenna_length, DMA_loc);
N = N_d*N_e;
L = N_d*N;

%% Channel
channel_vec = Do_Channels(DMA_element_loc, user_loc, boresight_gain, lambda);

%% Just creating a Lorentzian vector with a resolution for comparison in final stage
lorentz_vec = linspace(0, 2*pi, 5000);
lorentz_value = zeros(length(lorentz_vec), 1);
for i = 1:length(lorentz_vec)
    lorentz_value(i) = Q_DMA(lorentz_vec(i));
end

%% This is the procedure of forming the Q matrix
% initialize a random Q
phi_vec = rand(N_d,N_e); % imagine these are your phi values for the elements in range [0, 1]
Q = zeros(N, N_d); % Q structure
for i = 1:N_d
    for l = 1:N_e
        % Note that the input of Q_dma function is a value between 0, 2\pi
        Q((i - 1)*N_e + l, i) = Q_DMA(phi_vec(i,l) * 2 * pi); % construct Q using the values between [0, 1]
    end
end

%% Initial Move based on random Q, you will find the details of the equations in the paper and the algorithm
q = reshape(Q, [L, 1]);
idx_vec = find(q~=0);
a = channel_vec;
Bk = zeros(M, N_d, N_d);
for k = 1:M
    bk = (a(:, k))'*H*Q;
    Bk(k, :, :) = transpose(bk'*bk);
end
c = H*Q;
C = transpose(c'*c);
cvx_begin sdp 
    variable W(N_d,N_d) hermitian semidefinite
    minimize(real(trace(transpose(W)*C)))
    subject to
    for kk = 1:M
        real(trace(transpose(W)*reshape(Bk(kk, :, :), [N_d N_d]))) >= RF_Thr(kk);
    end
cvx_end
[W1, W2] = eig(W);
[sorted_eigval,  sorted_eigval_idx] = sort(diag(W2));
w = zeros(N_d, M);
for m = 1:M
    w(:, m) = sqrt(sorted_eigval(N_d - m + 1)).*W1(:, sorted_eigval_idx(N_d - m + 1));
end
Zmkbar = zeros(M, N, N);
Zmbar = zeros(N, N);
for k = 1:M
    for m = 1:M
        zmk = kron(transpose(w(:,m)), (a(:, k))'*H);
        zmkbar = zeros(1,N);
        for i = 1:length(idx_vec)
            zmkbar(i) = zmk(idx_vec(i));
        end
        Zmkbar(k, :, :) = Zmkbar(k, :, :) + reshape(zmkbar'*zmkbar, [1 N N]);
    end
end
for m = 1:M
    zm = kron(transpose(w(:,m)), H);
    zmbar = zeros(N,N);
    for i = 1:length(idx_vec)
        zmbar(:,i) = zm(:,idx_vec(i));
    end
    Zmbar = Zmbar + zmbar'*zmbar;
end
%% Now, the main optimization loop
% criterion and iterations over the following lines
cvx_begin sdp 
    variable Q2(N,N) hermitian semidefinite
    minimize(real(trace(Zmbar'*Q2)))
    subject to
        for kk = 1:M
            real(trace((reshape(Zmkbar(kk, :, :), [N N]))'*Q2)) >= RF_Thr(kk);
        end
cvx_end
[Q21, Q22] = eig(Q2);
[sorted_eigval,  sorted_eigval_idx] = sort(diag(Q22));
qstar = sqrt(sorted_eigval(end)).*Q21(:, sorted_eigval_idx(end));
qfinal = zeros(N, 1);
for r = 1:length(qfinal)
    [x,y] = find(abs(qstar(r) - lorentz_value) == min(abs(qstar(r) - lorentz_value)));
    qfinal(r) = lorentz_value(x(1));
end
Q = zeros(N, N_d);
for i = 1:N_d
    for l = 1:N_e
        Q((i - 1)*N_e + l, i) = qfinal((i - 1)*N_e + l);
    end
end
Bk = zeros(M, N_d, N_d);
for k = 1:M
    bk = (a(:, k))'*H*Q;
    Bk(k, :, :) = transpose(bk'*bk);
end
c = H*Q;
C = transpose(c'*c);
cvx_begin sdp 
    variable W(N_d,N_d) hermitian semidefinite
    minimize(real(trace(transpose(W)*C)))
    subject to
    for kk = 1:M
        real(trace(transpose(W)*reshape(Bk(kk, :, :), [N_d N_d]))) >= RF_Thr(kk);
    end
cvx_end
[W1, W2] = eig(W);
[sorted_eigval,  sorted_eigval_idx] = sort(diag(W2));
w = zeros(N_d, M);
for m = 1:M
    w(:, m) = sqrt(sorted_eigval(N_d - m + 1)).*W1(:, sorted_eigval_idx(N_d - m + 1));
end
Zmkbar = zeros(M, N, N);
Zmbar = zeros(N, N);
for k = 1:M
    for m = 1:M
        zmk = kron(transpose(w(:,m)), (a(:, k))'*H);
        zmkbar = zeros(1,N);
        for i = 1:length(idx_vec)
            zmkbar(i) = zmk(idx_vec(i));
        end
        Zmkbar(k, :, :) = Zmkbar(k, :, :) + reshape(zmkbar'*zmkbar, [1 N N]);
    end
end
for m = 1:M
    zm = kron(transpose(w(:,m)), H);
    zmbar = zeros(N,N);
    for i = 1:length(idx_vec)
        zmbar(:,i) = zm(:,idx_vec(i));
    end
    Zmbar = Zmbar + zmbar'*zmbar;
end