function AWR_CMRAC_Results = AWR_CMRAC

% Kim S, Lee H, Cho N, et al. 
% Data-Efficient Active Weighting Algorithm for Composite Adaptive Control 
% Systems[J]. 
% IEEE Transactions on Automatic Control, 2022.

% Cho N, Shin H S, Kim Y, et al. 
% Composite model reference adaptive control with parameter convergence 
% under finite excitation[J]. 
% IEEE Transactions on Automatic Control, 2017, 63(3): 811-818.


% plant dynamics
A = [0, 1; 0, 0];
B = [0; 1];

B_pinv = pinv(B);

% regressor
N = 5;
Psi = @(x)[x(1);
           x(2);
           abs(x(1))*x(2);
           abs(x(2))*x(2);
           x(1)^3];
% piece-wise constant parameters
Theta = @(t) [3.63; -8.58; 20.2; -21.9; -51.88] ...
             + [-22.22; 23.74; -82.66; 31.45; 73.33] * (t>=4 ) ...
             - [-22.22; 23.74; -82.66; 31.45; 73.33] * (t>=13);
         
% reference signal
t_r = [8, 17];
ref = @(t) 1.0*(t>=0 & t<=t_r(1)) + 0.0*(t>t_r(1) & t<=t_r(2)) ...
            + 1.0*(t>t_r(2));


Q_base = diag([2800, 1]);
R_base = 100;
[K, ~, ~] = lqr(A, B, Q_base, R_base);
Ar = A - B * K;
Kr = -5.2888;
Br = -B * Kr;
% Solve P for Ar' * P + P * Ar + Q = 0
Q = diag([100,10]);
P = lyap(Ar', Q);

dt = 1e-4;
t = 0 : dt : 25;

alpha = 10;  % filter constant
beta = 20;
sigma = 4;

lambda1 = 700;
lambda2 = 250;

% constant adaptation gain 
Gama1 = 1000;

Gama2 = zeros(1, length(t));
Gama2(1) = 1;

% system state
xt = zeros(2, length(t));
% reference model
xr = zeros(2, length(t));
% control
ut = zeros(1, length(t));
u_ad = zeros(1, length(t));
% tracking error
er = zeros(2, length(t));

% filter variables
er_f = zeros(2, length(t));
u_ad_f = zeros(1, length(t));
psi_f = zeros(N, length(t));
delta_f = zeros(1, length(t));
w = zeros(1, length(t));
Y = zeros(N, length(t));
U = zeros(N, length(t));
Omega = zeros(1, length(t));

% parameter estimation
theta_est = zeros(N, length(t));

% bound on the maximum eigenvalue of F
bF = 3.5;
% bound on the norm of weighted regressor
bh = 7;

% time step to observe new regressor
Dt = 0.01;

F = zeros(N, N);  
G = zeros(N, 1); 
h = 0;

% initial w0 and w1
w0 = 1 - dt * beta;
w1 = dt;

t_reset = 0;
i_reset = 1;

for i = 1 : length(t) - 1
    % reference system
    xr(:,i+1) = xr(:,i) + dt * (Ar * xr(:,i) + Br * ref(t(i)));
    % resetting filtered variables
    if find(abs(t_r-t(i)) < dt/5)
        t_reset = t(i);
        i_reset = i;
        er_f(:,i) = zeros(2,1);
        u_ad_f(i) = 0;
        psi_f(:,i) = random('uniform', -1e-50, 1e-50, [N, 1]);
        F = random('uniform', -1e-50, 1e-50, [N, N]); 
        G = random('uniform', -1e-50, 1e-50, [N, 1]);
        if t(i) == t_r(2)
            alpha = 1;
        end
    end                 
    er_f(:,i+1) = er_f(:,i) + dt * (-alpha * er_f(:,i) + er(:,i));
    u_ad_f(:,i+1) = u_ad_f(:,i) + dt * (-alpha * u_ad_f(:,i) + u_ad(:,i));
    psi_f(:,i+1) = psi_f(:,i) + dt * (-alpha * psi_f(:,i) + Psi(xt(:,i)));
    delta_f(:,i) = pinv(B) * (er(:,i) - alpha * er_f(:,i) ...
                   - exp(-alpha*(t(i+1) - t_reset)) * er(i_reset) ...
                   - Ar * er_f(:,i) + B * u_ad_f(:,i));
    F = (1 - dt * beta) * F + dt * psi_f(:,i) * psi_f(:,i)';
    G = (1 - dt * beta) * G + dt * psi_f(:,i) * delta_f(i);
    h = (1 - dt * beta) * h + dt * norm(psi_f(:,i+1));
    % find optimal (w0, w1) for every step of length Dt
    if mod(i, Dt/dt) == 0  
        [eig_vectorF, eig_valuesF] = eig(F); 
        % eigenvalues sorted in an increasing order
        [eig_ValF, indsVal] = sort(diag(eig_valuesF));
        eig_ValF(eig_ValF<0) = 0;
        eig_VecF = eig_vectorF(:,indsVal);
        r = length(find(eig_ValF==0));   % multiplicity
        r = r + 1*(r==0);
        % eta projected to the eigen directions
        psi_f_proj = eig_VecF' * psi_f(:,i+1);
        if r > 1
            % eigenvectors sorted using projected eta
            [~,indsVec] = sort(abs(psi_f_proj));
            eig_VecF = eig_VecF(:,indsVec);
            psi_f_proj = eig_VecF' * psi_f(:,i+1);
        end
        % minimum effective eigenvalues of F
        F_eig_min12 = eig_ValF(r:r+1);
        % maximum eigenvalues of F
        F_eig_max12 = eig_ValF(end-1:end);
        % compute optimal p and q
        LF = 1.5 * max(eig_ValF);
        Lh = 1.1 * h;
%         LF = (bF - max(eig_ValF)) / 2 * tanh(0.01*t(i) - 1) + (bF + max(eig_ValF)) / 2;
%         Lh = (bh - h) / 2 * tanh(0.01*t(i) - 1) + (bh + h) / 2;
        [w1, w0, fval] = Optim_w(r, psi_f_proj, F_eig_min12, F_eig_max12, LF, h, Lh);
        F = w0 * F + w1 * psi_f(:,i) * psi_f(:,i)';
        h = w0 * h + abs(w1) * norm(psi_f(:,i+1));
        G = w0 * G + w1 * psi_f(:,i) * delta_f(i);
    end
    % DREM prodecure
    w(i+1) = det(F);
    Y(:,i+1) = w(i+1) * lsqminnorm(F, G);
    U(:,i+1) = U(:,i) + dt * (exp(-sigma * (t(i+1) - t_reset)) * w(i) * Y(:,i));
    Omega(i+1) = Omega(i) + dt * (exp(-sigma * (t(i+1) - t_reset)) * w(i)^2);
    % update law
    theta_est(:,i+1) = theta_est(:,i) + dt * (Gama1 * Psi(xt(:,i)) * er(:,i)' * P * B...
                          + Gama2(i) * Omega(i) * (U(:,i) - Omega(i) * theta_est(:,i)));
    Gama2(i+1) = Gama2(i) + dt * (lambda1 * Gama2(i) - lambda2 * Gama2(i)^2 * Omega(i)^2);
    % adaptive control
    u_ad(:,i+1) = theta_est(:,i+1)' * Psi(xt(:,i));
    ut(:,i+1) = -K * xt(:,i) - Kr * ref(t(i)) - u_ad(:,i+1);
    % plant dynamics
    xt(:,i+1) = xt(:,i) + dt * (A * xt(:,i) + B * (ut(:,i+1) + Theta(t(i))'*Psi(xt(:,i))));
    % tracking error
    er(:,i+1) = xt(:,i+1) - xr(:,i+1);
end

AWR_CMRAC_Results.t = t;
AWR_CMRAC_Results.x = xt;
AWR_CMRAC_Results.u = ut;
AWR_CMRAC_Results.Omega = Omega;
AWR_CMRAC_Results.theta_est = theta_est;
AWR_CMRAC_Results.errs = [er; theta_est - Theta(t)];







































































































