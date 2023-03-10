function Glush_CMRAC_Results = Glush_CMRAC

% plant dynamics
A = [0, 1; 0, 0];
B = [0; 1];

% regressor
N = 5;
Psi = @(x)[x(1);
           x(2);
           abs(x(1))*x(2);
           abs(x(2))*x(2);
           x(1)^3];
% piece-wise constant parameters
Theta = @(t) [3.63; -8.58; 20.2; -21.9; -51.88] ...
             + [-22.22; 23.74; -82.66; 31.45; 73.33] * (t>=4) ...
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

alpha = 10;  % filter constant
sigma = 2;
beta = 1;

lambda1 = 1000;
lambda2 = 350;

Gama1 = 1500*eye(N);

dt = 1e-4;
t = 0 : dt : 25;

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
y = zeros(N, length(t));
phi = zeros(N, N, length(t));
w = zeros(1, length(t));
Y = zeros(N, length(t));
U = zeros(N, length(t));
Omega = zeros(1, length(t));

% parameter estimation
theta_est = zeros(N, length(t));

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
        psi_f(:,i) = random('uniform', -1e-30, 1e-30, [N, 1]);
        y(:,i) = random('uniform', -1e-30, 1e-30, [N, 1]); 
        phi(:,:,i) = random('uniform', -1e-30, 1e-30, [N, N]); 
%         U(:,i) = random('uniform', -1e-30, 1e-30, [N, 1]); 
%         Omega(i) = 0;
    end
    er_f(:,i+1) = er_f(:,i) + dt * (-alpha * er_f(:,i) + er(:,i));
    u_ad_f(:,i+1) = u_ad_f(:,i) + dt * (-alpha * u_ad_f(:,i) + u_ad(:,i));
    psi_f(:,i+1) = psi_f(:,i) + dt * (-alpha * psi_f(:,i) + Psi(xt(:,i)));
    delta_f(:,i) = pinv(B) * (er(:,i) - alpha * er_f(:,i) ...
                   - exp(-alpha*(t(i+1) - t_reset)) * er(i_reset) ...
                   - Ar * er_f(:,i) + B * u_ad_f(:,i)); 
%     delta_f(i) = Theta(t(i))' * psi_f(:,i);
    % DREM procedure
    y(:,i+1) = y(:,i) + dt * (-beta * y(:,i) + psi_f(:,i) * delta_f(i));
    phi(:,:,i+1) = phi(:,:,i) + dt * (-beta * phi(:,:,i) + psi_f(:,i) * psi_f(:,i)');
    w(i+1) = det(phi(:,:,i+1));
%     Y(:,i+1) = w(i+1) * Theta(t(i+1));
%     Y(:,i+1) = w(i+1) * (pinv(phi(:,:,i+1)) * y(:,i+1));
    Y(:,i+1) = w(i+1) * lsqminnorm(phi(:,:,i+1), y(:,i+1));
    U(:,i+1) = U(:,i) + dt * (exp(-sigma * (t(i+1) - t_reset)) * w(i) * Y(:,i));
    Omega(i+1) = Omega(i) + dt * (exp(-sigma * (t(i+1) - t_reset)) * w(i)^2);
    % adaptation law
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

Glush_CMRAC_Results.t = t;
Glush_CMRAC_Results.x = xt;
Glush_CMRAC_Results.u = ut;
Glush_CMRAC_Results.theta_est = theta_est;
Glush_CMRAC_Results.errs = [er; theta_est - Theta(t)];