function [w1_opt, w0_opt, maxF] = Optim_w(r, eta_proj, F_eig_min12, F_eig_max12, LF, h, Lh)
% find optimal weight p


F_min1 = F_eig_min12(1);
F_min2 = F_eig_min12(2);

F_max1 = F_eig_max12(2);
F_max2 = F_eig_max12(1);

% eta_proj, F_min1, F_min2, F_max1, F_max2, LF, h, Lh

% % qh as a function of p
% qh = @(p) (Lh - abs(p)*norm(eta_proj)) / h;

% % qF as a function of p
% b_p = @(p) p * F_max2 * norm(eta_proj)^2 ...
%             + p * (F_max1 - F_max2) * norm(eta_proj(1:end-1))^2 ...
%             - LF * (F_max1 - F_max2);
% c_p = @(p) LF^2 - p * LF * norm(eta_proj)^2;
% qF_1 = @(p) (- b_p(p) - sqrt(b_p(p).^2 - 4*F_max1*F_max2 * c_p(p))) /F_max2/F_max1/2;
% qF_2 = @(p) (- b_p(p) + sqrt(b_p(p).^2 - 4*F_max1*F_max2 * c_p(p))) /F_max2/F_max1/2;
% % note q > 0
% qF = @(p) qF_1(p).*(qF_1(p) > 0 & qF_1(p) < qF_2(p)) ...
%           + qF_2(p).*(qF_2(p) > 0 & qF_2(p) < qF_1(p)) ...
%           + qF_1(p).* (qF_1(p) > 0 & qF_2(p) < 0) ...
%           + qF_2(p).* (qF_2(p) > 0 & qF_1(p) < 0);

% % q = @(p)min(qh(p), qF(p))
% q_p = @(p) qh(p).*(qh(p) > 0 & qh(p) < qF(p)) ...
%             + qh(p).*(qh(p) > 0 & qF(p) == 0) ...
%             + qF(p).*(qF(p) < qh(p)) ...
%             + qF(p).*(qh(p) < 0);

% % objective to be maximized
% obj = @(p) -p/2 * norm(eta_proj(r:end))^2 - q_p(p)/2 * (F_min1 + F_min2) ...
%            + 1/2*sqrt((p * norm(eta_proj(r:end))^2 + q_p(p) * (F_min2 - F_min1)).^2 ...
%                 - 4 * eta_proj(r)^2 * (F_min2 - F_min1) * p * q_p(p));
   
% bounds on w1
w1_low = - Lh / norm(eta_proj);
w1_up = - w1_low;

options = optimset('MaxFunEvals',500,'MaxIter',500);

[w1_opt, fval] = fminbnd(@(w1)Obj_w1(w1, r, eta_proj, F_min1, F_min2, ...
                                     F_max1, F_max2, LF, h, Lh), ...
                         w1_low, w1_up, options);

w0_opt = w0_w1(w1_opt, eta_proj, F_max1, F_max2, LF, h, Lh);

maxF = -fval;

end

% w0 as a function of w1
function w0 = w0_w1(w1, eta_proj, F_max1, F_max2, LF, h, Lh)

% qh as a function of p
w0_h = (Lh - abs(w1)*norm(eta_proj)) / h;

% qF as a function of p
if  F_max2 == 0
    w0_F = (w1 * F_max1 * norm(eta_proj(1:end-1))^2 - LF * F_max1) ...
         / (w1 * LF * norm(eta_proj)^2 - LF^2);
else
    b_w1 = w1 * F_max2 * norm(eta_proj)^2 ...
           + w1 * (F_max1 - F_max2) * norm(eta_proj(1:end-1))^2 ...
           - LF * (F_max1 - F_max2);
    c_w1 = LF^2 - w1 * LF * norm(eta_proj)^2;
    w0_F_1 = (- b_w1 - sqrt(b_w1.^2 - 4*F_max1*F_max2 * c_w1)) /F_max2/F_max1/2;
    w0_F_2 = (- b_w1 + sqrt(b_w1.^2 - 4*F_max1*F_max2 * c_w1)) /F_max2/F_max1/2;
    w0_F = w0_F_1.*(w0_F_1 > 0 & w0_F_1 < w0_F_2) ...
          + w0_F_2.*(w0_F_2 > 0 & w0_F_2 < w0_F_1) ...
          + w0_F_1.* (w0_F_1 > 0 & w0_F_2 < 0) ...
          + w0_F_2.* (w0_F_2 > 0 & w0_F_1 < 0);
end

% w0 = @(w1)min(w0_h(w1), w0_F(w1))
w0 = w0_h.*(w0_h > 0 & w0_h < w0_F) ...
    + w0_h.*(w0_h > 0 & w0_F == 0) ...
    + w0_F.*(w0_F < w0_h) ...
    + w0_F.*(w0_h < 0);

end

% objective to be maximized
function obj = Obj_w1(w1, r, eta_proj, F_min1, F_min2, F_max1, F_max2, LF, h, Lh)

w0 = w0_w1(w1, eta_proj, F_max1, F_max2, LF, h, Lh);

obj = -w1/2 * norm(eta_proj(r:end))^2 - w0/2 * (F_min1 + F_min2) ...
      + 1/2 * norm([norm(eta_proj(r+1:end))^2-eta_proj(r)^2, F_min2-F_min1;...
                    2*norm(eta_proj(r+1:end))*norm(eta_proj(r)), 0] * [w1;w0]);
                
end
 


















