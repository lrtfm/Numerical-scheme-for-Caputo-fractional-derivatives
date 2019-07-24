% Test approximation formula of Caputo fractional derivatives
c_t = 3;
alpha = 0.8;
ufun = @(t)t^c_t;
dufun = @(t)t^(c_t-alpha)/gamma(c_t + 1 - alpha)*gamma(c_t+1);
qformula.alpha = alpha;
qformula.w = 1;

% first order
formula1 = @L1_formula;
formula2 = @Fast_L1_formula;

% second order
% single term
formula3 = @L2_1_sigma_single_term;
formula4 = @Fast_L2_1_sigma_single_term;

% uniform mesh
formula5 = @Fast_L2_1_sigma_uniform;
formula6 = @L1_2_formula_uniform;

T = 1;
m = 100;
k = 2; % k = 1 for uniform, k ~= 1 for non uniform
t_array = T*((0:m)/m).^k;
u0 = ufun(0);
fhp = formula1(qformula, t_array, u0); 
un = u0;
du_numerical = zeros(m, 1);
du_exact = zeros(m, 1);
left = zeros(m, 3);
for i = 1:m
    [fhp, uh] = fhp.update(i, un);
    % uh = fhp.get_history_array(i);
    un = ufun(fhp.get_tn());
    du_numerical(i) = fhp.get_wn(i)*un + uh;
    du_exact(i) = dufun(fhp.get_t());
end
plot(t_array(2:end), du_exact, '-*', t_array(2:end), du_numerical, '-o')
legend('du_exact', 'du_numerical');
