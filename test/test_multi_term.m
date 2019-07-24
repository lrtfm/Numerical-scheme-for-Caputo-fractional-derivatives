% Test approximation formula of Caputo fractional derivatives
function ret = test_multi_term(alpha, tol, T)
if nargin < 3
    T = 1;
end
if nargin < 2
    tol = 1e-10;
end
if nargin < 1
    alpha = [0.8, 0.9];
end

w = alpha; % ones(size(alpha));
c_t = 1.8;
% alpha = 0.8;
% tol = 1e-10;
ufun = @(t)t^c_t;
dufun = @(t)sum(w.*(t.^(c_t-alpha)./gamma(c_t + 1 - alpha)*gamma(c_t+1)));
% T = 20;

% first order
formula1 = @L1_formula;
formula2 = @Fast_L1_formula;

% second order
% single term
% formula3 = @L2_1_sigma_single_term;
% formula4 = @Fast_L2_1_sigma_single_term;

% uniform mesh
formula5 = @Fast_L2_1_sigma_uniform;
formula6 = @L1_2_formula_uniform;

k = 2;
f_handle = formula1;

n = 10;
dts = zeros(n, 1);
errs = zeros(n, 1);
for j = 1:n
    m = 2^(j+1);
    t_array = ((0:m)./m).^k*T;
    u0 = ufun(0);
    qformula.alpha = alpha;
    qformula.w = w;
    fhp = f_handle(qformula, t_array, u0, tol);
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
        left(i,:) = [uh, un, fhp.get_wn(i)]; 
    end
    % plot(t_array(2:end), du_exact, '-*', t_array(2:end), du_numerical, '-o')
    % legend('du_exact', 'du_numerical');
    dts(j) = m;
    errs(j) = abs(du_exact(end) - du_numerical(end));
end

rates = -log(errs(1:end-1)./errs(2:end))./log(dts(1:end-1)./dts(2:end));
ret = table(dts, abs(errs), [0; rates], 'VariableNames', {'dt', 'error', 'order'});
ret.Properties.Description =   func2str(f_handle);

