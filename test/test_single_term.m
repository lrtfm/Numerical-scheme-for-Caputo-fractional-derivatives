% Test approximation formula of Caputo fractional derivatives
function ret = test_single_term(f, k, nms)
if nargin < 1
    f = 3;
end
if nargin < 2
    k = 1;
end
if nargin < 3
    nms = [20, 40, 80, 160]';
end
nms = nms(:);

c_t = 0.6;
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

if f < 1 || f > 6
    error('f must be in 1 to 6')
end

eval(['formula = formula' int2str(f), ';']);

T = 1;

dts = zeros(size(nms));
errs = zeros(size(nms));
ts = zeros(size(nms));
for j = 1:length(nms)
    m = nms(j);
    % k = 2; % k = 1 for uniform, k ~= 1 for non uniform
    t_array = T*((0:m)/m).^k;
    u0 = ufun(0);
    fhp = formula(qformula, t_array, u0, 10^(-9)); 
    un = u0;
    du_numerical = zeros(m, 1);
    du_exact = zeros(m, 1);
    left = zeros(m, 3);
    s = tic;
    for i = 1:m
        [fhp, uh] = fhp.update(i, un);
        % uh = fhp.get_history_array(i);
        un = ufun(fhp.get_tn());
        du_numerical(i) = fhp.get_wn(i)*un + uh;
        du_exact(i) = dufun(fhp.get_t());
    end
    e = toc(s);
    % plot(t_array(2:end), du_exact, '-*', t_array(2:end), du_numerical, '-o')
    % legend('du_exact', 'du_numerical');
    dts(j) = m;
    errs(j) = du_exact(end) - du_numerical(end);
    ts(j) = e;
end


rates = -ConvergenceRate(dts, abs(errs));
ret = table(dts, abs(errs), rates, ts, 'VariableNames', {'dt', 'error', 'order', 'time'});
end
