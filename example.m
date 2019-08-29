% solve equation  D^alpha u = f

alpha = 0.8;
c = 3;
u_fun = @(t) t.^c;
f_fun = @(t) gamma(c+1)/gamma(c+1-alpha) * t.^(c-alpha); % source term

m = 100;
T = 1;
qformula.alpha = alpha;
qformula.w = 1;
t_array = T*(0:m)/m;
u0 = u_fun(0);

tol = 1e-8;
formula = @L1_formula;
fhp = formula(qformula, t_array, u0, tol);

un_pre = u0; % 
u = zeros(m, 1);
for i = 1:m
    % update fhp by add the previous value of u.  un_pre = u(i-1);
    [fhp, uh] = fhp.update(i, un_pre); % According to matlab's syntex, we must return the value of the object here.
    % uh = fhp.get_history_array(i);
    t = fhp.get_t(); % get the approximation point
    b = f_fun(t);
    c0 = fhp.get_wn(i); % get the coefficient of u(i);
    % we have fhp.get_wn(i) * u(i) + uh = f
    u(i) = (b - uh)/c0;
    un_pre = u(i);
end
u_exact = u_fun(t_array(2:end));
plot(t_array(2:end), u, t_array(2:end), u_exact);