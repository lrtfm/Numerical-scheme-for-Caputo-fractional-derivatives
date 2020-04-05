% Compare the L2_1_sigma and fast L2_1_sigma by
% solving fractional ordinary differential equation  
%                \[  D^alpha u = f  \]
% 
alpha = 0.75;
c = 0.9;
u_fun = @(t) t.^c;
f_fun = @(t) gamma(c+1)/gamma(c+1-alpha) * t.^(c-alpha); % source term

u0 = u_fun(0);

tol = 1e-8;

formulas = {@L2_1_sigma_single_term,  @Fast_L2_1_sigma_single_term};
ret = cell(2, 1);
for k = 1:length(formulas)
    formula = formulas{k};

    n = 5;
    errs = zeros(1, n);
    Ms = zeros(1, n);
    ts = zeros(1, n);
    figure
    for j = 1:n
        st = tic;
        m = 10*2.^j;
        Ms(j) = m;
        T = 1;
        qformula.alpha = alpha;
        qformula.w = 1;
        t_array = T*((0:m)/m).^3;
        fhp = formula(qformula, t_array, u0, tol);

        un_pre = u0; % 
        u = zeros(m, 1);
        for i = 1:m
            % update fhp by add the previous value of u.  un_pre = u(i-1);
            [fhp, uh] = fhp.update(i, un_pre); % According to matlab's syntex, we must return the value of the object here.

            % get the approximation point as follows
            t = fhp.get_t();

            b = f_fun(t);
            c0 = fhp.get_wn(i); % get the coefficient of u(i);

            u(i) = (b - uh)/c0;
            un_pre = u(i);
        end
        u_exact = u_fun(t_array(2:end));
        err_ = u(:) - u_exact(:);
        errs(j) = max(abs(err_));
        ts(j) = toc(st);

        plot(t_array(2:end), log10(abs(err_))); hold on
    end
    orders = -log(errs(2:end)./errs(1:end-1))./log(Ms(2:end)./Ms(1:end-1));
    ret{k} = [errs(:)'; 0, orders(:)'; ts(:)'];
end
M = Ms';
L2_1_sigma = array2table(ret{1}', 'VariableNames', { 'Error', 'Order', 'Time'});
Fast_L2_1_sigma = array2table(ret{2}', 'VariableNames', {'Error', 'Order', 'Time'});
result = table(M, L2_1_sigma, Fast_L2_1_sigma);
disp(result);
