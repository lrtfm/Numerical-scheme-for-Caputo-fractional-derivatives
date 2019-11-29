# Matlab Implementation of Some Numerical Schemes for Caputo Fractional Derivatives
## Scheme List
1. L1 scheme

    Reference:
    + Yanan Zhang, Zhizhong Sun and Honglin Liao. *Finite difference methods 
    for the time fractional diffusion equation on non-uniform meshes.*
    _Journal of Computational Physics_. 265 (2014) 195--210.
    _doi: 10.1016/j.jcp.2014.02.008_

2. Fast L1 scheme

    Reference:
    + Shidong Jiang, Jiwei Zhang, Qian Zhang2 and Zhimin Zhang. *Fast 
    Evaluation of the Caputo Fractional Derivative and its Applications to 
    Fractional Diffusion Equations*.
    _Communications in Computational Physics_. 21 (2017) 650--678.
    _doi: 10.4208/cicp.OA-2016-0136_

3. L1-2 scheme

    Ref: 
    + G. Gao et al., J. Comput. Phys., 259(2014) 33--50.

4. L2-1-<img src="./svgs/8cda31ed38c6d59d14ebefa440099572.svg" align=middle width=9.98290094999999pt height=14.15524440000002pt/> scheme

    Reference: 
    + Anatoly A.Alikhanov. *A new difference scheme for the time fractional diffusion
    equation*. _Journal of Computational Physics_. 280 (2015) 424-438
    _doi: 10.1016/j.jcp.2014.09.031_
    
    + Hong-lin Liao and William McLean and Jiwei Zhang. A Discrete Grönwall Inequality with 
    Applications to Numerical Schemes for Subdiffusion Problems. 
    _SIAM Journal on Numerical Analysis_. 57 (2019) 218-237.
    _doi: 10.1137/16m1175742_

5. Fast L2-1-<img src="./svgs/8cda31ed38c6d59d14ebefa440099572.svg" align=middle width=9.98290094999999pt height=14.15524440000002pt/> scheme

    Reference:
    + Fast Evaluation of the Caputo Fractional Derivative and its Applications to Fractional 
    Diffusion Equations: A Second-Order Scheme. 
    _Communications in Computational Physics_. 22 (04) 1028-1048.
    _doi: 10.4208/cicp.OA-2017-0019_

## Interface of the implementation

### Input
These formulas have same input parameters `qformula`, `t_array`, `u0`, `tol`.
+ `qformula`: This is a structure with two fileds `alpha` and `w`, where `alpha`
  is(are) the fractional order, and `w` is(are) the weight of the fractional term.
  It corresponds to the term <img src="./svgs/3f9faa76d988fc51a7c1cb332b7b841a.svg" align=middle width=34.82299094999999pt height=22.465723500000017pt/> of the equation.
+ `t_array`: This is the temporal mesh <img src="./svgs/12f4e4ef93be91774dee0fb3b3ceef7d.svg" align=middle width=114.04346085pt height=24.65753399999998pt/>.
+ `u0`: This is the initial value, it can be any shape: scale or vector or matrix.
+ `tol`: This parameter only used in fast algorithms. It is used to limit the error 
  of the SOE(sum-of-exponentials) approximation.
  
### Functions
There are some functions of the formulas.
+ `update`: update the data of the pre step by inputing current step number and pre step data in the for loop.
    If the current step number is `i` and the data value of step `i-1` is `pre_u`, 
    we call this function like this `[fhp, uh] = fhp.update(i, pre_u);` where `fhp` is the object handle 
    and `uh` ($u_{hist}$) is the history contribution for the fractional derivative of current step `i`. 
    See the examples for more details. This function must be call step by step from the beginning.
+ `get_wn`: the weight value of current step data. The product of the return value and the current data
    if the contribution of current step data for the fractional derivative, i.e. $D^\alpha_{t^{n-\sigma}} u = w u^n + u_{hist}$,
    where $w$ is the return value of `get_wn` and $u_{hist}$ is the second output of `update`.
+ `get_sigma`: get the sigma value which can be used to obtain the approximation point  <img src="./svgs/e4b4922007e8462a5e0e6b16f169c1e5.svg" align=middle width=32.36645114999999pt height=26.17730939999998pt/>.
    We can obtain it like this $t^{n-\sigma} = \sigma * t^(n-1) + (1 - \sigma) * t^n$;
+ `get_ti`: get the mesh value at step `i`;
+ `get_tn`: get the mesh value at current step.
+ `get_t`: get the approximation point. As the formula approximate the value at <img src="./svgs/46e8ac900863d40b0b44233229c9aaf0.svg" align=middle width=32.36645114999999pt height=26.17730939999998pt/>, it will
    return the value of $t^{n - sigma}$.

This is an example to solve the equation <img src="./svgs/b4931da6a1d4f221e6fe3cc1d3afd857.svg" align=middle width=64.57937639999999pt height=22.831056599999986pt/>, see file `example.m`. 
```  matlab
% solve fractional ordinary differential equation  
%                \[  D^alpha u = f  \]

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

formula1 = @L1_formula;
formula2 = @L1_2_formula_uniform;
formula3 = @L2_1_sigma_single_term;

formula4 = @Fast_L1_formula;
formula5 = @Fast_L2_1_sigma_uniform;
formula6 = @Fast_L2_1_sigma_single_term;

formula = formula6;
fhp = formula(qformula, t_array, u0, tol);

un_pre = u0; % 
u = zeros(m, 1);
for i = 1:m
    % update fhp by add the previous value of u.  un_pre = u(i-1);
    [fhp, uh] = fhp.update(i, un_pre); % According to matlab's syntex, we must return the value of the object here.
                 
    % get the approximation point as follows, or you can get the
    % approximation point just call  t = fhp.get_t();
    sigma = fhp.get_sigma();
    t_i_minus_1 = fhp.get_ti(i-1);  % get t^(i-1)
    t_i = fhp.get_ti(i);            % get t^i
    t = sigma * t_i_minus_1 + (1-sigma) * t_i; % t^{i-sigma} = sigma * t^(i-1) + (1 - sigma) * t^i;
    
    b = f_fun(t);
    c0 = fhp.get_wn(i); % get the coefficient of u(i);
    % we have fhp.get_wn(i) * u(i) + uh = f
    u(i) = (b - uh)/c0;
    un_pre = u(i);
end
u_exact = u_fun(t_array(2:end));
plot(t_array(2:end), u, '-', t_array(2:end), u_exact, 'x');
legend({'Numerical Result', 'Exact Solution'});
```

