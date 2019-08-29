classdef L1_formula
% The L1 approximation of Caputo fractional derivatives
%
% Input: 
%   qformula:  qformula.alpha: vector
%              qformula.w:     vector
%   t:         partation of [0, T] (can be non uniform)
%   u0:        initial value
%   tol:       EOS tolerence
%
%  Reference:
%    Yanan Zhang, Zhizhong Sun, Honglin Liao. *Finite difference methods 
%        for the time fractional diffusion equation on non-uniform meshes.*
%        Journal of Computational Physics. 265 (2014) 195â€?210.
%        _doi: 10.1016/j.jcp.2014.02.008_
%
%  Author: Zongze Yang
%  Email: yangzongze@gmail.com
%  Data: 2019-07-24
%
    properties (Access = private)
        alpha   % a column vector of fractional orders.
        w       % coefficience of the fractional operator.
        N       % the partation number of [0, T].
        t       % a row vector of the partation of [0, T].
                % partation:  t_0,  t_1,  t_2,  ..., t_N.
                % vector t:  [t(1)  t(2)  t(3)  ...  t(N+1)].
        tau     % a vector of interval length of the partation
                % \tau_n = t_{n} - t_{n-1}
                % tau = diff(t);
                % partation:  \tau_1, \tau_2, \tau_3, ..., \tau_N.
                % vector tau: [tau(1) tau(2)  tau(3)  ...  tau(N)].
        shape % the shape of original u, used to reshape the result
                % unh to its original shape.
        u0      % initial value of u
        up      % previous value of u
                % up = [u_0(:), u_1(:), u_2(:), ... ]
        a       % 
        wn      % weight of the u^n.
        un_h    % history part: w * D^alpha u = wn*u^n + un_h
        n       % current time level  n = 0, 1, 2, 3, ...
    end
    methods
        function obj = L1_formula(qformula, t, u0, ~)
            obj.alpha = qformula.alpha(:);
            obj.w = qformula.w(:);
            obj.t = t(:)';
            obj.tau = diff(t);
            obj.N = length(obj.tau);
            obj.shape = size(u0);
            obj.u0 = u0(:);
            obj.up = zeros(length(obj.u0), obj.N + 1);
            % obj.un_h = zeros(size(obj.u0));
            obj.n = 0;
            % obj.wn   % 
            % obj.a    % 
        end
        
        function [obj, ret] = update(obj, n, u_n_minus_one)
            assert(obj.n + 1 == n);
            obj.n = n;  % update the time level
            obj.up(:, n) = u_n_minus_one(:);

            t_n = obj.t(n+1);
            obj.a = ((t_n - obj.t(1:n)).^(1-obj.alpha) ...
                - (t_n - obj.t(2:n+1)).^(1-obj.alpha))...
               ./(obj.tau(1:n).*gamma(2-obj.alpha));  % Equ. (3.1) in Ref 1.
            %  a^n_k -> a(k)
            %
            %  Let 
            %   $ D_t^alpha u = \sum_{k=1}^n c^{(n)}_{n-k} (u^k - u^{k-1}) $
            %  Then  $c^{(n)}_{n-k} = a(k)$.
            %  Also
            %   \[ 
            %     D_t^alpha u = c^{(n)}_0 u^n + 
            %     \sum_{k=1}^{n-1} (c^{(n)}_{n-k} - c^{(n)}_{n-k-1}) u^k 
            %     - c^{(n)}_{n-1} u^{0}) 
            %   \]
            %  If we set c = a, then
            %   c_{n-1}, c_{n-2}, c_{n-3}      c_{0}
            %  [ c(1),    c(2),    c(3),   ..., c(n) ]
            %  Also we have
            %  \[
            %  D_t^alpha u 
            %   =  c(n)*u(n) + \sum_{k=1}^{n-1}(c(k)-c(k+1))*u(k) - c(1)*u(0) 
            %   =  -c(1)*u(0) + \sum_{k=1}^{n-1}(c(k)-c(k+1))*u(k) + c(n)*u(n)
            %  \] 
            c = obj.a;

            wc = obj.w'*c;
            obj.wn = wc(end);
            new_wc = [-wc(1), -diff(wc)];
            ret = reshape(obj.up(:, 1:n)*new_wc', obj.shape);
            obj.un_h = ret;
        end

        function ret = get_t(obj)
            ret = obj.t(obj.n+1);
        end
        
        function ret = get_tn(obj)
            ret = obj.t(obj.n+1);
        end
        
        function ret = get_sigma(obj)
            ret = 0;
        end

        function ret = get_sigma_same(obj)
            ret = 0;
        end
        
        function ret = get_history_array(obj, n)
            assert(obj.n == n);
            ret = obj.un_h;
        end
            
        function ret = get_wn(obj, n)
            assert(obj.n == n);
            ret = obj.wn;
        end
        
    end
    
end
