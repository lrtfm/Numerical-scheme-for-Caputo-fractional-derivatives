classdef L2_1_sigma_single_term
% This is an implementation of L2-1-sigma approximation of the Caputo fractional
% derivative on nonuniform grids.
%
% Input: 
%   qformula:  qformula.alpha: scalar
%              qformula.w:     scalar
%   t:         partation of [0, T] (can be non uniform)
%   u0:        initial value
%   tol:       EOS tolerence
%
% Reference: ???
% 
% Author: Zongze Yang
% Email: yangzongze@gmail.com
% Data: 2019-07-24

    properties (Access = public)
        alpha
        w
        N
        shape
        u0
        up
        upp
        t
        tau
        tol
        fun
        n
        wn
        rho
        sigma
        un_h
    end
    methods
        function obj = L2_1_sigma_single_term(qformula, t, u0, tol)
            if nargin < 4
                tol = 1e-10;
            end

            obj.tol = tol;
            alp = qformula.alpha;
            obj.alpha = alp;
            obj.w = qformula.w;
            obj.t = t(:)';
            obj.tau = diff(t(:)');
            obj.N = length(t)-1;
            obj.n = 0;
            obj.u0 = u0(:);
            obj.up = zeros(length(obj.u0), obj.N + 1);
            obj.shape = size(u0);
            obj.un_h = zeros(size(u0));
            
            sga = alp/2;
            obj.sigma = sga;
            obj.rho = obj.tau(1:end-1)./obj.tau(2:end);
            obj.fun = @(a, b, c, s)(s - (a + b)/2).*(c - s).^(-alp);
        end
        
        function [obj, ret] = update(obj, n, u_n_minus_one)
            assert(obj.n + 1 == n);
            obj.n = n;
            obj.up(:, n) = u_n_minus_one(:);
            t_n_sigma = obj.sigma*obj.t(obj.n) + (1-obj.sigma)*obj.t(obj.n+1);
            a = zeros(1, n);
            a(n) = (1-obj.sigma)^(1-obj.alpha)/((1-obj.alpha)*obj.tau(n)^(obj.alpha));
            if n == 1
                c_ = a;
            end
            if n > 1
                a(1:n-1) = (t_n_sigma - obj.t(1:n-1)).^(1-obj.alpha)./((1-obj.alpha)*obj.tau(1:n-1)) -...
                    (t_n_sigma - obj.t(2:n)).^(1-obj.alpha)./((1-obj.alpha)*obj.tau(1:n-1));
                b = arrayfun(@(a,b)quadgk(@(s)obj.fun(a,b,t_n_sigma,s), a, b),...
                    obj.t(1:n-1), obj.t(2:n));
                b = 2*b./(obj.tau(1:n-1).*(obj.tau(1:n-1)+obj.tau(2:n)));
                c_ = a;
                c_(2:end) = c_(2:end) + obj.rho(1:n-1).*b;
                c_(2:end-1) = c_(2:end-1) - b(2:end);
                c_(1) = c_(1) - b(1);
                
            end

            c_ = c_/gamma(1-obj.alpha);
            wc = obj.w*c_;

            obj.wn = wc(end);
            ret = obj.up(:, 2:n)*wc(1:n-1)' - obj.up(:,1:n)*wc';
            obj.un_h = ret;
        end
        
        function ret = get_t(obj)
            ret = obj.sigma*obj.t(obj.n) +...
                (1-obj.sigma)*obj.t(obj.n+1);
        end
        
        function ret = get_tn(obj)
            ret = obj.t(obj.n+1);
        end

        function ret = get_sigma(obj)
            ret = obj.sigma;
        end

        function ret = get_sigma_same(obj)
            ret = obj.sigma;
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
