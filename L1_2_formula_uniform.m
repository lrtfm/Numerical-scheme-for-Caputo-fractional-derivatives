classdef L1_2_formula_uniform
% The L1-2 formula on uniform grid for Caputo fractional derivatives
%
% Input: 
%   qformula:  qformula.alpha: vector
%              qformula.w:     vector
%   t:         partation of [0, T] (must be uniform)
%   u0:        initial value
%   tol:       EOS tolerence
%
% Ref: G. Gao et al., J. Comput. Phys., 259(2014), pp. 33--50.
%
%  Author: Zongze Yang
%  Email: yangzongze@gmail.com
%  Data: 2019-07-24

    properties (Access = private)
        alpha
        shape
        u0
        up
        u_hist
        tol
        N
        tau
        a
        b
        w
        wn
        un_h
        n
    end
    methods
        function obj = L1_2_formula_uniform(qformula, t, u0, tol)
            if nargin < 4
                tol = 1e-10;
            end

            tau_min = min(diff(t));
            tau_max = max(diff(t));
            assert(abs(tau_min - tau_max) < tol);% The partation must be uniform grid.
            tau = tau_min;

            obj.tol = tol;
            obj.alpha = qformula.alpha(:);
            obj.w = qformula.w(:);
            obj.tau = tau;
            obj.N = length(t) - 1;
            obj.n = 0;
            obj.u0 = u0(:);
            obj.up = zeros(length(obj.u0), obj.N + 1);
            obj.shape = size(u0);
            obj.un_h = zeros(size(u0));
            [a_, b_] = obj.init();
            obj.a = a_;
            obj.b = b_;
        end
        
        function [a, b] = init(obj)
            s = 1:obj.N;
            coef = obj.tau.^obj.alpha.*gamma(2-obj.alpha);
            a = s.^(1-obj.alpha) - (s-1).^(1-obj.alpha);
            b = (s.^(2-obj.alpha) - (s-1).^(2-obj.alpha))./(2-obj.alpha) ...
                - 0.5*(s.^(1-obj.alpha) + (s-1).^(1-obj.alpha));
            a(:, 1) = 1; % fix error for alpha = 1;
            b(:, 1) = obj.alpha./(4 - 2*obj.alpha); % fix error for alpha = 1;
            a = a./coef;
            b = b./coef;
        end
        
        function [obj, ret] = update(obj, n, u_n_minus_one)
            assert(obj.n + 1 == n);
            obj.n = n;
            obj.up(:, n) = u_n_minus_one(:);
            c = obj.a(:,1:n);
            if n > 1
                c(:, 1:(n-1)) = c(:, 1:(n-1)) + obj.b(:, 1:(n-1));
                c(:, 2:n) = c(:, 2:n) - obj.b(:, 1:(n-1));
            end
            wc = obj.w'*c;
            obj.wn = wc(1);
            nwc = [diff(wc), -wc(end)];
            nnwc = nwc(end:-1:1);
            ret = obj.up(:, 1:n)*nnwc';
            obj.un_h = reshape(ret, obj.shape);
        end

        function ret = get_t(obj)
            ret = obj.n * obj.tau;
        end
        
        function ret = get_tn(obj)
            ret = obj.n * obj.tau;
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

