classdef Fast_L1_formula
% The fast L1 approximation of Caputo fractional derivatives
% Input: 
%   qformula:  qformula.alpha: vector
%              qformula.w:     vector
%   t:         partation of [0, T] (can be non uniform)
%   u0:        initial value
%   tol:       EOS tolerence
%
%  Reference:
%    + Shidong Jiang, Jiwei Zhang, Qian Zhang2 and Zhimin Zhang. *Fast 
%        Evaluation of the Caputo Fractional Derivative and its Applications to 
%        Fractional Diffusion Equations*.
%        _Communications in Computational Physics_. 21 (2017) 650-678.
%        _doi: 10.4208/cicp.OA-2016-0136_
%
%  Author: Zongze Yang
%  Email: yangzongze@gmail.com
%  Data: 2019-07-24

    properties (Access = private)
        alpha
        w
        t
        tau
        shape
        u0
        up     % u^{n-1}
        upp    % u^{n-2}
        u_hist % saved history part to used in next time level
        un_h
        xs     % SOE approximation points
        ws     % SOE approximation weights
        n
        tol
    end
    methods
        function obj = Fast_L1_formula(qformula, t, u0, tol)
            if nargin < 4 
                tol = 1e-10;
            end
            obj.tol = tol;
            obj.alpha = qformula.alpha(:);
            obj.w = qformula.w(:);
            obj.t = t(:)';
            obj.tau = diff(obj.t);
            obj.n = 0;
            obj.shape = size(u0);
            obj.u0 = u0(:);

            len = length(obj.u0);
            m = length(obj.alpha);
            obj.xs = cell(m, 1);
            obj.ws = cell(m, 1);
            obj.u_hist = cell(m, 1);
            dt_min = min(obj.tau);
            for i = 1:m
                [tmpxs, tmpws] = SOEappr(obj.alpha(i) + 1, tol, dt_min, t(end));
                obj.xs{i} = tmpxs(:);
                obj.ws{i} = tmpws(:);
                obj.u_hist{i} = zeros(len, length(tmpws));
            end
            
        end
        
        function [obj, ret] = update(obj, n, u_n_minus_one)
            assert(obj.n + 1 == n);
            obj.n = n;
            if n > 1
                obj.upp = obj.up;
            end
            obj.up = u_n_minus_one(:);
            ret = zeros(size(obj.u0));
            for i = 1:length(obj.alpha)
                xs_ = obj.xs{i};
                ws_ = obj.ws{i};
                
                if n > 1
                    c = exp(-xs_*obj.tau(n))./(xs_.^2 * obj.tau(n-1));
                    tmp = xs_ * obj.tau(n-1);
                    a = exp(-tmp) - 1 + tmp;
                    b = 1 - exp(-tmp).*(1 + tmp);
                    d = exp(-xs_*obj.tau(n));
                    obj.u_hist{i} = obj.u_hist{i}*spdiags(d, 0, length(d), length(d)) + ...
                        obj.up * (c.*a)' + obj.upp * (c.*b)';
                end
                ret_tmp = obj.alpha(i)*obj.u_hist{i}*ws_ + ...
                    obj.u0/obj.t(n+1)^obj.alpha(i);
                ret_tmp = (1 - obj.alpha(i)) * ret_tmp + ...
                    obj.alpha(i)/obj.tau(n)^obj.alpha(i) * obj.up;
                ret = ret + obj.w(i) * ret_tmp/gamma(2-obj.alpha(i));
            end
            ret = reshape(ret, obj.shape);
            ret = - ret;   % make $ D^\alpha u = wn*u^n + ret $
            obj.un_h = ret;
        end
        
        function ret = get_t(obj)
            ret = obj.t(obj.n+1);
        end
        
        function ret = get_tn(obj)
            ret = obj.t(obj.n + 1);
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
            ret = obj.w'*(obj.tau(n).^(- obj.alpha)./gamma(2 - obj.alpha));
        end
        
        
    end
    
end
