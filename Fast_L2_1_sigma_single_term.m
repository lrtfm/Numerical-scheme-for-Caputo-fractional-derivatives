classdef Fast_L2_1_sigma_single_term
% This is an implementation of fast evaluation of the Caputo fractional
% derivative by FL2-1-sigma method.
%
% Input: 
%   qformula:  qformula.alpha: scalar
%              qformula.w:     scalar
%   t:         partation of [0, T] (can be non uniform)
%   u0:        initial value
%   tol:       EOS tolerence
%
% Ref: Y. Yan et al., Commun. Comput. Phys., 22(2017), pp. 1028-1048
%
%  Author: Zongze Yang
%  Email: yangzongze@gmail.com
%  Data: 2019-07-24
    properties (Access = public)
        alpha
        w
        shape
        u0
        up
        upp
        t
        u_hist
        tau
        tol
        n
        xs
        ws
        Ai
        Bi
        fun
        saveBi
        c
        lambda_n
        sigma
        un_h
    end

    methods
        function obj = Fast_L2_1_sigma_single_term(qformula, t, u0, tol)
            if nargin < 4
                tol = 1e-10;
            end

            obj.tol = tol;
            alp = qformula.alpha;
            obj.alpha = alp;
            obj.w = qformula.w;
            obj.t = t(:)';
            obj.tau = diff(t(:)');
            obj.n = 0;
            obj.u0 = u0(:);
            obj.shape = size(u0);
            obj.un_h = zeros(size(u0));
            len = length(obj.u0);
            
            dt_min = min(obj.tau);
            [tmpxs, tmpws] = SOEappr(alp, tol, dt_min, t(end));
            obj.xs = tmpxs(:);
            obj.ws = tmpws(:)/gamma(1-alp);
            obj.u_hist = zeros(len, length(tmpws));
            sga = alp/2;
            obj.sigma = sga;
            obj.fun = @(sii, c, d, s)(2*s - c).*exp(-sii.*(d - s));
        end
        
        function [obj, ret] = update(obj, n, u_n_minus_one)
            assert(obj.n + 1 == n);
            obj.n = n;
            if n > 1
                obj.upp = obj.up;
            end
            obj.up = u_n_minus_one(:);
            
            % tn_sigma = obj.sigma*obj.t(obj.n) + (1-obj.sigma)*obj.t(obj.n+1);
            obj.lambda_n = (1-obj.sigma)^(1-obj.alpha)/...
                (obj.tau(n)^obj.alpha*gamma(2-obj.alpha));
            if n > 1
                obj.c = exp(-obj.xs*(obj.sigma*obj.tau(n-1) +...
                    (1-obj.sigma)*obj.tau(n)));
                obj.Ai = - exp(-obj.xs*(1-obj.sigma)*obj.tau(n)).* ...
                    (exp(-obj.xs*obj.tau(n-1)).*(2 + ...
                    obj.xs*obj.tau(n-1)*(2+obj.tau(n)/obj.tau(n-1)))...
                    - 2 - obj.tau(n)/obj.tau(n-1)*obj.xs*obj.tau(n-1))./...
                    (obj.xs.^2*obj.tau(n-1)*(obj.tau(n) + obj.tau(n-1)));
               % Ai_ = - arrayfun(@(x)quadgk(@(s)obj.fun(x, obj.t(n+1) + obj.t(n),tn_sigma, s),...
               %      obj.t(n-1), obj.t(n)), obj.xs)/((obj.tau(n)+obj.tau(n-1))*obj.tau(n-1));
                obj.Bi = exp(-obj.xs*(1-obj.sigma)*obj.tau(n)).* ...
                    (exp(-obj.xs*obj.tau(n-1)).*(2 + ...
                    obj.xs*obj.tau(n-1)) - 2 + obj.xs*obj.tau(n-1))./...
                    (obj.xs.^2*obj.tau(n)*(obj.tau(n) + obj.tau(n-1)));
               % Bi_ = arrayfun(@(x)quadgk(@(s)obj.fun(x, obj.t(n-1) + obj.t(n),tn_sigma, s),...
               %      obj.t(n-1), obj.t(n)), obj.xs)/((obj.tau(n)+obj.tau(n-1))*obj.tau(n));
                if n == 2
                    obj.saveBi = zeros(size(obj.Bi));
                end
                tmp = obj.up * obj.saveBi'; % V^n = u_hist + tmp
                tmp2 = obj.up * obj.Bi'; % V^n = u_hist + tmp
                obj.u_hist = (obj.u_hist + tmp).* obj.c' ...
                             + (obj.up - obj.upp) * obj.Ai' - tmp2;
                obj.saveBi = obj.Bi;
            end
            
            ret = obj.u_hist * obj.ws - obj.lambda_n*obj.up;
            ret = obj.w*ret;
            obj.un_h = ret;
        end
        
        function ret = get_t(obj)
            ret = obj.sigma*obj.t(obj.n) + (1-obj.sigma)*obj.t(obj.n+1);
        end
        
        function ret = get_tn(obj)
            ret = obj.t(obj.n + 1);
        end
        
        function ret = get_history_array(obj, n)
            assert(obj.n == n);
            ret = obj.un_h;
        end

        function ret = get_ti(obj, i)
            ret = obj.t(i + 1);
        end
        
        function ret = get_sigma(obj)
            ret = obj.sigma;
        end
            
        function ret = get_wn(obj, n)
            ret = obj.lambda_n;
            if n > 1
                ret = ret + obj.Bi' * obj.ws;
            end
            ret = obj.w*ret;
        end

    end
    
end
