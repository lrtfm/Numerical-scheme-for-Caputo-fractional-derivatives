classdef Fast_L2_1_sigma_uniform
% This is an implementation of fast evaluation of the Caputo fractional
% derivative by FL2-1-sigma method for multi-term time fpdes.
%
% Input: 
%   qformula:  qformula.alpha: vector
%              qformula.w:     vector
%   t:         partation of [0, T] (must be uniform)
%   u0:        initial value
%   tol:       EOS tolerence
%
% Reference
%   1. Y. Yan et al., Commun. Comput. Phys., 22(2017), pp. 1028-1048
%   2. G. Gao et al., J. Sci. Comput., 73(2017), pp. 93-121
%
% Author: Zongze Yang
% Email: yangzongze@gmail.com
% Data: 2019-07-24
%
    properties (Access = public)
        alpha
        w
        tau
        T
        N
        shape
        u0
        up
        upp
        u_hist
        tol
        n
        xs
        ws
        Ai
        Bi
        c
        lambda_
        sigma
        un_h
    end
    methods
        function obj = Fast_L2_1_sigma_uniform(qformula, t, u0, tol)
            if nargin < 4
                tol = 1e-10;
            end

            tau_min = min(diff(t));
            tau_max = max(diff(t));
            assert(abs(tau_min - tau_max) < tol);% The partation must be uniform grid.
            tau_ = tau_min;

            obj.tol = tol;
            obj.alpha = qformula.alpha(:);
            obj.w = qformula.w(:);
            obj.tau = tau_;
            obj.N = length(t)-1;
            obj.n = 0;
            obj.u0 = u0(:);
            obj.shape = size(u0);
            obj.un_h = zeros(size(u0));
            len = length(obj.u0);
            
            obj.T = obj.tau*obj.N;
            obj.sigma = obj.init_sigma();
            m = length(obj.alpha);
            obj.Ai = cell(m, 1);
            obj.Bi = cell(m, 1);
            obj.c = cell(m, 1);
            obj.xs = cell(m, 1);
            obj.ws = cell(m, 1);
            obj.u_hist = cell(m, 1);
            sigma_ = obj.sigma;
            obj.lambda_ = zeros(m, 1);
            for i = 1:m
                alpha_ = obj.alpha(i);
                [tmpxs, tmpws] = SOEappr(alpha_, tol, obj.tau, obj.T);
                obj.xs{i} = tmpxs(:);
                obj.ws{i} = tmpws(:)/gamma(1-alpha_);
                obj.u_hist{i} = zeros(len, length(tmpws));

                w_ = obj.xs{i}*tau_;
                obj.lambda_(i) = sigma_^(1-alpha_)/(tau_^alpha_*gamma(2-alpha_));
                Ai_ = (exp(-w_*sigma_).*(2+w_) - exp(-w_*(1+sigma_)).*(2+3*w_))./(2*w_.^2);
                Bi_ = (exp(-w_*sigma_).*(-2+w_) + exp(-w_*(1+sigma_)).*(2+w_))./(2*w_.^2);
                obj.Ai{i} = Ai_';
                obj.Bi{i} = Bi_';
                obj.c{i} = exp(-w_)';
            end
        end
        
        % The sigma in the reference paper is used as $D^{\alpha} u^{n- (1 - \sigma)}$
        % Here, we return the value (1- \sigma), which
        % is used in Fast_L2_1_sigma_single_term, to get a uniform interface.
        % So, the return value sigma is used as $D^{\alpha} u^{n - \sigma}$
        function sigma = get_sigma(obj)
            sigma = 1 - obj.sigma;
        end
        
        function sigma = init_sigma(obj)
            sigma = 0;
            sigma_n = 1 - min(obj.alpha)/2;
            while abs(sigma_n - sigma) > obj.tol
                sigma = sigma_n;
                [f, df] = get_value_F(obj, sigma);
                sigma_n = sigma - f/df;
            end
            sigma = sigma_n;
        end

        function [f, df] = get_value_F(obj, sigma)
            alpha_ = obj.alpha;
            w_ = obj.w;
            tau_ = obj.tau;
            f = w_.*sigma.^(1-alpha_).*(sigma - (1 - alpha_/2)) ...
                .*tau_.^(2-alpha_)./gamma(3-alpha_);
            f = sum(f);
            df =  w_.*sigma.^(-alpha_).*(sigma - (1 - alpha_)/2) ...
                .*tau_.^(2-alpha_)./gamma(2-alpha_);
            df = sum(df);
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
                if n > 1
                    tmp = obj.up * obj.Bi{i}; % V^n = u_hist + tmp
                    obj.u_hist{i} = (obj.u_hist{i} + tmp).* obj.c{i} ...
                                 + (obj.up - obj.upp) * obj.Ai{i} - tmp;
                end
            
                ret_tmp = obj.u_hist{i} * obj.ws{i} - obj.lambda_(i)*obj.up;
                ret = ret + obj.w(i)*ret_tmp;
            end
            obj.un_h = ret;
        end
        
        function ret = get_t(obj)
            ret = (obj.n - 1 + obj.sigma)*obj.tau;
        end
        
        function ret = get_tn(obj)
            ret = obj.tau*obj.n;
        end

        function ret = get_ti(obj, i)
            ret = obj.tau*i;
        end
        
        function ret = get_history_array(obj, n)
            assert(obj.n == n);
            ret = obj.un_h;
        end
            
        function ret = get_wn(obj, n)
            ret = obj.lambda_'*obj.w;
            if n > 1
                for i = 1:length(obj.alpha)
                    ret = ret + obj.w(i)*(obj.Bi{i} * obj.ws{i});
                end
            end
            % ret = obj.coeff*ret;
        end

    end
    
end
