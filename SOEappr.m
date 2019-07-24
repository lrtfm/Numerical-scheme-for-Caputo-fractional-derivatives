function [xs, ws, nexp] = SOEappr(beta, reps, dt, Tfinal)
% Author: Shidong Jiang, Department of Mathematical Sciences
% New Jersey Institute of Technology, Newark, NJ
% Email: shidong.jiang@njit.edu
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Please kindly cite the paper:
%
% Shidong Jiang, Jiwei Zhang, Qian Zhang and Zhimin Zhang.
% Fast Evaluation of the Caputo Fractional Derivative
% and its Applications to Fractional Diffusion Equations.
% Commun. Comput. Phys. Vol. 21, No. 3, pp. 650-678, 2017.
%
%
% Change Log:
% [2018-11-10] by Zongze Yang (yangzongze@gmail.com)
%   1) Change code format. 
%   2) Correct some syntax warning.
%
%
% This function return sum-of-exponentials approximation for 1/t^beta for 
% t in [dt, Tfinal] with relative error bounded by reps, that is,
% 
% 1/t^beta = \sum_{i=1}^\nexp ws(i) e^{-xs(i) t}
%
% Output parameters:
% xs : nodes
% ws : weights
% nexp : number of exponentials
%
% Input parameters:
% beta : the power of the power function 1/t^beta
% reps : desired relative error
% [dt, Tfinal] : the interval on which the power function is approximated.

    delta = dt/Tfinal;

    h = 2*pi/(log(3)+beta*log(1/cos(1))+log(1/reps));

    tlower = 1/beta*log(reps*gamma(1+beta));

    if beta>=1
        tupper = log(1/delta)+log(log(1/reps))+log(beta)+1/2;
    else
        tupper = log(1/delta)+log(log(1/reps));
    end

    M = floor(tlower/h);
    N = ceil(tupper/h);

    n1 = M:-1;
    xs1 = -exp(h*n1);
    ws1 = h/gamma(beta)*exp(beta*h*n1);
    [ws1new,xs1new] = prony(xs1,ws1);

    n2= 0:N;
    xs2 = -exp(h*n2);
    ws2 = h/gamma(beta)*exp(beta*h*n2);

    xs = [-real(xs1new); -real(xs2.')];
    ws = [real(ws1new); real(ws2.')];

    xs = xs/Tfinal;
    ws = ws/Tfinal^beta;

    nexp = length(ws);

    if 1==2 % change it to 1==2 if the testing is not needed
        % test the accuracy on the interval [dt, Tfinal]
        m = 10000;

        estart = log10(dt);
        eend = log10(Tfinal);
        texp = linspace(estart,eend,m);
        t = 10.^texp;

        ftrue = 1./t.^beta;
        fcomp = zeros(size(ftrue));

        for i = 1:m
            fcomp(i) = sum(ws.*exp(-t(i)*xs));
        end

        fcomp = real(fcomp);

        rerr = norm((ftrue-fcomp)./ftrue,Inf);
        fprintf(1, 'The actual relative L_inf error is %0.5g\n',rerr);
    end

    return

end
%
%
%
function [wsnew, xsnew] = prony(xs,ws)
    M = length(xs);
    errbnd = 1d-12;

    h = xs.^((0:2*M-1)')*ws';
    C=h(1:M);
    R=h(M:2*M-1);
    H=hankel(C,R);

    b=-h;

    q = myls2(H, b, errbnd);

    r = length(q);

    Coef = [1; flipud(q)];

    xsnew=roots(Coef);

    A = xsnew'.^((0:2*M-1)');
    wsnew = myls(A,h,errbnd);

    ind = find(real(xsnew)>=0);
    p = length(ind);
    assert(sum(abs(wsnew(ind))<1d-15) == p)

    ind = find(real(xsnew)<0);
    xsnew = xsnew(ind);
    wsnew = wsnew(ind);

end
%
%
%
function x=myls2(A,b,tol)
% solve the rank deficient least squares problem by SVD
% x is the LS solution, res is the residue

    [m,~]=size(A);
    [Q,R]=qr(A,0);

    if nargin < 3
        tol=1e-13;
    end
    s=diag(R);
    r=sum(abs(s)>tol);

    Q = Q(:, 1:r);
    R = R(1:r,1:r);
    b1 = b(r+1:m+r);

    x= R\(Q.'*b1);

end
%
%
%
function [x,res]=myls(A,b,eps)
% solve the rank deficient least squares problem by SVD
% x is the LS solution, res is the residue
    [~,n]=size(A);
    [U,S,V]=svd(A,0);

    if nargin < 3
        eps=1e-12;
    end
    s=diag(S);
    r=sum(s>eps);

    x=zeros(n,1);
    for i=1:r
        x=x+(U(:,i)'*b)/s(i)*V(:,i);
    end

    if (nargout>1)
        res = norm(A*x-b)/norm(b);
    end

end