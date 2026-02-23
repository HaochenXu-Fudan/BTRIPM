function [x,y,hist,time] = augnewton(x0,y0,F,f,diff,opts)
% check and set essential parameters
if (nargin < 5); error('[x,y,error,hist_error] = bvfim(x0,y0,F,f,diff,opts)'); end
if (nargin < 6); opts = []; end
if ~isfield(opts,'epsilon1');     opts.epsilon1 = 1e-5;   end
if ~isfield(opts,'epsilon2');     opts.epsilon2 = 5e-3;  end
if ~isfield(opts,'u');     opts.u = @(k) 4/1.2^k;   end
if ~isfield(opts,'tau0');     opts.tau0 = 1;    end
if ~isfield(opts,'eta');     opts.eta = 1.08;    end
if ~isfield(opts,'flag');  opts.flag = 0;    end
if ~isfield(opts,'L');   opts.L = 1;  end
if ~isfield(opts,'K');   opts.K = 200;   end
if ~isfield(opts,'Tz');    opts.Tz = 50;   end
if ~isfield(opts,'s1');  opts.s1 = 0.1;   end

% copy the parameters from the struct
epsilon1 = opts.epsilon1;   epsilon2 = opts.epsilon2;
u = opts.u; tau0 = opts.tau0;   eta = opts.eta; flag = opts.flag;
L = opts.L; K = opts.K; Tz = opts.Tz;   s1 = opts.s1;  
% iterations
x = x0; y = y0; z = y0; tau = tau0; flg = 0;
hist = struct();    hist.x = zeros(K,length(x));    hist.y = zeros(K,length(y));
hist.f = zeros(K,1);    hist.F = zeros(K,1);   
m = length(x); n = length(y);
tic;
for k = 1:K
    mu = u(k);
%     mu is a constant when mu<epsilon1
    if mu<epsilon1
        mu = epsilon1;
        flg = flag;
    end
%     stop if the constraint is violated
    if k>1 && fkl-f(xp,yp)<=0
        fprintf('augnewton:when k = %d, the constraint is false\n',k);
        x = xp;
        y = yp;
        break;
    end
     for l=1:L     
        tau = tau/eta;
        if tau<epsilon2||flg == 1
            tau = tau*eta;
        end
        for i=1:Tz
            z = z - s1*(diff.fy(x,z));
        end
        xp = x;
        yp = y;
        fkl = f(xp,z)+mu;
        A = zeros(m+2*n+1,m+2*n+1);
        A(1:m,1:m) = diff.Fxx(xp,yp)+tau*(-diff.fxx(xp,z)+diff.fxx(xp,yp))/(fkl-f(xp,yp));
        A(1:m,m+1:m+n) = diff.Fxy(xp,yp)+tau*diff.fxy(xp,yp)/(fkl-f(xp,yp));
        A(1:m,m+n+1) = diff.fx(xp,z)-diff.fx(xp,yp);
        A(1:m,m+n+2:m+2*n+1) = -diff.fxy(xp,z);
        A(m+1:m+n,1:m) = A(1:m,m+1:m+n)';
        A(m+1:m+n,m+1:m+n) = diff.Fyy(xp,yp)+tau*diff.fyy(xp,yp)/(fkl-f(xp,yp));
        A(m+1:m+n,m+n+1) = -diff.fy(xp,yp);
        A(m+n+1,1:m) = -tau/(fkl-f(xp,yp))*A(1:m,m+n+1)';
        A(m+n+1,m+1:m+n) = -tau/(fkl-f(xp,yp))*A(m+1:m+n,m+n+1)';
        A(m+n+1,m+n+1) = fkl - f(xp,yp);
        A(m+n+2:m+2*n+1,1:m) = -tau/(fkl-f(xp,yp))*A(1:m,m+n+2:m+2*n+1)';
        A(m+n+2:m+2*n+1,m+n+2:m+2*n+1) = diff.fyy(xp,z); 
        b1 = diff.Fx(xp,yp)+tau*(diff.fx(xp,yp)-diff.fx(xp,z))/(fkl-f(xp,yp));
        b2 = diff.Fy(xp,yp)+tau*diff.fy(xp,yp)/(fkl-f(x,y));
        b = [b1;b2;zeros(n+1,1)];
        sol = A\b;
        x = xp - sol(1:m);
        y = yp - sol(m+1:m+n);
    end  
    hist.x(k,:) = x; hist.y(k,:) = y; hist.f(k) = f(x,y); hist.F(k) = F(x,y);
end
time = toc;
end
