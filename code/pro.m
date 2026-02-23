function [x,y,hist,time] = pro(x0,y0,F,f,diff,opts)
% check and set essential parameters
if (nargin < 5); error('[x,y,error,hist_error] = bvfim(x0,y0,F,f,diff,opts)'); end
if (nargin < 6); opts = []; end
if ~isfield(opts,'epsilon1');     opts.epsilon1 = 2e-3;   end
if ~isfield(opts,'epsilon2');     opts.epsilon2 = 5e-2;  end
if ~isfield(opts,'u');     opts.u = @(k) 60/1.04^k;   end
if ~isfield(opts,'tau0');     opts.tau0 = 1;    end
if ~isfield(opts,'eta');     opts.eta = 1.02;    end
if ~isfield(opts,'flag');  opts.flag = 0;    end
if ~isfield(opts,'L');   opts.L = 1;  end
if ~isfield(opts,'K');   opts.K = 500;   end
if ~isfield(opts,'Tz');    opts.Tz = 50;   end
if ~isfield(opts,'s1');  opts.s1 = 0.1;   end
if ~isfield(opts,'s2');  opts.s2 = 0.01;    end
if ~isfield(opts,'alpha');  opts.alpha = 0.01;    end

% copy the parameters from the struct
epsilon1 = opts.epsilon1;   epsilon2 = opts.epsilon2;
u = opts.u; tau0 = opts.tau0;   eta = opts.eta;
L = opts.L; K = opts.K; Tz = opts.Tz;   flag = opts.flag;
s1 = opts.s1;   s2 = opts.s2;   alpha = opts.alpha;

% iterations
x = x0; y = y0; z = y0; tau = tau0; flg = 0; mu = u(1);%% modified 10/28
hist = struct();    hist.x = zeros(K,length(x));    hist.y = zeros(K,length(y));
hist.f = zeros(K,1);    hist.F = zeros(K,1);   
tic;
for k = 1:K
%     mu = u(k);
%     mu is a constant when mu<epsilon1
    if mu<epsilon1||flg==1
%         mu = epsilon1;
        flg = flag;
    else
        mu = u(k);
    end
%     stop if the constraint is violated
    if k>1 && fkl-f(xp,yp)<=0
        fprintf('pro:when k = %d, the constraint is false\n',k);
        x = xp;
        y = yp;
        break;
    end
    for l = 1:L
        tau = tau/eta;
        if tau < epsilon2||flg==1
            tau=tau*eta;
            flg = flag;
        end
        for i=1:Tz
            z = z - s1*(diff.fy(x,z));
        end
        fkl = f(x,z)+mu;
        xp = x;
        yp = y;
        y = yp - s2*(diff.Fy(xp,yp)+tau*diff.fy(xp,yp)/(fkl-f(xp,yp)));
        x = xp - alpha*(diff.Fx(xp,yp)+tau*(diff.fx(xp,yp)-diff.fx(xp,z))/(fkl-f(xp,yp)));
    end   
    hist.x(k,:) = x; hist.y(k,:) = y; hist.f(k) = f(x,y); hist.F(k) = F(x,y);
end
time = toc;
end
