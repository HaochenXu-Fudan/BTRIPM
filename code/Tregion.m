function [x,y,hist,time] = Tregion(x0,y0,F,f,diff,opts)
% check and set essential parameters
if (nargin < 5); error('[x,y,error,hist_error] = bvfim(x0,y0,F,f,diff,opts)'); end
if (nargin < 6); opts = []; end
if ~isfield(opts,'epsilon1');     opts.epsilon1 = 1e-8;   end
if ~isfield(opts,'epsilon2');     opts.epsilon2 = 1e-4;  end
if ~isfield(opts,'u');     opts.u = @(k) 2/1.4^k;   end
if ~isfield(opts,'tau0');     opts.tau0 = 1;    end
if ~isfield(opts,'eta');     opts.eta = 1.02;    end
if ~isfield(opts,'flag');  opts.flag = 1;    end
if ~isfield(opts,'L');   opts.L = 1;  end
if ~isfield(opts,'K');   opts.K = 100;   end
if ~isfield(opts,'Tz');    opts.Tz = 50;   end
if ~isfield(opts,'s1');  opts.s1 = 0.1;   end
if ~isfield(opts,'maxit');  opts.maxit = 5;   end
if ~isfield(opts,'hessian');  opts.hessian = 'accurate';   end
if ~isfield(opts,'zj');  opts.zj = 0;   end

% copy the parameters from the struct
epsilon1 = opts.epsilon1;   epsilon2 = opts.epsilon2; hessian = opts.hessian;
u = opts.u; tau0 = opts.tau0;   eta = opts.eta; flag = opts.flag; zj=opts.zj;
L = opts.L; K = opts.K; Tz = opts.Tz;   s1 = opts.s1;  maxit = opts.maxit;
ini_Delta=0; %decide whether initialize Delta by hand.
if isfield(opts,'Delta'); Delta=opts.Delta; ini_Delta=1; end

% iterations
m = length(x0); n = length(y0);
x = x0; y = y0; z = y0; tau = tau0; flg = 0; mu = u(1);
hist = struct();    hist.x = zeros(K,m);    hist.y = zeros(K,n);
hist.f = zeros(K,1);    hist.F = zeros(K,1);   
tic;
for k = 1:K
%     mu = u(k);
%     mu is a constant when mu<epsilon1
    if mu<epsilon1||flg==1
%          mu=epsilon1;
        flg = flag;
    else
         mu = u(k);
    end
    if zj
        mu = f(x,y);
    end
%     stop if the constraint is violated
    if k>1 && fkl-f(xp,yp)<=0
        fprintf('Tregion:when k = %d, the constraint is false\n',k);
        x = xp;
        y = yp;
        break;
    end
     for l=1:L     
        tau = tau/eta;
        if tau<epsilon2||flg == 1
            tau = tau*eta;
            flg=flag;
        end
        for i=1:Tz
            z = z - s1*(diff.fy(x,z));
        end
        fkl = f(x,z)+mu; 
        xp = x;
        yp = y;
        opts = struct();
        opts.Delta = 2;
%         opts.maxit = 4;
        opts.maxit = maxit;
        opts.m = m;
        opts.n = n;
        opts.difffy = diff.fy;
        opts.Tz = Tz;
        opts.s1 = s1;
        opts.verbose = 0;
        opts.record = 0;
        if ini_Delta 
            opts.Delta = Delta; 
        end
        fun = @(w) prefun(w,F,f,diff,tau,mu,Tz,s1,m,n);
%         hess = @(w,z,v) prehess(w,f,diff,tau,mu,v,Tz,s1,m,n);
        hess = @(w,z,v) prehess(w,z,f,diff,tau,mu,v,m,n,hessian);
        [w,~] = fminTR([x;y],fun,hess,opts);
        x = w(1:m);
        y = w(m+1:m+n);
    end  
    hist.x(k,:) = x; hist.y(k,:) = y; hist.f(k) = f(x,y); hist.F(k) = F(x,y);
end
time = toc;
end
function [g,grad]=prefun(w,F,f,diff,tau,mu,Tz,s,m,n)
x = w(1:m);
y = w(m+1:m+n);
z = y;
for i=1:Tz
    z = z - s*(diff.fy(x,z));
end
fkl = f(x,z)+mu; 
g = F(x,y)-tau*log(fkl-f(x,y));
grad = [diff.Fx(x,y);diff.Fy(x,y)]-tau*[diff.fx(x,z)-diff.fx(x,y);-diff.fy(x,y)]/(fkl-f(x,y));
end
function H = prehess(w,z,f,diff,tau,mu,v,m,n,hessian)
x = w(1:m);
y = w(m+1:m+n);
% here we need to improve
% HF = [diff.Fxx(x,y),diff.Fxy(x,y);diff.Fxy(x,y)',diff.Fyy(x,y)]*v;
% Hf1 = -tau/(fkl-f(x,y))*[diff.fxx(x,z)-diff.fxy(x,z)/diff.fyy(x,z)*diff.fxy(x,z)'-diff.fxx(x,y),-diff.fxy(x,y);-diff.fxy(x,y)',-diff.fyy(x,y)]*v;
% Hf2 = tau/(fkl-f(x,y))^2*[(diff.fx(x,z)-diff.fx(x,y))*(diff.fx(x,z)-diff.fx(x,y))',(diff.fx(x,y)-diff.fx(x,z))*diff.fy(x,y)';diff.fy(x,y)*(diff.fx(x,y)-diff.fx(x,z))',diff.fy(x,y)*diff.fy(x,y)']*v;
% H = HF+Hf1+Hf2;
% to avoid repeadly calculate the hard.
% fxy_xz = diff.fxy(x,z);
% 
% fyx_xyv1 = diff.fyx(x,y,v1);
v1 = v(1:m); v2 = v(m+1:m+n);
if strcmp(hessian,'accurate')
    fkl = f(x,z)+mu; 
    f_xy = f(x,y);
    fxfx = -diff.fx(x,z)+diff.fx(x,y);
    fy_xy = diff.fy(x,y);
    HF = [diff.Fxx(x,y,v1)+diff.Fxy(x,y,v2);diff.Fyx(x,y,v1)+diff.Fyy(x,y,v2)];
    Hf1 = -tau/(fkl-f_xy)*[diff.fxx(x,z,v1)-diff.fxx(x,y,v1)-diff.fxy(x,y,v2);-diff.fyx(x,y,v1)-diff.fyyv(x,y,v2)];
    % Hf11 = -tau/(fkl-f_xy)*[-fxy_xz*(diff.fyy(x,z)\(fxy_xz'*v(1:m)));zeros(n,1)];
    Hf11 = -tau/(fkl-f_xy)*[-diff.fxy(x,z,diff.fyy(x,z)\diff.fyx(x,z,v1));zeros(n,1)];
    % Hf2 = tau/(fkl-f_xy)^2*[fxfx*fxfx',fxfy;fxfy',fy_xy*fy_xy']*v;
    fxfxv = fxfx'*v1;
    fyv = fy_xy'*v2;
    Hf2 = tau/(fkl-f_xy)^2*[fxfx*(fxfxv+fyv);fy_xy*(fxfxv+fyv)];
    H = HF+Hf1+Hf11+Hf2;
elseif strcmp(hessian,'nonaccurate')
    HF = [diff.Fxx(x,y,v1)+diff.Fxy(x,y,v2);diff.Fyx(x,y,v1)+diff.Fyy(x,y,v2)];
    H = HF;
elseif strcmp(hessian,'guassian')
     fkl = f(x,z)+mu; 
    f_xy = f(x,y);
    fxfx = -diff.fx(x,z)+diff.fx(x,y);
    fy_xy = diff.fy(x,y);
     fxfxv = fxfx'*v1;
    fyv = fy_xy'*v2;
    Hf2 = tau/(fkl-f_xy)^2*[fxfx*(fxfxv+fyv);fy_xy*(fxfxv+fyv)];
    HF = [diff.Fxx(x,y,v1)+diff.Fxy(x,y,v2);diff.Fyx(x,y,v1)+diff.Fyy(x,y,v2)];
     H = HF+Hf2;
end
end