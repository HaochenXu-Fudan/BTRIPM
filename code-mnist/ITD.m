function [w,lbd,hist,time] = ITD(w0,lbd0,F,f,diff,opts)
% check and set essential parameters
if (nargin < 5); error('[x,y,error,hist_error] = cg(x0,y0,F,f,diff,opts)'); end
if (nargin < 6); opts=struct(); end
if ~isfield(opts,'s1');  opts.s1 = 0.01;   end
if ~isfield(opts,'J');    opts.J = 200;   end
if ~isfield(opts,'epsilon1');     opts.epsilon1 = 1e-4;   end
if ~isfield(opts,'alpha'); opts.alpha=0.01; end
if ~isfield(opts,'t'); opts.t=10; end
if ~isfield(opts,'K'); opts.K=50; end
% if ~isfield(opts,'s2');  opts.s2 = 0.01;    end
% if ~isfield(opts,'K');   opts.K = 1;   end
% if ~isfield(opts,'T');    opts.T = 100;   end
% if ~isfield(opts,'epsilon2');     opts.epsilon2 = 1e-3;   end

w=w0; lbd=lbd0; 
t=opts.t; epsilon1=opts.epsilon1; J=opts.J; s1=opts.s1; K=opts.K;
alpha=opts.alpha;
hist = struct();    hist.x = zeros(K,length(lbd));    hist.y = zeros(K,length(w));
hist.f = zeros(K,1);    hist.F = zeros(K,1);   
tic;

for k=1:K
    for i=1:t
        w=w-alpha*diff.fy(w,lbd);
    end
%     wplbd=-diff.fyy(w,lbd)*diff.fxy(w,lbd);
    q=sd_solve(diff,epsilon1,lbd,w,J);
    p=diff.Fx(w,lbd)-q'*diff.Fy(w,lbd);
    lbd=lbd-s1*p;
    hist.x(k,:)=lbd;
    hist.y(k,:)=w;
    hist.f(k)=f(lbd,w);
    hist.F(k)=F(lbd,w);
end
time=toc;
end


