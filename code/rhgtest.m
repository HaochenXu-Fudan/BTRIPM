function [x,y,hist,time] = rhgtest(x0,y0,F,f,diff,opts)
% check and set essential parameters
if (nargin < 5); error('[x,y,error,hist_error] = rhg(x0,y0,F,f,diff,opts)'); end
if (nargin < 6); opts = []; end
if ~isfield(opts,'s');  opts.s = 0.05;   end
if ~isfield(opts,'s1');  opts.s1 = 0.1;   end
if ~isfield(opts,'s2');  opts.s2 = 0.01;   end
if ~isfield(opts,'K');   opts.K = 500;   end
if ~isfield(opts,'T');   opts.T = 100;   end
if ~isfield(opts,'Ty');    opts.Ty = 100;   end
if ~isfield(opts,'epsilon1');     opts.epsilon1 = 1e-6;   end
if ~isfield(opts,'epsilon2');     opts.epsilon2 = 1e-6;   end
if ~isfield(opts,'epsilon3');     opts.epsilon3 = 1e-6;   end

% copy the parameters from the struct
x = x0; y = y0; T=opts.T;  K=opts.K; Ty=opts.Ty;
s=opts.s;s1=opts.s1;s2=opts.s2;
epsilon1=opts.epsilon1;
ystr=zeros(T+1,length(y));lstr=zeros(T,length(y));
hist = struct();    hist.x = zeros(K,length(x));    hist.y = zeros(K,length(y));
hist.f = zeros(K,1);    hist.F = zeros(K,1);   
tic;

% iterations
for k=1:K
    ystr(1,:)=y;
    for t=1:Ty
        ystr(t+1,:)=ystr(t,:)-s2*diff.fy(x,ystr(t,:)')';
    end
    g=diff.Fx(x,ystr(T+1,:)');
    lstr(T,:)=diff.Fy(x,ystr(T+1,:)')';
    for t=T:-1:2
        g=g+s*(-diff.fxy(x,ystr(t,:)',lstr(t,:)'));
        lstr(t-1,:)=lstr(t,:)-s*diff.fyyv(x,ystr(t,:)',lstr(t,:)')';
    end
    x=x-s1*g;
    y=ystr(T+1,:)';
    %y1=y;
    %for t=1:T
        %y=y-s2*diff.fy(x,y')';
        %if abs(y-y1)<epsilon1
         %   break;
        %end
    %end
    hist.x(k,:) = x';
    hist.y(k,:) = ystr(T+1,:)';
    hist.f(k) = f(x,y)';
    hist.F(k) = F(x,y)';
end
time=toc;
end