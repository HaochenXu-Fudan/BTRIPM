function [x,y,hist,time] = cgtest(x0,y0,F,f,diff,opts)
% check and set essential parameters
if (nargin < 5); error('[x,y,error,hist_error] = cg(x0,y0,F,f,diff,opts)'); end
if (nargin < 6); opts = struct(); end
if ~isfield(opts,'s1');  opts.s1 = 0.01;   end
if ~isfield(opts,'s2');  opts.s2 = 0.01;    end
if ~isfield(opts,'K');   opts.K = 1;   end
if ~isfield(opts,'T');    opts.T = 100;   end
if ~isfield(opts,'J');    opts.J = 200;   end
if ~isfield(opts,'epsilon1');     opts.epsilon1 = 1e-4;   end
if ~isfield(opts,'epsilon2');     opts.epsilon2 = 1e-3;   end


% copy the parameters from the struct
x = x0; y = y0; K=opts.K;
T=opts.T;
epsilon1=opts.epsilon1;epsilon2=opts.epsilon2;
s1=opts.s1;s2=opts.s2;J=opts.J;
hist = struct();    hist.x = zeros(K,length(x));    hist.y = zeros(K,length(y));
hist.f = zeros(K,1);    hist.F = zeros(K,1);   
tic;

% iterations
for k = 1:K
    %A=diff.fyy(x,y);
    %b=diff.Fy(x,y);
    %q=cg_linear(A,b,epsilon2);
    y1=y;
    for t=1:T
        y=y-s2*diff.fy(x,y);
        if abs(y-y1)<epsilon2
            break;
        end
    end
    q=cg_linear(diff,epsilon1,x,y,J);
    p=diff.Fx(x,y)-diff.fxy(x,y,q);
    x=x-s1*p;
%     y1=y;
%     for t=1:T
%         y=y-s2*diff.fy(x,y);
%         if abs(y-y1)<epsilon2
%             break;
%         end
%     end
    hist.x(k,:)=x;
    hist.y(k,:) = y;
    hist.f(k) = f(x,y);
    hist.F(k) = F(x,y);
end
time=toc;
end


