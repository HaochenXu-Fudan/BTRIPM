%% IMBLO
% This is a file to demo various interior methods on bilevel optimization
% problems.
%% Function Initialization
clear;
experiment = 1; 
switch experiment
    case 1
        funname = '20news';x0 = ones(26214*20,1);y0 = ones(26214*20,1)/20*0.1;
end
[F,f,diff,others] = fun_read(funname);
%% cg
opts = Opts(funname,'cg');
opts.K = 50;
[x_cg,y_cg,hist_cg,time_cg] = cgtest(x0,y0,F,f,diff,opts);
%% rhg & trhg
opts = Opts(funname,'rhg');
opts.K = 50;
[x_rhg,y_rhg,hist_rhg,time_rhg] = rhgtest(x0,y0,F,f,diff,opts);
opts = Opts(funname,'trhg');
opts.K = 25;
[x_trhg,y_trhg,hist_trhg,time_trhg] = rhgtest(x0,y0,F,f,diff,opts);
%% ITD
% opts = Opts(funname,'ITD');
% opts.K = 1;
% [x_ITD,y_ITD,hist_ITD,time_ITD] = ITD(x0,y0,F,f,diff,opts);
%% bvfim
opts = Opts(funname,'bvfim');
opts.K=100;
[x_bvfim,y_bvfim,hist_bvfim,time_bvfim] = bvfim(x0,y0,F,f,diff,opts);
%% Trust region
opts = Opts(funname,'Tregion');
opts.K=10;
[x_Tregion,y_Tregion,hist_Tregion,time_Tregion] = Tregion(x0,y0,F,f,diff,opts);
%% save result
saldir = '.\result\';
savePath = [saldir funname '.mat'];
save(savePath,'hist_Tregion','hist_bvfim','hist_rhg','hist_cg','time_Tregion','time_bvfim','time_rhg','time_cg');
%% plot
plotyes = 1;
if plotyes==1
        F0 = F(x0,y0);
        f0 = f(x0,y0);
        plot_val(funname,hist_bvfim,hist_Tregion,hist_rhg,hist_cg,F0,f0)
end