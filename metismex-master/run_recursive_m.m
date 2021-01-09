clear;
testcase='m14b';
matfilename=strcat('./',testcase,'.mat');

load(matfilename,'LG','b');
disp('Graph read finished');

M = recursive_m(LG, 3, @isLeafPred);
% tic;[x,flag,relres,iter,RESVEC]=pcg(LG,b,1e-3,1000);toc;
% disp(['MLR: ',num2str(iter)]);
