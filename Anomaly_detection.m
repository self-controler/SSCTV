clear;
clc

addpath(genpath('comparison_methods'));
addpath(genpath('dataset'));
addpath(genpath('prox_operator'));
addpath(genpath('tensor_toolbox'));
addpath(genpath('tensor_SVD'));
addpath(genpath('utils'));

Method_list = {'RPCA','CTV','CSSTV'};
Run_RPCA   = 1;
Run_CTV   = 1;
Run_SSCTV   = 1;
en_list = [];

    load Urban.mat
    M=Urban; clear Urban
    M=M(1:80,189:288,:);
%     M=M(1:150,139:288,:);
    [no_lines,no_rows, no_bands]=size(M);
    load UGt.mat
    GT=UGt;
%     GT(70:end, 1:10)=0;
%     GT(70:end,90:end)=0;
    sa=80;
    sb=100;
%     h=100;
%     w=80;  0.8919

%M_3D= M./max(M(:));
ndim = size(M);
M = hyperConvert2D(M);
M=M./repmat(sqrt(sum(M.^2)),size(M,1),1);
M_3D= hyperConvert3D(M,no_lines,no_rows);
GT2=zeros(no_lines, no_rows);
X = M;
GT2(1:sa,end-sb+1:end)=GT;
X=X';
X=reshape(X,no_lines,no_rows,no_bands);

%% RPCA
if Run_RPCA == 1
    lambda = 1/sqrt(no_rows*no_bands);
    [~,E,~] = rpca_m(M');
    E = reshape(E,no_lines,no_rows,no_bands);
    re{1}=reshape(sqrt(sum(E.^2,3)),[no_lines,no_rows]);
    [tpr{1},fpr{1},thresholds1] = roc(GT2(:)',re{1}(:)');
    AUC(1)=trapz(fpr{1},tpr{1});%figure
    fprintf('AUC of RPCA:%f\n',AUC(1));
    en_list = [en_list,1];
end 
%% CTV
if Run_CTV == 1
    opts.lambda = 1/sqrt(no_rows*no_bands);
    [~,E] = ctv_rpca(M_3D,opts);
    E = reshape(E,no_lines,no_rows,no_bands);
    re{2}=reshape(sqrt(sum(E.^2,3)),[no_lines,no_rows]);
    [tpr{2},fpr{2},thresholds2] = roc(GT2(:)',re{2}(:)');
    AUC(2)=trapz(fpr{2},tpr{2});%figure
    fprintf('AUC of CTV-RPCA:%f\n',AUC(2));
    en_list = [en_list,2];
end
%% SSCTV
if Run_SSCTV == 1
opts.lambda = 1/sqrt(no_rows*no_bands);
opts.rho=1.25;%learning rate of ADMM
    [~,E] = csstv_rpca(M_3D,opts);
    E = reshape(E,no_lines,no_rows,no_bands);
    re{3}=reshape(sqrt(sum(E.^2,3)),[no_lines,no_rows]);
    [tpr{3},fpr{3},thresholds3] = roc(GT2(:)',re{3}(:)');
    AUC(3)=trapz(fpr{3},tpr{3});%figure
    fprintf('AUC of SSCTV-RPCA:%f\n',AUC(3));
    en_list = [en_list,3];
end
    




