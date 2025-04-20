clear;
clc;
addpath(".\ddebiftool\")
% addpath(".\ddebiftool\hom_dem\")
% addpath(".\ddebiftool\sd_demo\")

fid=fopen('.\origBifurJ2tau0p5c.txt','wt');

for i=1:1
stst.kind='stst';
stst.parameter=[0.05 2 1.0005 0.5*i];
stst.x=[1 -0.6667]';
method=df_mthod('stst');
[stst,success]=p_correc(stst,[],[],method.point)
method.stability.minimal_real_part=-22;
stst.stability=p_stabil(stst,method.stability)
figure(2); clf;
p_splot(stst);

branch1=df_brnch(3,'stst')
branch1.parameter
branch1.parameter.min_bound
branch1.parameter.min_bound(1,:)=[3 0];
branch1.parameter.max_bound(1,:)=[3 2];
branch1.parameter.max_step(1,:)=[3 0.004];

branch1.point=stst;
stst.parameter(3)=stst.parameter(3)+0.001;
[stst,success]=p_correc(stst,[],[],method.point)
branch1.point(2)=stst;
branch1.method.continuation.plot=0;
[branch1,s,f,r]=br_contn(branch1,600)
branch1=br_rvers(branch1);
[branch1,s,f,r]=br_contn(branch1,600)
branch1.method.stability.minimal_real_part=-22;
branch1=br_stabl(branch1,0,0);
[xm,ym]=df_measr(1,branch1);
figure(3); clf;
br_plot(branch1,xm,ym,'b');
plot([0 2],[0 0],'-.');
axis([0 2 -0.5 0.5]);

ym.subfield='l0';
figure(4);clf;
br_plot(branch1,[],ym,'b');
br_plot(branch1,[],ym,'b.');
plot([0,1200],[0,0],'-.')

  for i = 1:length(branch1.point)-1
    if (max(real(branch1.point(i).stability.l0))*max(real(branch1.point(i+1).stability.l0)))<0
        break;
    end
    numcritical = i+1;
  end

hopf=p_tohopf(branch1.point(numcritical));
method=df_mthod('hopf');
method.stability.minimal_real_part=-22;
[hopf,success]=p_correc(hopf,3,[],method.point)

first_hopf=hopf;
hopf.stability=p_stabil(hopf,method.stability);
figure(5); clf;
p_splot(hopf);
branch2=df_brnch([4 3],'hopf');
branch2.parameter.min_bound(1:2,:)=[[4 0.0001]' [3 -1.5]']';
branch2.parameter.max_bound(1:2,:)=[[4 20]' [3 3.5]']';
branch2.parameter.max_step(1:2,:)=[[4 0.05]' [3 0.005]']';
branch2.point=hopf;
hopf.parameter(4)=hopf.parameter(4)+0.005;
[hopf,success]=p_correc(hopf,3,[],method.point)

% branch2.point(2)=hopf;
% figure(6); clf;
% [branch2,s,f,r]=br_contn(branch2,1000)
% branch2=br_rvers(branch2);
% [branch2,s,f,r]=br_contn(branch2,950)
% branch2=br_stabl(branch2,0,0);
% figure(7); clf;
% [xm,ym]=df_measr(1,branch2);
% ym.subfield='l0';
% br_plot(branch2,[],ym,'c');
% ym.subfield='l1';
% br_plot(branch2,[],ym,'b');

% ll=length(branch2.point)
% aa=zeros(2,ll);
% for i=1:ll
%     aa(2,i)=branch2.point(i).parameter(1,3);
%     aa(1,i)=branch2.point(i).parameter(1,4);
% end

% fid=fopen('F:\hopf.txt','wt');
% fprintf(fid,'%12.8f %12.8f\n',aa);
% fclose(fid);

intervals=40;
degree=8;
[psol,stepcond]=p_topsol(first_hopf,1e-4,degree,intervals);
method=df_mthod('psol');
[psol,success]=p_correc(psol,3,stepcond,method.point)

% clear(branch3);
branch3=df_brnch(3,'psol');
branch3.parameter.min_bound(2,:)=[3,0];
branch3.parameter.max_bound(1,:)=[3,20];
branch3.parameter.max_step(1,:)=[3,0.005];
branch3.method.continuation.steplength_growth_factor=1.01
deg_psol=p_topsol(first_hopf,0,degree,intervals);
branch3.point=deg_psol;
branch3.point(2)=psol;
figure(6);clf;
[branch3,s,f,r]=br_contn(branch3,1720)
point=branch3.point(end);
p_ampl=max(point.profile(1,:))-min(point.profile(1,:));
plot(point.parameter(3),p_ampl,'o')

ll=length(branch3.point);

% 提取乘子
% branch31=br_stabl(branch3,0,0)
% max(abs(branch31.point(10).stability.mu))

% 提取极限环的周期
figure(7);clf;
[xm,ym]=df_measr(0,branch3);
ym
ym.field='period';
ym.col=1;

psol=branch3.point(5);
psol.stability=p_stabil(psol,method.stability);
intervals=40;
degree=8;
psol=p_remesh(psol,degree,intervals);
method.point.adapt_mesh_after_correct=1;
method.point.newton_max_iterations=7;
method.point.newton_nmon_iterations=2;
[psol,success]=p_correc(psol,[],[],method.point)

% clear(branch4);
branch4=df_brnch(3,'psol');
branch4.parameter=branch3.parameter;
branch4.point=psol;
psol.parameter(3)=psol.parameter(3)-0.0001;
[psol,success]=p_correc(psol,[],[],method.point,1)
branch4.point(2)=psol;
branch4.method=method;
[xm,ym]=df_measr(0,branch4);
ym.field='period';
ym.col=1;
figure(8); axis auto; hold on;
branch4.method.continuation.plot_measure.x=xm;
branch4.method.continuation.plot_measure.y=ym;
branch4=br_contn(branch4,2000);

% 将周期和乘子保存至同一个文件,首先要将branch4的各点乘子计算出来
% clear(branch41)
branch41=br_stabl(branch4,0,0)

lln=length(branch4.point)

Ac=zeros(6,lln);
  for i=1:lln
    Ac(1,i)=branch4.point(i).parameter(3);
    Ac(2,i)=branch4.point(i).parameter(4);
    Ac(3,i)=max(abs(branch4.point(i).profile(1,:)));
    Ac(4,i)=branch4.point(i).period;
    Ac(5,i)=max(abs(branch41.point(i).stability.mu));
    Ac(6,i)=max(branch4.point(i).profile(1,:))-min(branch4.point(i).profile(1,:))
  end
fprintf(fid,'%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n',Ac);

end

fclose(fid);



fid=fopen('.\forJ2tau05.txt','wt');

index=220;
lll=length(branch4.point(index).profile);
B=zeros(3,lll);
  B(1,:)= branch4.point(index).profile(1,:);
  B(2,:)= branch4.point(index).profile(2,:);
  B(3,:)= 1;
fprintf(fid,'%12.8f %12.8f %12.8f\n',B);
fclose(fid);

rmpath(".\ddebiftool\")