mv=size(data,1);
for i =1:size(data,1)
    W1 = constructW_PKN(data{i}',5, 1);  %   X : d*n  / 5 neighbors /  is symmetric
    D = diag(sum(W1));
    L =D-W1;
    [F{i}, ~, ~]=eig1(L,length(unique(gt)),0); 
end
n=size(data{1},1);
c=length(unique(gt));
SS=zeros(n,n);
for i=1:mv
    SS=SS+eye(n)-F{i}*F{i}';
end
[Fstar,temp,ev]=eig1(SS,c,0);
Y=eye(n,c);
Y(:,1) = 1;

para1=[0.001,0.01,0.1,1,10];
para2=[0.001,0.01,0.1,1,10];
para3=[1,100,1000];

best_acc = 0;
best_params = [];

for i=1:length(para1)
    for j=1:length(para2)
        for k = 1:length(para3)             
              [result] = MOMSC(data, gt,F,Y,Fstar,para1(i), para2(j), para3(k));
               acc = result(1);
               nmi = result(2);
               ARI = result(3);
               f = result(4);
               e=result(5);
               if acc > best_acc 
                     best_acc = acc;
                     best_nmi = nmi;
                     best_ARI = ARI;
                     best_f = f;
                     best_e = e;
               end
        end
    end
end
 disp(best_acc);
 disp(best_nmi);
 disp(best_ARI);
 disp(best_f);
  disp(best_e);