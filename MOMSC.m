function [result]=MOMSC(X,s,F,Y,Fstar,alpha,beta,lambda)
% s is the true class label.
mv=size(X,1);
n=size(X{1},1);
c=length(unique(s));
wv=ones(mv,1)/mv;
C=zeros(n);
B=repmat(C,[1,1,mv]);
R = eye(c);
for ii=1:100
    Ysold = Y;
    CC=zeros(n);
    for i=1:mv
        L1=diag(B(:,:,i)*ones(n,1))-B(:,:,i);
        [U,temp,ev]=eig1(L1,c,0);
        W{i}=U*U';
        W{i}(find(W{i}<0))=0;
        for ij=1:n
            d=distance(F{i},n,ij);
            xx=X{i}*X{i}';
            B1=B(:,:,i);
            Z(:,ij)=(xx'+lambda*eye(n))\(xx(:,ij)+lambda*B1(:,ij) - 1/4*d'); 
            Z(find(Z<0))=0;
            Z= (Z+Z')/2;
        end
        Dv = diag(sum(Z));
        Dv=inv(sqrtm(Dv));
        lv = eye(n)-Dv*Z*Dv;
        Lmx= wv(i)*(Fstar*Fstar');
        A=lv-Lmx;
        F{i}=eig1(A,c,0);
        CC=CC+wv(i)*F{i}*F{i}';
        L2=diag(W{i})*ones(n,1)'-W{i};
        Q=Z-alpha/(2*lambda)*L2;
        Q=Q-diag(diag(Q));
        AA=(Q+Q)/2;
        AA(find(AA<0))=0;
        B(:,:,i)=AA;
        wv(i)=1/norm(F{i}*F{i}'-Fstar*Fstar','fro');
    end
    for rep = 1:100
        M = 4* CC * Fstar+beta*Y*R';
        [Um,~,Vm] = svd(M,'econ');
        Fstar = Um*Vm';  
        fobj(rep+1) = 2*trace( Fstar' * CC * Fstar ) + beta * trace(R'* Fstar' * Y );
        if rep>4 && ((fobj(rep)-fobj(rep-1))/fobj(rep)<1e-3)
            break;
        end
    end
  
    [Ur,~,Vr] = svd(Fstar'*Y,'econ');
    R = Ur*Vr';
    Y = zeros(n,c);
    Fstar1=Fstar*R;
    for i = 1:n
        [~, max_index] = max(Fstar1(i, :));
        Y(i, max_index) = 1;
    end
    if  ii>4 && (norm(Y - Ysold)/norm(Ysold)<1e-5)
        break
    end
end

Y_column_vector = [];
for i = 1:size(Y, 1)
    nonzero_indices = find(Y(i, :) ~= 0);
    Y_column_vector = [Y_column_vector; nonzero_indices(:)];
end
result=Measures(s,Y_column_vector);
end
function [all]=distance(F,n,ij)
  for ji=1:n
            all(ji)=(norm(F(ij,:)-F(ji,:)))^2;
  end
end   
       
