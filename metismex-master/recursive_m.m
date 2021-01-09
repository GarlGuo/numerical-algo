% 多级方法，递归处理
function M = recursive_m(M, level, isLeafPred)
    isLeaf = isLeafPred(level);
    nparts = 2;
    [~, n] = size(M);
    [i,j,k]=find(tril(M,-1));
    [p,q1,q2]=mexmetis(i,j,k,[n],[nparts]);
    
    if isLeaf
        P = 1:size(M, 1); 
        P(1:end) = amd(M);
        R = chol(M(P, P));
        tmp = R * R';
        M = tmp \ M;
    else
        m1=q1(2)-1;
        M=M(p,p);
        %分治
        M_l = recursive_m(M(1:m1, 1:m1), level-1, isLeafPred);
        M_r = recursive_m(M(m1+1:end, m1+1:end), level-1, isLeafPred);
        M(1:m1) = M_l; 
        M(m1+1:end) = M_r;
        
        W=-M(q2(1):q1(2)-1,q2(2):q1(3)-1);
        [m, ~]=size(W);
        X1=zeros(m,1);
         
     for i=1:m
         X1(i)=sqrt(norm(W(i,:)));
         W(i,:)=W(i,:)/X1(i);
     end
         X1=spdiags(X1,[0],m,m);
         X2=W;
         E=[zeros(q2(1)-q1(1),m);X1;zeros(q2(2)-q1(2),m);X2'];
         B=M+E*E';
         
        B11=B(1:m1,1:m1);B22=B(m1+1:end,m1+1:end);
        p2=1:n;
        P=1:size(B11,1);P(1,1:q2(1)-q1(1))=amd(B11(1:q2(1)-1,1:q2(1)-1));R1=chol(B11(P,P));p2(1,1:m1)=P;
        P=1:size(B22,1);P(1,1:q2(2)-q1(2))=amd(B22(1:q2(2)-q1(2),1:q2(2)-q1(2)));R2=chol(B22(P,P));p2(1,m1+1:end)=P+m1;
        M=M(p2,p2);
        
        % 低秩修正，rank=2
        k=2;
        [U,S,V]=svds(@(b,tflag)Afun(b,tflag,E,R1,R2,m1),size(E),k);
        U=U*S;
        V=V(:,1:k);
        H=inv(eye(k)-U'*E*V);
        
        M = MLRsolveNonLeaf(M,R1,R2,m1,U,H);
    end
    
    
end

function y=Afun(b,tflag,E,R1,R2,m1)
    if strcmp(tflag,'notransp')
        x=E*b;
        x1=x(1:m1,:);y1=R1\(R1'\x1);
        x2=x(m1+1:end,:);y2=R2\(R2'\x2);
        y=[y1;y2];
    else
        x1=b(1:m1,:);y1=R1\(R1'\x1);
        x2=b(m1+1:end,:);y2=R2\(R2'\x2);
        y=E'*[y1;y2];
    end
end

% 2.17 method. handle leaf


function y=MLRsolveNonLeaf(x,R1,R2,m1,U,H)
    y=zeros(size(x,1),1);
    y(1:m1,:)=R1\(R1'\x(1:m1,:));
    y(m1+1:end,:)=R2\(R2'\x(m1+1:end,:));
    y0=U*(H*(U'*x));
    y=y+y0;
end
