function [makespan]=flowshop(y,t)
[R x]= sort(y,2,'descend');
[n,m]=size(t);
C(x(1),1)=t(x(1),1);
for i=2:n %baris
    C(x(i),1)=C(x(i-1),1)+t(x(i),1); 
end
for j=2:m %kolom
    C(x(1),j)=C(x(1),j-1)+t(x(1),j); 
end
for i=2:n
    for j=2:m
        C(x(i),j)=max(C(x(i-1),j),C(x(i),j-1))+t(x(i),j);
        
    end
end
s=C(x(n),m);
makespan=s;
end

