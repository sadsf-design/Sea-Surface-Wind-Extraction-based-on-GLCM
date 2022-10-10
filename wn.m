function result = wn(m, N)
result = zeros(m, N);    
    for i = 1 : m
        x1 = 1e7*rand;
        y=zeros(1,N);
        k = floor(1e3*rand);
        for j=1:N
           x1=mod(k*x1,10^7);
           y(1,j)=x1/1e7;
        end
        y1=zeros(1,1);n=1;
        while n<=N
            c=randperm(numel(y));b=y(c(1:12));
            b1=b(1:6);b2=b(7:12);
            d=sum(b1)-sum(b2);
            y1(:,n)=d;
            n=n+1;
        end
        result(i, :) = y1;
    end
end