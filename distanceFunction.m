function d = distfun(X0,X)
    d = zeros(1, size(X,1));
    for iter = 1:size(X,1)
        sum = 0.0;
        i = 1;
        while(i < size(X,2))
            sum_train_vi = degtorad(abs(X0(i))) + degtorad(abs(X0(i+1)));
            sum_test_vi = degtorad(abs(X(iter,i))) + degtorad(abs(X(iter, i+1)));
            sum = sum + exp(abs(sum_test_vi - sum_train_vi));
            i = i + 2;
        end    
        d(iter)=sum;
    end
    d = d';
end