function d = expDiff(v1,v2)
    sum = 0.0;
    i = 1;
    while(i < size(v1,2))
        sum_train_vi = degtorad(abs(v1(i))) + degtorad(abs(v1(i+1)));
        sum_test_vi = degtorad(abs(v2(i))) + degtorad(abs(v2(i+1)));
        sum = sum + exp(abs(sum_test_vi - sum_train_vi));
        i = i + 2;
    end    
    d=sum;
end