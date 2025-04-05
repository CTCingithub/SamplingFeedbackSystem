function result = MatExpTrunc(Mat, rank_trunc)
    result=zeros(length(Mat));
    if rank_trunc == -1
        result = expm(Mat);
    else
        for i=0:rank_trunc
            result = result + mpower(Mat, i) / factorial(i);
        end
    end
end

