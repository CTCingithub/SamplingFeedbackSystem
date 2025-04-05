function result = MatExpJordan(Mat, trunc_rank)
    [P, J] = jordan(Mat);
    result = P * MatExpTruc(J, trunc_rank) / P;
end

