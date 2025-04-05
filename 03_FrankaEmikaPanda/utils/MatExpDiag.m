function result = MatExpDiag(Mat)
    [V,D] = eig(Mat);
    result = V * exp(D) / V;
end