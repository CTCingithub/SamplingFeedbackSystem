function InUnitCircle = IsInUnitCircle(ComplexNum)
    if real(ComplexNum)^2+imag(ComplexNum)^2 <= 1
        InUnitCircle=1;
    else
        InUnitCircle=0;
    end
end

