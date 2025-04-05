function [real_part, imag_part] = complex_seperate(complex_polynomial, re_id, im_id, highest_order)
    % complex_polynomial should be converted into the form containing re_id
    % and im_id, such as a + b i being transformed to re_id*a+ im_id*b.
    
    % Noting that i=i, i^2=-1, i^3=-i, i^4=1. Therefore, the value of i^j
    % can be justified with mod(j,4).
    
    % For real_part, substitute re_id with 1, im_id^j with 1 if mod(j,4)=0,
    % -1 if mod(j,4)=2, 0 if mod(j,4)=1 or 3.
    
    % For imag_part, substitute re_id with 0, im_id^j with 0 if mod(j,4)=0
    % or 2, 1 if mod(j,4)=1, -1 if mod(j,4)=3.
    real_part = subs(complex_polynomial, re_id, 1);
    imag_part = subs(complex_polynomial, re_id, 0);
    
    for j = highest_order : -1 : 1
        if mod(j,4) == 0
            real_part = subs(real_part, im_id^j, 1);
            imag_part = subs(imag_part, im_id^j, 0);
        elseif mod(j,4) == 1
            real_part = subs(real_part, im_id^j, 0);
            imag_part = subs(imag_part, im_id^j, 1);
        elseif mod(j,4) == 2
            real_part = subs(real_part, im_id^j, -1);
            imag_part = subs(imag_part, im_id^j, 0);
        elseif mod(j,4) == 3
            real_part = subs(real_part, im_id^j, 0);
            imag_part = subs(imag_part, im_id^j, -1);
        end
    end
end

