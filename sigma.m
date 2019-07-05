%sigma->neural activity check fn
function y=sigma(x)
    if x<0
        y=0;
    elif x>0 & x<1
        y=x;
    else
        y=1;
    end
end
