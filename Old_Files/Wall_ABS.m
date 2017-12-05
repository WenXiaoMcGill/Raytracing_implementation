function P = Wall_ABS(P,beta,n,f_concrete,f_wood,f_lime)
if beta(n) == 1;
    P = step(f_concrete,P);
else if beta(n) == 2;
        P = step(f_wood,P);
    else if beta(n) == 3;
            P = step(f_lime,P);
        end
    end
end