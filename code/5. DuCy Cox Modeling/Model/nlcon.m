function [c,ceq] = nlcon(x)

    global SS_W SS_B;
    global TR_B TR_W;


    ceq = SS_W - TR_W;
    c = 0.85 - 1 + SS_B;


end
