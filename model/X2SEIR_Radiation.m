function VARS = X2SEIR_Radiation(X1)

    X = reshape(X1, [30, length(X1) / 30])';
   
    VARS.S = X(1:end, 1:3);
    VARS.E = X(1:end, 4:6);
    VARS.I = X(1:end, 7:9);
    VARS.J = X(1:end, 10:12);
    VARS.V = X(1:end, 13:15);
    VARS.IV = X(1:end, 16:18);
    VARS.VR = X(1:end, 19:21);
    VARS.M = X(1:end, 22:24);
    VARS.R = X(1:end, 25:27);
    VARS.RI = X(1:end, 28:30);

end
