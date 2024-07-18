test_a = [-0.416147 -0.909297         0;
    0.909297 -0.416147         0;
    0         0         1]
logm_test_a = logm(test_a)

test_b = [-0.999999 -0.00159255           0;
    0.00159255   -0.999999           0;
    0           0           1]
logm_test_b = logm(test_b)

test_c = [ 1 -0.381773   1.30117;
    1.30117         1 -0.381773;
    -0.381773   1.30117         1]
logm_test_c = logm(test_c)

yaw_d = 3.14;
ea_d = [0;0;yaw_d];
expm(hatMap(ea_d))
vR_d = reshape(expm(hatMap(ea_d)),[9,1])