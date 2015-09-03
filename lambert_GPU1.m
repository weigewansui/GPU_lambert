function [V1_1, V1_2, V1_3, V2_1, V2_2, V2_3, extremal_distances, exitflag] = lambert_GPU(...
        r1vec_1, r1vec_2, r1vec_3, r2vec_1, r2vec_2, r2vec_3, tf, m, muC)
        
    r1vec = [r1vec_1, r1vec_2, r1vec_3];
    r2vec = [r2vec_1, r2vec_2, r2vec_3];
    
    [V1, V2, extremal_distances, exitflag] = lambert(...
        r1vec, r2vec, tf, m, muC);
    
    V1_1 = V1(1);
    V1_2 = V1(2);
    V1_3 = V1(3);
    V2_1 = V2(1);
    V2_2 = V2(2);
    V2_3 = V2(3);
    
 end