                powers_1 = y_tmp^1;
                powers_2 = y_tmp^2;
                powers_3 = y_tmp^3;
                powers_4 = y_tmp^4;
                powers_5 = y_tmp^5;
                powers_6 = y_tmp^6;
                powers_7 = y_tmp^7;
                powers_8 = y_tmp^8;
                powers_9 = y_tmp^9;
                powers_10 = y_tmp^10;
                powers_11 = y_tmp^11;
                powers_12 = y_tmp^12;
                powers_13 = y_tmp^13;
                powers_14 = y_tmp^14;
                powers_15 = y_tmp^15;
                powers_16 = y_tmp^16;
                powers_17 = y_tmp^17;
                powers_18 = y_tmp^18;
                powers_19 = y_tmp^19;
                powers_20 = y_tmp^20;
                powers_21 = y_tmp^21;
                powers_22 = y_tmp^22;
                powers_23 = y_tmp^23;
                powers_24 = y_tmp^24;
                powers_25 = y_tmp^25;


                 % sigma itself
                sig = 4/3 ...
            + powers_1*an_1 ...
            + powers_2*an_2 ...
            + powers_3*an_3 ...
            + powers_4*an_4 ...
            + powers_5*an_5 ...
            + powers_6*an_6 ...
            + powers_7*an_7 ...
            + powers_8*an_8 ...
            + powers_9*an_9 ...
            + powers_10*an_10 ...
            + powers_11*an_11 ...
            + powers_12*an_12 ...
            + powers_13*an_13 ...
            + powers_14*an_14 ...
            + powers_15*an_15 ...
            + powers_16*an_16 ...
            + powers_17*an_17 ...
            + powers_18*an_18 ...
            + powers_19*an_19 ...
            + powers_20*an_20 ...
            + powers_21*an_21 ...
            + powers_22*an_22 ...
            + powers_23*an_23 ...
            + powers_24*an_24 ...
            + powers_25*an_25; 
    
    % dsigma / dx (derivative)
%     dsigdx = ( (1:25).*[1, powers(1:24)] ) * an;
            dsigdx = 1 * 1 * an_1 ...
             + 2 * powers_1 * an_2 ...
             + 3 * powers_2 * an_3 ...
             + 4 * powers_3 * an_4 ...
             + 5 * powers_4 * an_5 ...
             + 6 * powers_5 * an_6 ... 
             + 7 * powers_6 * an_7 ...
             + 8 * powers_7 * an_8 ...
             + 9 * powers_8 * an_9 ...
             + 10 * powers_9 * an_10 ...
             + 11 * powers_10 * an_11 ...
             + 12 * powers_11 * an_12 ...
             + 13 * powers_12 * an_13 ...
             + 14 * powers_13 * an_14 ...
             + 15 * powers_14 * an_15 ...
             + 16 * powers_15 * an_16 ...
             + 17 * powers_16 * an_17 ...
             + 18 * powers_17 * an_18 ...
             + 19 * powers_18 * an_19 ...
             + 20 * powers_19 * an_20 ...
             + 21 * powers_20 * an_21 ...
             + 22 * powers_21 * an_22 ...
             + 23 * powers_22 * an_23 ...
             + 24 * powers_23 * an_24 ...
             + 25 * powers_24 * an_25; 
    
    % d2sigma / dx2 (second derivative)
%     d2sigdx2 = ( (1:25).*(0:24).*[1/y, 1, powers(1:23)] ) * an;
    
   d2sigdx2  = 1 * 0 * 1/y_tmp * an_1 ...
             + 2 * 1 * 1 * an_2 ...
             + 3 * 2 * powers_1 * an_3 ...
             + 4 * 3 * powers_2 * an_4 ...
             + 5 * 4 * powers_3 * an_5 ...
             + 6 * 5 * powers_4 * an_6 ...
             + 7 * 6 * powers_5 * an_7 ...
             + 8 * 7 * powers_6 * an_8 ...
             + 9 * 8 * powers_7 * an_9 ...
             + 10 * 9 * powers_8 * an_10 ...
             + 11 * 10 * powers_9 * an_11 ...
             + 12 * 11 * powers_10 * an_12 ...
             + 13 * 12 * powers_11 * an_13 ...
             + 14 * 13 * powers_12 * an_14 ...
             + 15 * 14 * powers_13 * an_15 ...
             + 16 * 15 * powers_14 * an_16 ...
             + 17 * 16 * powers_15 * an_17 ...
             + 18 * 17 * powers_16 * an_18 ...
             + 19 * 18 * powers_17 * an_19 ...
             + 20 * 19 * powers_18 * an_20 ...
             + 21 * 20 * powers_19 * an_21 ...
             + 22 * 21 * powers_20 * an_22 ...
             + 23 * 22 * powers_21 * an_23 ...
             + 24 * 23 * powers_22 * an_24 ...
             + 25 * 24 * powers_23 * an_25;
         
    % d3sigma / dx3 (third derivative)
%     d3sigdx3 = ( (1:25).*(0:24).*(-1:23).*[1/y/y, 1/y, 1, powers(1:22)] ) * an;
    
    d3sigdx3 = 1 * 0 * -1 * 1/y_tmp/y_tmp * an_1 ...
             + 2 * 1 * 0 *  1/y_tmp * an_2 ...
             + 3 * 2 * 1 * 1 * an_3 ...
             + 4 * 3 * 2 * powers_1 * an_4 ...
             + 5 * 4 * 3 * powers_2 * an_5 ...
             + 6 * 5 * 4 * powers_3 * an_6 ...
             + 7 * 6 * 5 * powers_4 * an_7 ...
             + 8 * 7 * 6 * powers_5 * an_8 ...
             + 9 * 8 * 7 * powers_6 * an_9 ...
             + 10 * 9 * 8 * powers_7 * an_10 ...
             + 11 * 10 * 9 * powers_8 * an_11 ...
             + 12 * 11 * 10 * powers_9 * an_12 ...
             + 13 * 12 * 11 * powers_10 * an_13 ...
             + 14 * 13 * 12 * powers_11 * an_14 ...
             + 15 * 14 * 13 * powers_12 * an_15 ...
             + 16 * 15 * 14 * powers_13 * an_16 ...
             + 17 * 16 * 15 * powers_14 * an_17 ...
             + 18 * 17 * 16 * powers_15 * an_18 ...
             + 19 * 18 * 17 * powers_16 * an_19 ...
             + 20 * 19 * 18 * powers_17 * an_20 ...
             + 21 * 20 * 19 * powers_18 * an_21 ...
             + 22 * 21 * 20 * powers_19 * an_22 ...
             + 23 * 22 * 21 * powers_20 * an_23 ...
             + 24 * 23 * 22 * powers_21 * an_24 ...
             + 25 * 24 * 23 * powers_22 * an_25;