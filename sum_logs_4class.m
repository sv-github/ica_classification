%%% --- nested function to test examples  --- %%%
function [max_class idx_class] = sum_logs_4class(X_vec, Prob_c1,Prob_c2,Prob_c3,Prob_c4, ...
    attrib_prob_c1, attrib_prob_c2, attrib_prob_c3, attrib_prob_c4)

    Sum1_f = log(Prob_c1);    % class 1
    for attr = 1:size(X_vec,2)
        Sum1_f = Sum1_f + log(attrib_prob_c1(X_vec(1,attr),attr));
    end

    Sum2_f = log(Prob_c2);    % class 2
    for attr = 1:size(X_vec,2)
        Sum2_f = Sum2_f + log(attrib_prob_c2(X_vec(1,attr),attr));
    end

    Sum3_f = log(Prob_c3);    % class 3
    for attr = 1:size(X_vec,2)
        Sum3_f = Sum3_f + log(attrib_prob_c3(X_vec(1,attr),attr));
    end

    Sum4_f = log(Prob_c4);    % class 4
    for attr = 1:size(X_vec,2)
        Sum4_f = Sum4_f + log(attrib_prob_c4(X_vec(1,attr),attr));
    end

    [max_class idx_class] = max([Sum1_f Sum2_f Sum3_f Sum4_f]);

end %%% --- nested function for testing examples  --- %%%