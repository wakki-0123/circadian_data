function sleep_soukan(e_all_array_0411)
    % 暫定　被験者の対応はまた詳しく調べる
    test_score = [0.886324022, 0.866727644, 0.911264649, 0.877716281, 0.885269186, 0.855124121, 0.871734582, 0.859731051, 0.900972038, 0.854969546];
    
    % Preallocate the col array for efficiency
    col = zeros(1, 10);
    
    for i = 1:10
        % Compute the correlation between the test score and the corresponding row in e_all_array_0411
        col(i) = corr(test_score(i) * ones(1000, 1), e_all_array_0411(i, 1:1000)');
    end
    
    % Display the correlation values
    disp(col);
end
