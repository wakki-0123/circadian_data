function [p_values,t_values] = sleep_scatter(e_all_array_0411) %upper
    % 暫定 被験者の対応はまた詳しく調べる(どうやら合っていそう)
    test_score = [0.886324022, 0.866727644, 0.911264649, 0.877716281, 0.885269186, 0.855124121, 0.871734582, 0.859731051, 0.900972038, 0.854969546];
    entropy = e_all_array_0411;
    num = 10;
    scale = 200;
    
    % E_person を初期化
    E_person = zeros(num, scale);

    % 各スケールのエントロピー値を E_person に格納
    for i = 1:scale
        E_person(:, i) = entropy(:, i); % 10×1 が1000個
    end

    % 相関係数の計算と表示
    correlation_coefficients = zeros(1, scale);
    for i = 1:scale
        [correlation_coefficients(i),p_values(i)] = corr(E_person(:, i), test_score','Type','Spearman');
        n= num; % サンプル数
        t_values(i) = correlation_coefficients(i) * sqrt((n - 2) / (1 - correlation_coefficients(i)^2));

    end

    % 相関係数を表示
    %disp('相関係数:');
    %disp(correlation_coefficients);
    disp(p_values)
    %disp(t_values)
    
    % 図の作成
    figure;
    plot(1:scale, correlation_coefficients, 'b-', 'LineWidth', 2);
    xlabel('Time Scale');
    ylabel('Spearman Correlation Coefficient');
    title('Correlation Coefficients over Different Time Scales');
    grid on;







% 
%     % 相関係数の絶対値が大きい上位10個を取り出す
%     [sorted_coeffs, sorted_indices] = sort(abs(correlation_coefficients), 'descend');
%     top_10_indices = sorted_indices(1:upper); %元は10(右辺)
%     top_10_coeffs = correlation_coefficients(top_10_indices);
% 
%     % 結果を表示
%     disp("上位"+ upper+"個の相関係数:");
%     for i = 1:upper
%         disp(['Index: ', num2str(top_10_indices(i)), '　Time Scale[s]: ',num2str(top_10_indices(i)*5) , '　Correlation Coefficient: ', num2str(top_10_coeffs(i))]);
%     end
% % 上位10個の相関係数に対応する散布図を作成
%     figure;
%     for i = 1:upper
%         subplot(upper/2, 2, i); % 2列5行のサブプロットを作成
%         scatter(E_person(:, top_10_indices(i)), test_score, 'filled');
%         title(['Index: ', num2str(top_10_indices(i)), '　Time Scale[s]: ',num2str(top_10_indices(i)*5),', Corr: ', num2str(top_10_coeffs(i))]);
%         xlabel('Fuzzy Entropy');
%         ylabel('Sleep Efficiency');
%         grid on;
%         % 各点に番号を付ける
%         for j = 1:num
%             text(E_person(j, top_10_indices(i)), test_score(j), num2str(j), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
%         end
%     end
% 
%     % グラフを表示
%     hold off;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function sleep_scatter(e_all_array_0411)
% % 暫定　被験者の対応はまた詳しく調べる
%     test_score = [0.886324022, 0.866727644, 0.911264649, 0.877716281, 0.885269186, 0.855124121, 0.871734582, 0.859731051, 0.900972038, 0.854969546];
%     entropy = e_all_array_0411;
%     num = 10;
%     scale = 1000;
%     % E_person を初期化
%     E_person = zeros(num, scale);
% 
%     for i=1:scale
%         E_person(:,i) = entropy(:,i); %10×1 が1000個
%     end
% 
% 
% 
% % 相関係数の計算と表示
%     correlation_coefficients = zeros(1, scale);
%     for i = 1:scale
%         correlation_coefficients(i) = corr(E_person(:, i), test_score');
%     end
% 
%     % 相関係数を表示
%     disp('相関係数:');
%     disp(correlation_coefficients);
