function sleep_t_FDR_BHanarsis(p_atai, t_atai)
   [p_hc, p_hc_I] = sort(p_atai);

% q値を計算する
q = zeros(length(p_hc), 1);

% 参照値: pi
pi = 0.05;

m = length(p_hc);
%threshold_p_hc = 0;

for i = m:-1:1
    a = (i / m) * pi;
    try
        q(i) = a;
        if p_hc(i) <= q(i)
            %threshold_p_hc = a;
            fprintf(2, 'HC: p値のしきい値: %f, %f; %f\n', a, p_hc(i), t_atai(p_hc_I(i)));
                        break
        end
    catch
        % 例外を処理する（必要に応じて）
        disp("NO")
    end
end

end