function sleep_t_FDR_BHanarsis_0823(p_atai, t_atai)
   [p_hc, p_hc_I] = sort(p_atai);

   % q値を計算する
   q = zeros(length(p_hc), 1);

   % 参照値: pi
   pi = 0.05;

   m = length(p_hc);

   for i = m:-1:1
       a = (i / m) * pi;
       q(i) = a;
       fprintf('i: %d, a: %f, q(i): %f, p_hc(i): %f\n', i, a, q(i), p_hc(i));
       
       if p_hc(i) <= q(i)
           fprintf(2, 'HC: p値のしきい値: %f, %f; %f\n', a, p_hc(i), t_atai(p_hc_I(i)));
           break;
       end
   end

end
