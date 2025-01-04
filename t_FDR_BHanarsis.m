function t_FDR_BHanarsis(p_atai, t_atai)

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
        %disp('No')
    end
end

end





% function t_FDR_BHanarsis(p_atai)
% 
% 
% 
% [p_hc, p_hc_I] = sort(p_atai);
% 
% % calculate q-values
% q=zeros(length(p_hc),1);
% 
% t_hc=t_atai(p_hc_I);
% 
% %reference value:pi
% pi=0.05;
% 
% %length of p value:
% m = length(p_hc);
% threshold_p_hc=0;
% for i=m:-1:0
%     a=(i/m)*pi;
%     %fprintf(2,'%d %f %f\n',i,a,p_ad(i));
%     try
%         q(i)=a;
%         if p_hc(i)<=q(i)
%             threshold_p_hc=a;
%             fprintf(2,'HC: Threshold of p value: %f, %f; %f\n',a,p_hc(i),t_hc(i));
%             break
%         end
%     catch
%     end
% end
% 
% %%%%%%%
% 
% 
