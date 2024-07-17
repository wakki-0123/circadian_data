function [s, r] = IAAFT(data, c, maxiter)
    % IAAFT: Iteratively Adjusted Amplitude and Frequency Transform
    
    % データの長さを取得
    data_length = length(data);
    
    % サロゲートデータ配列をGPU上で初期化
    s = gpuArray.zeros(data_length, c);
    
    % c個のサロゲートデータセットを生成
    for j = 1:c
        % データをランダムにシャッフル
        r = data(randperm(data_length));
        
        % オリジナルデータの振幅スペクトルを計算
        A = abs(fft(data));  % fft関数を使用してgpuArrayを引数にする
        
        for iter = 1:maxiter
            % サロゲートデータをフーリエ変換
            S = fft(r);
            
            % 振幅をオリジナルデータのものに置き換え
            S = A .* exp(1i * angle(S));
            
            % 逆フーリエ変換して時間領域に戻す
            r = real(ifft(S));
            
            % 順位付け
            [~, I] = sort(r);
            [~, J] = sort(data);
            r(I) = data(J);
        end
        
        % サロゲートデータを保存
        s(:, j) = r;
    end
end
