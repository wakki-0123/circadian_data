function s = coarsegraining(inputSignal, scaleFactor)
    % 入力シグナルの長さを取得
    signalLength = length(inputSignal);

    % スケールファクターに基づいて新しいサンプル数を計算
    newLength = floor(signalLength / scaleFactor);

    % 粗視化されたシグナルを格納する配列を初期化
    s = zeros(1, newLength);

    % コースグレーニング処理を実行
    for i = 1:newLength
        % 各粗視化されたデータポイントを計算
        startIndex = (i - 1) * scaleFactor + 1;
        endIndex = i * scaleFactor;
        s(i) = mean(inputSignal(startIndex:endIndex));
    end
end
