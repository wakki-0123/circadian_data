%function data_cell = data_to_cell(data1,data2,data3,data4,data5,data6,data7,data8,data9,data10,data11,data12,data13,data14,data15)
% サイズが 1x10 の空のセル配列を作成
% data_cell = cell(1, 14);
% data_cell{1} = data1;
% data_cell{2} = data2;
% data_cell{3} = data3;
% data_cell{4} = data4;
% data_cell{5} = data5;
% data_cell{6} = data6;
% data_cell{7} = data7;
% data_cell{8} = data8;
% data_cell{9} = data9;
% data_cell{10} = data10;
% data_cell{11} = data11;
% data_cell{12} = data12;
% data_cell{13} = data13;
% data_cell{14} = data14;
% data_cell{15} = data15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%size 1*2 version
function data_cell = data_to_cell(data1,data2)
data_cell = cell(1, 2);
data_cell{1} = data1;
data_cell{2} = data2;
%data_cell{3} = data3;


% function data_cell = data_to_cell(data1)
% data_cell = cell(1, 1);
% data_cell{1} = data1;

