function [res] = rd(path, num)
    % 文件路径
    fid = fopen(path, 'r');
    if fid == -1
        error('文件无法打开，请检查路径是否正确。');
    end

    % 读取文件的所有行到 cell 数组中
    allLines = textscan(fid, '%s', 'Delimiter', '\n', 'CommentStyle', '%');
    fclose(fid);

    % 从 textscan 返回的 cell 数组中提取行
    % textscan 函数返回一个 cell 数组，其中第一个元素包含所有行
    lines = allLines{1};

    % 获取行数
    numLines = numel(lines);

    % 确保文件至少有10行
    if numLines < num
        error('文件行数少于10行。');
    end

    % 读取最后10行到 cell 数组
    res = lines(end-(num - 1):end);
end