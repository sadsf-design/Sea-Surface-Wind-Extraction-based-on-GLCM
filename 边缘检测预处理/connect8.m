function [number, stack, label] = connect8( mIn )
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[m, n] = size(mIn);%��ȡ����ͼƬ����������
label = ones(m, n);%labels����������������ǩ�Ĳ���
stack = zeros(2, 10000);%��ջ���飬�����������Ҫ���
top = 0;%��ջָ��
count = 0;%������
amount = zeros(10000, 1);
TuCou = 0; %�����С������
for i = 2 : m - 2
    for j = 2 : n - 2
        TuCou = 0;
        if (mIn(i, j) == 0) && (label(i,j) == 1)
            count = count + 1;
            top = top + 1;
            stack(1, top) = i;
            stack(2, top) = j;
            label(i, j) = 0;
            TuCou = TuCou + 1;
        end
        in = i;
        jn = j;
        while top ~= 0
            in = stack(1, top);
            jn = stack(2, top);
            if in <= 0 || in - 1 <= 0
                break;
            end
            if jn <= 0 || jn - 1 <= 0
                break;
            end
            if jn > n || jn + 1 > n
                break;
            end
            if in > m || in + 1 > m
                break;
            end
            top = top - 1;
            if (mIn(in - 1, jn) == 0) && (label(in - 1, jn) == 1)
                top = top + 1;
                stack(1, top) = in - 1;
                stack(2, top) = jn;
                label(in - 1, jn) = 0;
                TuCou = TuCou + 1;
            end
            if (mIn(in + 1, jn) == 0) && (label(in + 1, jn) == 1)
                top = top + 1;
                stack(1, top) = in + 1;
                stack(2, top) = jn;
                label(in + 1, jn) = 0;
                TuCou = TuCou + 1;
            end
            if (mIn(in, jn - 1) == 0) && (label(in, jn - 1) == 1)
                top = top + 1;
                stack(1, top) = in;
                stack(2, top) = jn - 1;
                label(in, jn - 1) = 0;
                TuCou = TuCou + 1;
            end
            if (mIn(in - 1, jn - 1) == 1) && (label(in - 1, jn - 1) == 1)
                top = top + 1;
                stack(1, top) = in - 1;
                stack(2, top) = jn - 1;
                label(in - 1, jn - 1) = 0;
                TuCou = TuCou + 1;
            end
            if (mIn(in + 1, jn - 1) == 0) && (label(in + 1, jn - 1) == 1)
                top = top + 1;
                stack(1, top) = in + 1;
                stack(2, top) = jn - 1;
                label(in + 1, jn - 1) = 0;
                TuCou = TuCou + 1;
            end
            if (mIn(in, jn + 1) == 0) && (label(in, jn + 1) == 1)
                top = top + 1;
                stack(1, top) = in;
                stack(2, top) = jn + 1;
                label(in, jn + 1) = 0;
                TuCou = TuCou + 1;
            end
            if (mIn(in - 1, jn + 1) == 0) && (label(in - 1, jn + 1) == 1)
                top = top + 1;
                stack(1, top) = in - 1;
                stack(2, top) = jn + 1;
                label(in - 1, jn + 1) = 0;
                TuCou = TuCou + 1;
            end
            if (mIn(in + 1, jn + 1) == 0) && (label(in + 1, jn + 1) == 1)
                top = top + 1;
                stack(1, top) = in + 1;
                stack(2, top) = jn + 1;
                label(in + 1, jn + 1) = 0;
                TuCou = TuCou + 1;
            end
        end
        if count ~= 0 && TuCou > amount(count)
            amount(count) = TuCou;
        end
    end
end
number = count;


end