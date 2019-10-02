xV = [];
xNS = [];
xEW = [];
for i = 1:10000 %30 mins
    ii = i * 3 - 2;
    xV = vertcat(xV, A(ii).d);
    iii = i * 3 - 1;
    xNS = vertcat(xNS, A(iii).d);
    iiii = i * 3;
    xEW = vertcat(xEW, A(iiii).d);
    disp(i)
end


