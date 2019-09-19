function P = opt(y)
    [~,length] = size(y);
    half = floor(length/2);
    y_1 = y((half+1):length);
    y_2 = flip(y(1:half),2);
    Pmat = zeros(1,71);
    for j = 1:71
        k = j+39;
        Pmat(j) = sum(y_1(1:k:half))+sum(y_2(k:k:half));
        Pmat(j) = Pmat(j)*k;
    end
    maximum = max(Pmat);
    [~,P] = find(Pmat==maximum);
    P = P(1)+39;
end