function P = opt2(y)
    [~,length] = size(y);
    half = floor(length/2);
    y_1 = y((half+1):length);
    y_2 = flip(y(1:half),2);
    lower = 19;
    upper = 109;
    Pmat = zeros(1,upper+1);
    for j = 1:(upper-lower)
        k = j+lower;
        Pmat(j) = sum(y_1(1:k:half))+sum(y_2(k:k:half));
        Pmat(j) = Pmat(j)*k;
    end
    maximum = max(Pmat);
    [~,P] = find(Pmat==maximum);
    P = P(1)+lower-1;
end