function F = myfun(x)
    global As
    F=0;
    Wk = 200:600;
    for i=1:50
        for j =1:length(Wk)
            H=(x(1)+x(2)*i+x(3)*i^2+x(4)*i^3+x(5)*i^4)+(x(6)+x(7)*i+x(8)*i^2+x(9)*i^3+x(10)*i^4)*Wk(j)+...
            (x(11)+x(12)*i+x(13)*i^2+x(14)*i^3+x(15)*i^4)*Wk(j)^2+(x(16)+x(17)*i+x(18)*i^2+x(19)*i^3+x(20)*i^4)*Wk(j)^3+...
            (x(21)+x(22)*i+x(23)*i^2+x(24)*i^3+x(25)*i^4)*Wk(j)^4;
            F=F+(As(i,j)-(H))^2;
        end
    end 
end
        
