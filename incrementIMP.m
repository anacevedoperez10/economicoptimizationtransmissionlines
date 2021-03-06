function F = myfun(x)
global a b rho n w60 R0 wmax density L Lt kI Tmax dCGdR dTdR dTdW dCGdW dCIdW NT NPH bw Eog
for k=1:n% Derivatives
    aa=(a(1)+a(6)*(x(k+n))+a(11)*(x(k+n))^2+a(16)*(x(k+n))^3+a(21)*(x(k+n))^4);
ab=(a(2)+a(7)*(x(k+n))+a(12)*(x(k+n))^2+a(17)*(x(k+n))^3+a(22)*(x(k+n))^4);
ac=(a(3)+a(8)*(x(k+n))+a(13)*(x(k+n))^2+a(18)*(x(k+n))^3+a(23)*(x(k+n))^4);
ad=(a(4)+a(9)*(x(k+n))+a(14)*(x(k+n))^2+a(19)*(x(k+n))^3+a(24)*(x(k+n))^4);
ae=(a(5)+a(10)*(x(k+n))+a(15)*(x(k+n))^2+a(20)*(x(k+n))^3+a(25)*(x(k+n))^4);
Tx(k)=(density(k))*((a(1)+a(6)*(x(k+n))+a(11)*(x(k+n))^2+a(16)*(x(k+n))^3+a(21)*(x(k+n))^4)+...
                    (a(2)+a(7)*(x(k+n))+a(12)*(x(k+n))^2+a(17)*(x(k+n))^3+a(22)*(x(k+n))^4)*x(k)+...
                    (a(3)+a(8)*(x(k+n))+a(13)*(x(k+n))^2+a(18)*(x(k+n))^3+a(23)*(x(k+n))^4)*x(k)^2+...
                    (a(4)+a(9)*(x(k+n))+a(14)*(x(k+n))^2+a(19)*(x(k+n))^3+a(24)*(x(k+n))^4)*x(k)^3+...
                    (a(5)+a(10)*(x(k+n))+a(15)*(x(k+n))^2+a(20)*(x(k+n))^3+a(25)*(x(k+n))^4)*x(k)^4);     
NL(k)=density(k)*(bw)/10;
Irayo(k)=31*((0.6*NL(k))/Tx(k)-1)^(1/2.6);
Rb(k)=sqrt((rho(k)*Eog*x(k)^2)/(rho(k)*Eog-2*pi*Irayo(k)*x(k)^2));   
CG(k)=NT(k)*b(1)*(rho(k)/Rb(k))^(b(2));
CI(k)=NT(k)*NPH(k)*kI*(x(k+n)-w60);
dTdR(k)=(density(k))*((a(2)+a(7)*(x(k+n))+a(12)*(x(k+n))^2+a(17)*(x(k+n))^3+a(22)*(x(k+n))^4)+...
                     2*(a(3)+a(8)*(x(k+n))+a(13)*(x(k+n))^2+a(18)*(x(k+n))^3+a(23)*(x(k+n))^4)*x(k)+...
                     3*(a(4)+a(9)*(x(k+n))+a(14)*(x(k+n))^2+a(19)*(x(k+n))^3+a(24)*(x(k+n))^4)*x(k)^2+...
                     4*(a(5)+a(10)*(x(k+n))+a(15)*(x(k+n))^2+a(20)*(x(k+n))^3+a(25)*(x(k+n))^4)*x(k)^3);
A1(k)=dTdR(k)*.6*NL(k)/Tx(k)^2;
A2(k)=31*(1/2.6)*(0.6*NL(k)/Tx(k)-1)^(inv(2.6)-1)*A1(k);
Ax(k)=A2(k)*x(k)^2+2*Irayo(k)*x(k);
 A3(k)=(Eog*rho(k)*x(k)*(pi*A2(k)*x(k)^3 + Eog*rho(k)))/((Eog*rho(k) -...
     2*pi*Irayo(k)*x(k)^2)^2*((Eog*rho(k)*x(k)^2)/(Eog*rho(k) - 2*pi*Irayo(k)*x(k)^2))^(1/2));
A3(k)=(1-(2*pi*Irayo(k)*x(k)^2)/(rho(k)*Eog))^(-1/2)+...
    0.5*x(k)*((1-(2*pi*Irayo(k)*x(k)^2)/(rho(k)*Eog))^-(3/2))*((2*pi)/(rho(k)*Eog)*Ax(k));
dCGdRx(k)=-NT(k)*b(2)*b(1)*((rho(k)/Rb(k))^(b(2)-1))*rho(k)*inv(Rb(k)^2)*A3(k); 

dCGdR(k)=-(Eog*NT(k)*b(1)*b(2)*rho(k)*x(k)*(rho(k)/((Eog*rho(k)*x(k)^2)/...
    (Eog*rho(k) - 62*pi*x(k)^2*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x(k) +...
    5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/...
    (5*density(k)*(aa + ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(5/13)))^(1/2))^b(2)*(13*...
    Eog*aa^2*density(k)*rho(k)*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x(k) +...
    5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*...
    (aa + ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(8/13) - 186*pi*NL(k)*ac*x(k)^4 -...
    279*pi*NL(k)*ad*x(k)^5 - 372*pi*NL(k)*ae*x(k)^6 - 93*pi*NL(k)*ab*x(k)^3 +...
    13*Eog*ab^2*density(k)*rho(k)*x(k)^2*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x(k) +...
    5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa +...
    ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(8/13) + 13*Eog*ac^2*density(k)*rho(k)*x(k)^4*(-(5*aa*density(k) -...
    3*NL(k) + 5*ab*density(k)*x(k) + 5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa +...
    ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(8/13) +...
    13*Eog*ad^2*density(k)*rho(k)*x(k)^6*(-(5*aa*density(k) - 3*NL(k) +...
    5*ab*density(k)*x(k) + 5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 +...
    5*ae*density(k)*x(k)^4)/(5*density(k)*(aa + ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 +...
    ab*x(k))))^(8/13) + 13*Eog*ae^2*density(k)*rho(k)*x(k)^8*(-(5*aa*density(k) - 3*NL(k) +...
    5*ab*density(k)*x(k) + 5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 +...
    5*ae*density(k)*x(k)^4)/(5*density(k)*(aa + ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(8/13) +...
    26*Eog*aa*ac*density(k)*rho(k)*x(k)^2*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x(k) +...
    5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa +...
    ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(8/13) + 26*Eog*aa*ad*density(k)*rho(k)*x(k)^3*(-(5*aa*density(k) -...
    3*NL(k) + 5*ab*density(k)*x(k) + 5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa +...
    ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(8/13) + 26*Eog*ab*ac*density(k)*rho(k)*x(k)^3*(-(5*aa*density(k) - 3*NL(k) +...
    5*ab*density(k)*x(k) + 5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa +...
    ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(8/13) + 26*Eog*aa*ae*density(k)*rho(k)*x(k)^4*(-(5*aa*density(k) - 3*NL(k) +...
    5*ab*density(k)*x(k) + 5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa + ac*x(k)^2 + ad*x(k)^3 +...
    ae*x(k)^4 + ab*x(k))))^(8/13) + 26*Eog*ab*ad*density(k)*rho(k)*x(k)^4*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x(k) +...
    5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa + ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 +...
    ab*x(k))))^(8/13) + 26*Eog*ab*ae*density(k)*rho(k)*x(k)^5*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x(k) +...
    5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa + ac*x(k)^2 + ad*x(k)^3 +...
    ae*x(k)^4 + ab*x(k))))^(8/13) + 26*Eog*ac*ad*density(k)*rho(k)*x(k)^5*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x(k) +...
    5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa + ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 +...
    ab*x(k))))^(8/13) + 26*Eog*ac*ae*density(k)*rho(k)*x(k)^6*(-(5*aa*density(k) - 3*NL(k) +...
    5*ab*density(k)*x(k) + 5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa + ac*x(k)^2 +...
    ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(8/13) + 26*Eog*ad*ae*density(k)*rho(k)*x(k)^7*(-(5*aa*density(k) - 3*NL(k) +...
    5*ab*density(k)*x(k) + 5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa +...
    ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(8/13) + 26*Eog*aa*ab*density(k)*rho(k)*x(k)*(-(5*aa*density(k) - 3*NL(k) +...
    5*ab*density(k)*x(k) + 5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa + ac*x(k)^2 +...
    ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(8/13)))/(13*density(k)*(Eog*rho(k) - 62*pi*x(k)^2*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x(k) +...
    5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa + ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 +...
    ab*x(k))))^(5/13))^2*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x(k) + 5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 +...
    5*ae*density(k)*x(k)^4)/(5*density(k)*(aa + ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(8/13)*((Eog*rho(k)*x(k)^2)/(Eog*rho(k) -...
    62*pi*x(k)^2*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x(k) + 5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 +...
    5*ae*density(k)*x(k)^4)/(5*density(k)*(aa + ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(5/13)))*(aa + ac*x(k)^2 +...
    ad*x(k)^3 + ae*x(k)^4 + ab*x(k))^2);

dCIdW(k)=kI*NT(k)*NPH(k);
 dTdW(k)=density(k)*((a(6)+a(7)*x(k)+a(8)*x(k)^2+a(9)*x(k)^3+a(10)*x(k)^4)+...
      2*(a(11)+a(12)*x(k)+a(13)*x(k)^2+a(14)*x(k)^3+a(15)*x(k)^4)*(x(k+n))^1+...
      3*(a(16)+a(17)*x(k)+a(18)*x(k)^2+a(19)*x(k)^3+a(20)*x(k)^4)*(x(k+n))^2+...
      4*(a(21)+a(22)*x(k)+a(23)*x(k)^2+a(24)*x(k)^3+a(25)*x(k)^4)*(x(k+n))^3);
A4(k)=dTdW(k)*.6*NL(k)/Tx(k)^2;
A5(k)=31*(1/2.6)*(0.6*NL(k)/Tx(k)-1)^(inv(2.6)-1)*A4(k);
A6(k)=-.5*x(k)*((1-2*pi*Irayo(k)*x(k)^2/(rho(k)*Eog))^-1.5)*A5(k)*(2*pi*x(k)^2)/(rho(k)*Eog);
dCGdW(k)=-NT(k)*b(2)*b(1)*((rho(k)/Rb(k))^(b(2)-1))*rho(k)*inv(Rb(k)^2)*A6(k);
end
for k=1:n
F(k)=dCGdR(k)+x(6*n+1)*(L(k)/Lt)*dTdR(k)+x(k+4*n)-x(k+5*n);
F(k+n)=dCGdW(k)+dCIdW(k)+x(6*n+1)*(L(k)/Lt)*dTdW(k)-x(k+2*n)+x(k+3*n);
F(k+2*n)=x(k+2*n)*(w60-(x(k+n)));
F(k+3*n)=x(k+3*n)*((x(k+n))-wmax);
F(k+4*n)=x(k+4*n)*(x(k)-R0);
F(k+5*n)=x(k+5*n)*(0-x(k));
end
T=0;
for k=1:n
T=T+Tx(k)*L(k)/Lt;
end
F(6*n+1)=T-Tmax;
% F
% pause