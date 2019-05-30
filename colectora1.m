clear all, close all, clc
%% Backflashover Rate - Base Case
%% 345 kV Transmission Line - Double Circuit
format long 
ft = 0.3048;
Vnom = 500; %kV
Vo = sqrt(2)*Vnom/sqrt(3);  %kV
% coordinates;
X = [-15 15 -15 -25 -15 15 25 15]*ft; %m
Y = [172 172 1345 104 73 135 104 73]*ft;  %m
R = [0.49 0.49 2.03 2.03 2.03 2.03 2.03 2.03]/100;    % m radius
hs = 7;   %m midspan
Dbundle = 0.5;    %m distance in the bundle
Eo = 1500;    % kV/m 
W = 6; % m Insulatorlength
span = 435; %m 
Po = 760; %mmHg, atmospheric pressure at sea level
H = 0;    %m, Altitude
Ta = 20;  %Celsius degrees, air temperature
C = 0.003671; %Thermal expansion coefficient 1/C
P = Po/10^(H/(18400*(1+C*Ta))); %Presure mmHg
Kt = 0.3855*P/(273+Ta);   % Correction pressure and temperature
T = 70;   %keraunic level tunderstorm days per year
TR = 10;
Rg = 20;    % [Ohm]
W = 2.63; % m Insulatorlength
span = 435; %m 
topcrossarm = [2.7 9.3 15.3];   %[m]
TR = 10;    % [m] Tower Base Radius
%% Methods
Tepri = EpriRed(Vnom,X,Y,topcrossarm,R,Dbundle,span,W,Rg,T);
Tieee = IEEE(Vnom,X,Y,topcrossarm,TR,R,Dbundle,span,W,Rg,T);
Tcigre = cigre_method(Vnom,X,Y,R,span,W,TR,Rg,T);
%% Parameterization Tk
global Rk Wk
Rk = 1:1:50;
Wk = 200:1:600;
BFR = zeros(length(Rk),length(Wk));

for i = Rk
    for j = 1:length(Wk)
        BFR(i,j) = EpriRed(Vnom,X,Y,topcrossarm,R,Dbundle,span,Wk(j)/100,Rk(i),T);
        %BFR(i,j) = IEEE(Vnom,X,Y,topcrossarm,TR,R,Dbundle,span,Wk(j)/100,Rk(i),T);
        %BFR(i,j) = cigre_method(Vnom,X,Y,R,span,Wk(j)/100,TR,Rk(i),T);
    end 
end
%% Least Squares Curve Fitting
global As
As = BFR;
x0 = zeros(1,25);
x0=[0.00127387716756549,0.00165619285997361,0.0501984848515261,-0.000666394404053768,3.27921469183686e-06,0.0267835198865763,-0.00240106761891375,-0.000154216621538726,2.50027968617140e-08,1.48996169830103e-09,-0.000117009056697966,-1.42476859893615e-06,1.03734050242485e-07,4.53242185074383e-09,5.84121171540187e-11,-3.92608112656691e-08,7.10829960549131e-08,-5.95016215426245e-10,-1.56324995533515e-11,2.68569725170256e-13,2.97245501222870e-10,-1.01254540577334e-10,-1.84545110304144e-12,1.64514834646120e-13,-3.05583879847322e-15]
x0=[-1.71149621145269,0.490799068230792,0.0218434491980167,-0.000364492331303949,1.07363893000548e-06,0.0224518619761490,-0.00244450025630907,-0.000150229233577477,3.22061468504681e-06,-5.79298589601012e-10,-0.000154025941103151,-4.17140918253837e-07,6.43482743392851e-08,-7.41094277894265e-10,-1.32683541058916e-10,3.95929047628967e-07,3.74198167521024e-08,1.32426785647594e-10,-4.46306920908083e-12,2.88256453610159e-13,-2.58185594968066e-10,-8.00442535471314e-11,1.08427958984056e-12,-2.51995700981102e-14,1.07994136916595e-16]
x0=[-1.71149621145269,0.490799068230792,0.0218434491980167,-0.000364492331303949,1.07363893000548e-06,0.0224518619761490,-0.00244450025630907,-0.000150229233577477,3.22061468504681e-06,-5.79298589601012e-10,-0.000154025941103151,-4.17140918253837e-07,6.43482743392851e-08,-7.41094277894265e-10,-1.32683541058916e-10,3.95929047628967e-07,3.74198167521024e-08,1.32426785647594e-10,-4.46306920908083e-12,2.88256453610159e-13,-2.58185594968066e-10,-8.00442535471314e-11,1.08427958984056e-12,-2.51995700981102e-14,1.07994136916595e-16]
x0=[-1.45030012019566,0.428067890854809,0.0150884723865904,0.000291950000328748,-1.22233147267906e-06,0.00243374062044858,-0.00135984191881789,-0.000162355976372466,-2.56909445892791e-06,4.66245528215599e-10,-3.27220724452856e-08,1.52357864545020e-08,4.02652906127660e-07,1.35189722703924e-08,3.00000892718842e-11,2.78483675988418e-08,2.89551269657220e-09,-1.14415792424063e-10,-2.88590075510244e-11,-1.89546892365479e-13,3.05853911499560e-13,-1.14480691252984e-11,2.96549313861846e-13,2.08151002854526e-15,4.52562844419792e-16]
x0=[-1.45528955910841,0.411583484978611,0.00673244224134971,0.000630132906844289,-2.45619313096732e-06,0.00870868239382704,-0.00241868914604116,-1.23209021466055e-05,-6.74931419244360e-06,-5.34050067942074e-10,-1.95808098353684e-07,2.22224602140542e-08,2.09134783296954e-08,2.38759828650788e-08,1.51812186941811e-10,-2.33827753742979e-08,1.33985470878761e-08,1.12230674272334e-10,-3.90086156848372e-11,-5.61151418636535e-13,-3.20059498983315e-11,-9.92470931458885e-12,-5.51604208848438e-13,3.09303022430254e-14,5.45999984808206e-16]
x0=[-1.28175558864981,0.486944667197670,0.0162796195535683,0.000532963745035843,-4.58743527407805e-06,-0.00257244763122654,-0.00217378944306810,-0.000225848465843316,-2.77631650965398e-06,9.85519110811828e-11,6.04840431163072e-05,8.39093620561482e-09,1.06015895012050e-06,1.44126295697381e-09,1.99574887878166e-10,-1.23489709617257e-07,1.92354277520006e-09,-1.52571279467304e-09,1.31698204400833e-12,-6.53329005755942e-13,1.34459015048611e-11,1.15617361851027e-11,2.76311950873532e-14,1.20022942353116e-14,5.48956923599540e-16]
x0=[-1.13615366783882,0.386536822859161,0.0204284173766757,0.000433074018244550,-3.38570488638837e-06,0.00190855168535316,-0.00186594524895810,-0.000219071330569411,-3.10846122687952e-06,-3.64513418500652e-11,2.59294707158351e-05,-9.69834519257203e-08,9.38249478537298e-07,6.13988974115746e-09,1.64834420727788e-10,-6.09620966523212e-08,6.43517787742177e-09,-1.50179193194743e-09,-5.41130437702587e-12,-5.83272713956387e-13,-1.09552568580519e-12,-3.25043414930232e-13,4.85882150016661e-13,7.93810986933548e-15,5.51410878327862e-16]
x0=[2.57231675736281,0.566763470620598,-0.0211415847464837,0.00132754997012515,-6.90564669516609e-06,-0.0282302258911667,-0.00277623219147546,9.32454124366149e-05,-9.84225329363900e-06,3.24813321488654e-08,7.98515313315166e-05,2.36707283117884e-06,1.10851811524053e-07,2.56278640133477e-08,-1.50449952480972e-11,-2.63260328319338e-08,1.21452089513314e-09,-5.56347734612423e-10,-3.71980763694403e-11,-8.53561721539972e-16,-8.44254839789147e-11,1.30368006073539e-12,2.56519748323408e-13,2.98751430521554e-14,-9.32699670171218e-17]
x0=[2.53465889877699,-0.0421984346409335,0.0329135788486635,-0.000121131228230497,2.22308658013595e-06,-0.0293594552236457,0.00144959563800640,-0.000269667791049930,6.10458538029075e-07,-2.74577733014864e-08,0.000113860544485835,-6.62543354113905e-06,6.94451527373884e-07,4.71612658577582e-09,7.34932025513316e-11,-1.67732974490062e-07,8.04308954701348e-09,-3.87888740412993e-10,-2.75358302159655e-11,3.09796613446275e-15,7.20436077724228e-11,-1.98023641080176e-13,-3.90645666907296e-13,3.39622709903678e-14,-1.21708482531817e-16]
x0=[2.80909216300600,0.191541267236160,0.0197937271787954,-6.79887179821098e-07,5.80666682735882e-07,-0.0289747440900898,-0.000865355083709777,-0.000159035447515910,-8.67118855549825e-08,-7.15385358740062e-09,0.000124670661468869,-3.06300552260099e-06,6.94389418324282e-07,-2.92518098160765e-09,4.91604309985369e-11,-1.84913358522983e-07,6.79703456120779e-09,-9.11084211173999e-10,-1.02505906506708e-12,-1.23204977852173e-14,6.40510007566537e-11,1.30949301916804e-12,8.56998633768115e-14,1.24259017267877e-14,-1.06995894873053e-16]
x0=[2.86132215864575,0.0940309051555633,0.0338392132818588,6.33621112475851e-07,-4.77823458377046e-06,-0.0300456681899995,-0.000587441861446656,-0.000270817519462979,-3.16825622659198e-07,4.52468621592318e-08,0.000125003549153690,-1.53597680526391e-07,9.10338663273354e-07,5.29126772202983e-10,-1.49040145787628e-10,-2.27479698855147e-07,4.23036329320093e-09,-1.41621984417427e-09,6.02343790812110e-13,2.14420660358133e-13,1.49746039793975e-10,-4.63315463061784e-12,8.21803054081679e-13,-8.10308299680270e-16,-1.18426037546083e-16]
x0=[2.83722063369395,0.120015970492283,0.0336925587413547,-5.69555056020303e-06,-4.84302283904610e-06,-0.0309818985710153,-0.000725633280848452,-0.000270095625333877,-4.21368399896559e-07,4.86884467783004e-08,0.000126682876854041,4.80794328476922e-07,8.99160536738524e-07,1.38057634757084e-09,-1.67844049556811e-10,-2.24575249826300e-07,3.08018484864899e-09,-1.41965422425315e-09,-2.88721405292986e-14,2.37685813831418e-13,1.44491007128855e-10,-4.09703727715999e-12,8.64529562235261e-13,-1.87387527165836e-15,-1.17479352359148e-16]
x0=[2.76535120168782,0.106433913248196,0.0336764786011796,3.67061176362950e-06,-4.87654718694623e-06,-0.0301445034962532,-0.000593020176645579,-0.000271286700913909,-4.62192296448163e-07,4.85234067367623e-08,0.000124351934807619,-4.12070073097699e-08,9.08307647871324e-07,1.38315681916561e-09,-1.66261298487578e-10,-2.23368311169653e-07,3.99104030456489e-09,-1.43747673290984e-09,1.43369618451030e-14,2.35790223210480e-13,1.45700882999004e-10,-4.67934767461175e-12,8.74991214438357e-13,-1.85610704801436e-15,-1.17486574333425e-16]
x0=[2.75864337646668,0.034473380242091,0.0391386176339167,-0.000135079724660226,-3.68208033382066e-06,-0.0274796655915574,-0.000304579335465087,-0.000299406839002771,2.66967337875967e-07,4.21744725622123e-08,0.000107791385380334,3.90616246573419e-08,9.39431473513265e-07,5.25767619004378e-10,-1.58395413837482e-10,-1.92402195548960e-07,3.11390001845305e-09,-1.43447516555475e-09,-3.48301700454578e-14,2.35461996758345e-13,1.28502799776365e-10,-4.33898849263444e-12,8.78003392865210e-13,-1.84211459547052e-15,-1.17445365958245e-16]
x0=[2.75920441950744,0.0408449302445427,0.0387687599054210,-0.000123962449357726,-3.80591473772760e-06,-0.0271644009119494,-0.000378556517931010,-0.000296974182370828,2.10059898427073e-07,4.28556594932559e-08,0.000105363503768229,3.55960249825564e-07,9.34284348359103e-07,5.68068902212560e-10,-1.59238003264415e-10,-1.85830368468068e-07,2.43495029564606e-09,-1.42606754215499e-09,7.42219720375152e-15,2.35517732542581e-13,1.22736768250301e-10,-3.78657556131978e-12,8.69626872183993e-13,-1.84743886029577e-15,-1.17440107426252e-16]
x0=[2.75631319982440,0.0302003465929522,0.0373644165928017,-6.27719807178712e-05,-4.32681214350192e-06,-0.0288456470543954,3.77239740185839e-05,-0.000298910663470697,-1.96881803810503e-08,4.44643006505737e-08,0.000113827677706462,-1.72885978630840e-06,9.70446110528340e-07,6.59282576364944e-10,-1.57139635512471e-10,-1.97076248914955e-07,5.72795808571968e-09,-1.48200297687518e-09,-1.11887506138640e-14,2.28126269674276e-13,1.25165425315408e-10,-5.23746450026062e-12,8.77186920351755e-13,-1.36550392174854e-15,-1.15257921235044e-16]
% Least squares curve fitting
tol1 = 1e-5;
tol2 = 1e-5;
tol3 = 1e-5;
LB = [];
UB = [];
A = [];
B = [];
Aeq = [];
Beq = [];
options=optimset('Algorithm','active-set','Display','iter','TolFun',tol1,'TolCon',tol2,'TolX',tol3,'MaxIter',300,'MaxFunEvals',10000);
%FMINCON Call objetive function and non-linear constraints 
[x,fval,exitflag,output,lambda,grad,hessian]=fmincon('minq4',x0,A,B,Aeq,Beq,LB,UB,'restr',options);
%% Incremental costs of transmission lightning protection systems
% Four cases, n=1, n=3, n=12 and n=300 sections
global a b kI rho Tmax density CG CI L n Lt Tx w60  wmax R0 dCGdR dCGdW dCIdW dTdR dTdW NPH NT bw Eog Irayo Rb NL
% System Data
% T parameters (Table 1) BFOR against density, string length and impulse
% resistance Red Book Method
% a=(1/3.6)*[2.75631319982440 0.0302003465929522 0.0373644165928017 -6.27719807178712e-05...
%    -4.32681214350192e-06 -0.0288456470543954 3.77239740185839e-05 -0.000298910663470697...
%    -1.96881803810503e-08  4.44643006505737e-08 0.000113827677706462 -1.72885978630840e-06...
%    9.70446110528340e-07 6.59282576364944e-10 -1.57139635512471e-10 -1.97076248914955e-07...
%    5.72795808571968e-09 -1.48200297687518e-09 -1.11887506138640e-14 2.28126269674276e-13...
%   1.25165425315408e-10 -5.23746450026062e-12 8.77186920351755e-13 -1.36550392174854e-15 -1.15257921235044e-16];

a = x/3.6;
% Ground cost parameters Figure 5 
b =[55.52532416 1.221891523];
kI =25;  % incremental tower insulation cost $/cm
bw = 201.5938;    %m shadow
Eog = 500;    % kV/m ionization
Lt = 100; % km line lenght
Tmax = 1.075; % BF tripout rate 
w60 = 263;    % 60Hz and switching insulation
wmax = 445;   % [cm] max string length
R0 = 50000000;    % [ohms] default tower resistance
%% Case study definition
%% Case 1 One section
% n = 1;%number of sections
% rho = [1833];
% density = [3.6];

%% Case 2: 3 sections
% n = 3;
% rho = [1000 1500 3000];
% density = [4.6 3.6 2.6];

%% Case 3: 12 towers
n = 12;
rho=[850 750 1250 1150 975 1125 1500 2400 4000 4200 2000 1800];
density=[4.6 4.6 4.6 4.6 3.6 3.6 3.6 3.6 2.6 2.6 2.6 2.6];
%% Case 4: 300 towers
% n=300;
% rhox = [850 750 1250 1150 975 1125 1500 2400 4000 4200 2000 1800];
% densityx = [4.6 4.6 4.6 4.6 3.6 3.6 3.6 3.6 2.6 2.6 2.6 2.6 ];
% j=0;
% e=.10;
% e2=0.02;
% 
% for i=1:8  
%     for k=1:25
%         j=j+1;
%         rho(j)=(1-e/2)*rhox(i)+(k-1)*(e)*rhox(i)/24;      %0.01*rhox(i)*(-1 + (2)*rand(1,1));
%         density(j)=densityx(i);%+0.01*densityx(i)*(-1 + (2)*rand(1,1));
%     end
% end
% for i=9:12  
%     for k=1:25
%         j=j+1;
%         rho(j)=(1.0+e2/2)*rhox(i)-(k-1)*(e2)*rhox(i)/24;      %0.01*rhox(i)*(-1 + (2)*rand(1,1));
%         density(j)=densityx(i);%+0.01*densityx(i)*(-1 + (2)*rand(1,1));
%     end
% end
% %close all
% figure
% plot(rho)
% e3=0;
% j=0;
% for i=1:12  
%     for k=1:25
%         j=j+1;
%         density(j)=(1-e3/2)*densityx(i)+(k-1)*(e3)*densityx(i)/24;      %0.01*rhox(i)*(-1 + (2)*rand(1,1));
%         %density(j)=densityx(i);%+0.01*densityx(i)*(-1 + (2)*rand(1,1));
%     end
% end

NPH = ones(1,n)*6;% Number of phases per section
NT = ones(1,n)*300/n; %Number of towers per section
for k=1:n
    L(k)=Lt/n;%section length km
end																								
  
%% Setup the optimization problem
x0a = 19.9626*ones(1,n);
x0b = 263.0000*ones(1,n);
tol1 = 1e-6;  
tol2 = 1e-6;
tol3 = 1e-6;
LB = [];
UB = [];
A = [];
B = [];
Aeq = [];
Beq = [];
x0 = horzcat(x0a,x0b);
options=optimset('Algorithm','interior-point','Display','iter','TolFun',tol1,'TolCon',tol2,'TolX',tol3,'MaxIter',30000000,'MaxFunEvals',10000000);
%FMINCON Call objetive function and non-linear constraints 
[x1,fval,exitflag,output,lambda,grad,hessian]=fmincon('minCostimp',x0,A,B,Aeq,Beq,LB,UB,'restrCost2imp',options);

%Print Results
for k=1:n
    Tx(k)=(density(k))*((a(1)+a(6)*(x1(k+n))+a(11)*(x1(k+n))^2+a(16)*(x1(k+n))^3+a(21)*(x1(k+n))^4)+...
                    (a(2)+a(7)*(x1(k+n))+a(12)*(x1(k+n))^2+a(17)*(x1(k+n))^3+a(22)*(x1(k+n))^4)*x1(k)+...
                    (a(3)+a(8)*(x1(k+n))+a(13)*(x1(k+n))^2+a(18)*(x1(k+n))^3+a(23)*(x1(k+n))^4)*x1(k)^2+...
                    (a(4)+a(9)*(x1(k+n))+a(14)*(x1(k+n))^2+a(19)*(x1(k+n))^3+a(24)*(x1(k+n))^4)*x1(k)^3+...
                    (a(5)+a(10)*(x1(k+n))+a(15)*(x1(k+n))^2+a(20)*(x1(k+n))^3+a(25)*(x1(k+n))^4)*x1(k)^4);    
    Rimpulse(k)=x1(k);
    wstringlength(k)=x1(k+n);
    NL(k)=density(k)*(bw)/10;
    Irayo(k)=31*((0.6*NL(k))/Tx(k)-1)^(1/2.6);
    Rb(k)=sqrt(rho(k)*Eog*x1(k)^2/(rho(k)*Eog-2*pi*Irayo(k)*x1(k)^2));
    Ig(k)=rho(k)*Eog/(2*pi*Rb(k)^2);
    CG(k)=NT(k)*b(1)*(rho(k)/Rb(k))^(b(2));
    CI(k)=NT(k)*NPH(k)*kI*(x1(k+n)-w60);
    count(k)=k;
end
Ctot = sum(CG)+sum(CI);
plot(count,CG+CI,'b-',count,CG,'g-',count,CI,'c--')
hold on
yyaxis right
plot(count,Tx,'r-',count,Tmax,'k-')

%% otlgimpulse vaying Tmax from 0.1 to 1.075
% b=[	55.20367932 1.22425597];
% b=[	55.52532427 1.221891523];

b = [55.52532416 1.221891523];
kI = 25;    %incremental tower insulation cost $/cm
bw = 201.5938;  % [m]
Eog = 500;  % kV/m
Lt = 100; % [km]    Line Length
for u=1:110
    Tmax=.1+u/100;1.075;%BF tripout rate 
    w60=263;%60Hz and switching insulation
    wmax=445;%cm
    R0=50000000;%ohms default tower resistance
%% % Case 1 One section
%n=1;
%  rho=[1833];
%  density=[3.6];

% 
% %Case 2 3 sections
%n=3;%number of sections
% rho=[1000 1500 3000];
% density=[4.6 3.6 2.6];

%% Case 3 12 towers
% n=12;%number of sections
% rho=[	850	750	1250	1150  975        1125        1500        2400	 4000	4200	2000	1800	];
% density=[	4.6	4.6	4.6	4.6		3.6	3.6	3.6	3.6	2.6	2.6	2.6	2.6 ];
%% Case 4 300 towers
%  n=300;%number of sections
%  rhox=[	850	750	1250	1150  975        1125        1500        2400	 4000	4200	2000	1800	];
%  densityx=[	4.6	4.6	4.6	4.6		3.6	3.6	3.6	3.6	2.6	2.6	2.6	2.6 ];
% j=0;
% e=.10;
% e2=0.02;
% for i=1:8  
% for k=1:25
%     j=j+1;
%     rho(j)=(1-e/2)*rhox(i)+(k-1)*(e)*rhox(i)/24;      %0.01*rhox(i)*(-1 + (2)*rand(1,1));
%     density(j)=densityx(i);%+0.01*densityx(i)*(-1 + (2)*rand(1,1));
% end
% end
% for i=9:12  
% for k=1:25
%     j=j+1;
%     rho(j)=(1.0+e2/2)*rhox(i)-(k-1)*(e2)*rhox(i)/24;      %0.01*rhox(i)*(-1 + (2)*rand(1,1));
%     density(j)=densityx(i);%+0.01*densityx(i)*(-1 + (2)*rand(1,1));
% end
% end
% %close all
% %plot(rho)
% e3=0;
% j=0;
% for i=1:12  
%     for k=1:25
%         j=j+1;
%         density(j)=(1-e3/2)*densityx(i)+(k-1)*(e3)*densityx(i)/24;      %0.01*rhox(i)*(-1 + (2)*rand(1,1));
%         %density(j)=densityx(i);%+0.01*densityx(i)*(-1 + (2)*rand(1,1));
%     end
% end

NPH=ones(1,n)*6;% Number of phases per section
NT=ones(1,n)*300/n; %Number of towers per section

for k=1:n
    L(k)=Lt/n;%section length km
end																								
    
x0a = 19.9626*ones(1,n);
x0b = 263.0000*ones(1,n);
tol1 = 1e-6;  
tol2 = 1e-6;
tol3 = 1e-6;
LB=[];
UB=[];
A=[];
B=[];
Aeq=[];
Beq=[];
x0=horzcat(x0a,x0b);

%for u=1:1;  
%Tmax=u/10;
%options=optimset('Display','iter','LargeScale','on','ActiveConstrTol',1);
options=optimset('Algorithm','interior-point','Display','iter','TolFun',tol1,'TolCon',tol2,'TolX',tol3,'MaxIter',30000000,'MaxFunEvals',10000000);
%FMINCON Call objetive function and non-linear constraints 
[x2,fval,exitflag,output,lambda,grad,hessian]=fmincon('minCostimp',x0,A,B,Aeq,Beq,LB,UB,'restrCost2imp',options);
%Cost structure
    for k=1:n
        Tx(k)=(density(k))*((a(1)+a(6)*(x2(k+n))+a(11)*(x2(k+n))^2+a(16)*(x2(k+n))^3+a(21)*(x2(k+n))^4)+...
                            (a(2)+a(7)*(x2(k+n))+a(12)*(x2(k+n))^2+a(17)*(x2(k+n))^3+a(22)*(x2(k+n))^4)*x2(k)+...
                            (a(3)+a(8)*(x2(k+n))+a(13)*(x2(k+n))^2+a(18)*(x2(k+n))^3+a(23)*(x2(k+n))^4)*x2(k)^2+...
                            (a(4)+a(9)*(x2(k+n))+a(14)*(x2(k+n))^2+a(19)*(x2(k+n))^3+a(24)*(x2(k+n))^4)*x2(k)^3+...
                            (a(5)+a(10)*(x2(k+n))+a(15)*(x2(k+n))^2+a(20)*(x2(k+n))^3+a(25)*(x2(k+n))^4)*x2(k)^4);     
        NL(k)=density(k)*(bw)/10;
        Irayo(k)=31*((0.6*NL(k))/Tx(k)-1)^(1/2.6);
        Rb(k)=sqrt(rho(k)*Eog*x2(k)^2/(rho(k)*Eog-2*pi*Irayo(k)*x2(k)^2));
        Ig(k)=rho(k)*Eog/(2*pi*Rb(k)^2);
        CG(k)=NT(k)*b(1)*(rho(k)/Rb(k))^(b(2));
        CI(k)=NT(k)*NPH(k)*kI*(x2(k+n)-w60);
        count(k)=k;
    end

%% Determine incremental costs
xf=x2;
mu=lambda.ineqnonlin';
lmd=lambda.eqnonlin;
x0=horzcat(xf,mu,lmd(1));
x2=x0;
    
for k=1:n% Derivatives
aa = (a(1)+a(6)*(x2(k+n))+a(11)*(x2(k+n))^2+a(16)*(x2(k+n))^3+a(21)*(x2(k+n))^4);
ab = (a(2)+a(7)*(x2(k+n))+a(12)*(x2(k+n))^2+a(17)*(x2(k+n))^3+a(22)*(x2(k+n))^4);
ac = (a(3)+a(8)*(x2(k+n))+a(13)*(x2(k+n))^2+a(18)*(x2(k+n))^3+a(23)*(x2(k+n))^4);
ad = (a(4)+a(9)*(x2(k+n))+a(14)*(x2(k+n))^2+a(19)*(x2(k+n))^3+a(24)*(x2(k+n))^4);
ae = (a(5)+a(10)*(x2(k+n))+a(15)*(x2(k+n))^2+a(20)*(x2(k+n))^3+a(25)*(x2(k+n))^4);
Tx(k) = (density(k))*((a(1)+a(6)*(x2(k+n))+a(11)*(x2(k+n))^2+a(16)*(x2(k+n))^3+a(21)*(x2(k+n))^4)+...
                    (a(2)+a(7)*(x2(k+n))+a(12)*(x2(k+n))^2+a(17)*(x2(k+n))^3+a(22)*(x2(k+n))^4)*x2(k)+...
                    (a(3)+a(8)*(x2(k+n))+a(13)*(x2(k+n))^2+a(18)*(x2(k+n))^3+a(23)*(x2(k+n))^4)*x2(k)^2+...
                    (a(4)+a(9)*(x2(k+n))+a(14)*(x2(k+n))^2+a(19)*(x2(k+n))^3+a(24)*(x2(k+n))^4)*x2(k)^3+...
                    (a(5)+a(10)*(x2(k+n))+a(15)*(x2(k+n))^2+a(20)*(x2(k+n))^3+a(25)*(x2(k+n))^4)*x2(k)^4);     
NL(k) = density(k)*(bw)/10;
Irayo(k)=31*((0.6*NL(k))/Tx(k)-1)^(1/2.6);
Rb(k)=sqrt((rho(k)*Eog*x2(k)^2)/(rho(k)*Eog-2*pi*Irayo(k)*x2(k)^2));   

CG(k)=NT(k)*b(1)*(rho(k)/Rb(k))^(b(2));
CI(k)=NT(k)*NPH(k)*kI*(x2(k+n)-w60);

dTdR(k)=(density(k))*((a(2)+a(7)*(x2(k+n))+a(12)*(x2(k+n))^2+a(17)*(x2(k+n))^3+a(22)*(x2(k+n))^4)+...
                     2*(a(3)+a(8)*(x2(k+n))+a(13)*(x2(k+n))^2+a(18)*(x2(k+n))^3+a(23)*(x2(k+n))^4)*x2(k)+...
                     3*(a(4)+a(9)*(x2(k+n))+a(14)*(x2(k+n))^2+a(19)*(x2(k+n))^3+a(24)*(x2(k+n))^4)*x2(k)^2+...
                     4*(a(5)+a(10)*(x2(k+n))+a(15)*(x2(k+n))^2+a(20)*(x2(k+n))^3+a(25)*(x2(k+n))^4)*x2(k)^3);

dCGdR(k)=-(Eog*NT(k)*b(1)*b(2)*rho(k)*x2(k)*(rho(k)/((Eog*rho(k)*x2(k)^2)/...
    (Eog*rho(k) - 62*pi*x2(k)^2*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x2(k) +...
    5*ac*density(k)*x2(k)^2 + 5*ad*density(k)*x2(k)^3 + 5*ae*density(k)*x2(k)^4)/...
    (5*density(k)*(aa + ac*x2(k)^2 + ad*x2(k)^3 + ae*x2(k)^4 + ab*x2(k))))^(5/13)))^(1/2))^b(2)*(13*...
    Eog*aa^2*density(k)*rho(k)*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x2(k) +...
    5*ac*density(k)*x2(k)^2 + 5*ad*density(k)*x2(k)^3 + 5*ae*density(k)*x2(k)^4)/(5*density(k)*...
    (aa + ac*x2(k)^2 + ad*x2(k)^3 + ae*x2(k)^4 + ab*x2(k))))^(8/13) - 186*pi*NL(k)*ac*x2(k)^4 -...
    279*pi*NL(k)*ad*x2(k)^5 - 372*pi*NL(k)*ae*x2(k)^6 - 93*pi*NL(k)*ab*x2(k)^3 +...
    13*Eog*ab^2*density(k)*rho(k)*x2(k)^2*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x2(k) +...
    5*ac*density(k)*x2(k)^2 + 5*ad*density(k)*x2(k)^3 + 5*ae*density(k)*x2(k)^4)/(5*density(k)*(aa +...
    ac*x2(k)^2 + ad*x2(k)^3 + ae*x2(k)^4 + ab*x2(k))))^(8/13) + 13*Eog*ac^2*density(k)*rho(k)*x2(k)^4*(-(5*aa*density(k) -...
    3*NL(k) + 5*ab*density(k)*x2(k) + 5*ac*density(k)*x2(k)^2 + 5*ad*density(k)*x2(k)^3 + 5*ae*density(k)*x2(k)^4)/(5*density(k)*(aa +...
    ac*x2(k)^2 + ad*x2(k)^3 + ae*x2(k)^4 + ab*x2(k))))^(8/13) +...
    13*Eog*ad^2*density(k)*rho(k)*x2(k)^6*(-(5*aa*density(k) - 3*NL(k) +...
    5*ab*density(k)*x2(k) + 5*ac*density(k)*x2(k)^2 + 5*ad*density(k)*x2(k)^3 +...
    5*ae*density(k)*x2(k)^4)/(5*density(k)*(aa + ac*x2(k)^2 + ad*x2(k)^3 + ae*x2(k)^4 +...
    ab*x2(k))))^(8/13) + 13*Eog*ae^2*density(k)*rho(k)*x2(k)^8*(-(5*aa*density(k) - 3*NL(k) +...
    5*ab*density(k)*x2(k) + 5*ac*density(k)*x2(k)^2 + 5*ad*density(k)*x2(k)^3 +...
    5*ae*density(k)*x2(k)^4)/(5*density(k)*(aa + ac*x2(k)^2 + ad*x2(k)^3 + ae*x2(k)^4 + ab*x2(k))))^(8/13) +...
    26*Eog*aa*ac*density(k)*rho(k)*x2(k)^2*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x2(k) +...
    5*ac*density(k)*x2(k)^2 + 5*ad*density(k)*x2(k)^3 + 5*ae*density(k)*x2(k)^4)/(5*density(k)*(aa +...
    ac*x2(k)^2 + ad*x2(k)^3 + ae*x2(k)^4 + ab*x2(k))))^(8/13) + 26*Eog*aa*ad*density(k)*rho(k)*x2(k)^3*(-(5*aa*density(k) -...
    3*NL(k) + 5*ab*density(k)*x2(k) + 5*ac*density(k)*x2(k)^2 + 5*ad*density(k)*x2(k)^3 + 5*ae*density(k)*x2(k)^4)/(5*density(k)*(aa +...
    ac*x2(k)^2 + ad*x2(k)^3 + ae*x2(k)^4 + ab*x2(k))))^(8/13) + 26*Eog*ab*ac*density(k)*rho(k)*x2(k)^3*(-(5*aa*density(k) - 3*NL(k) +...
    5*ab*density(k)*x2(k) + 5*ac*density(k)*x2(k)^2 + 5*ad*density(k)*x2(k)^3 + 5*ae*density(k)*x2(k)^4)/(5*density(k)*(aa +...
    ac*x2(k)^2 + ad*x2(k)^3 + ae*x2(k)^4 + ab*x2(k))))^(8/13) + 26*Eog*aa*ae*density(k)*rho(k)*x2(k)^4*(-(5*aa*density(k) - 3*NL(k) +...
    5*ab*density(k)*x2(k) + 5*ac*density(k)*x2(k)^2 + 5*ad*density(k)*x2(k)^3 + 5*ae*density(k)*x2(k)^4)/(5*density(k)*(aa + ac*x2(k)^2 + ad*x2(k)^3 +...
    ae*x2(k)^4 + ab*x2(k))))^(8/13) + 26*Eog*ab*ad*density(k)*rho(k)*x2(k)^4*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x2(k) +...
    5*ac*density(k)*x2(k)^2 + 5*ad*density(k)*x2(k)^3 + 5*ae*density(k)*x2(k)^4)/(5*density(k)*(aa + ac*x2(k)^2 + ad*x2(k)^3 + ae*x2(k)^4 +...
    ab*x2(k))))^(8/13) + 26*Eog*ab*ae*density(k)*rho(k)*x2(k)^5*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x2(k) +...
    5*ac*density(k)*x2(k)^2 + 5*ad*density(k)*x2(k)^3 + 5*ae*density(k)*x2(k)^4)/(5*density(k)*(aa + ac*x2(k)^2 + ad*x2(k)^3 +...
    ae*x2(k)^4 + ab*x2(k))))^(8/13) + 26*Eog*ac*ad*density(k)*rho(k)*x2(k)^5*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x2(k) +...
    5*ac*density(k)*x2(k)^2 + 5*ad*density(k)*x2(k)^3 + 5*ae*density(k)*x2(k)^4)/(5*density(k)*(aa + ac*x2(k)^2 + ad*x2(k)^3 + ae*x2(k)^4 +...
    ab*x2(k))))^(8/13) + 26*Eog*ac*ae*density(k)*rho(k)*x2(k)^6*(-(5*aa*density(k) - 3*NL(k) +...
    5*ab*density(k)*x2(k) + 5*ac*density(k)*x2(k)^2 + 5*ad*density(k)*x2(k)^3 + 5*ae*density(k)*x2(k)^4)/(5*density(k)*(aa + ac*x2(k)^2 +...
    ad*x2(k)^3 + ae*x2(k)^4 + ab*x2(k))))^(8/13) + 26*Eog*ad*ae*density(k)*rho(k)*x2(k)^7*(-(5*aa*density(k) - 3*NL(k) +...
    5*ab*density(k)*x2(k) + 5*ac*density(k)*x2(k)^2 + 5*ad*density(k)*x2(k)^3 + 5*ae*density(k)*x2(k)^4)/(5*density(k)*(aa +...
    ac*x2(k)^2 + ad*x2(k)^3 + ae*x2(k)^4 + ab*x2(k))))^(8/13) + 26*Eog*aa*ab*density(k)*rho(k)*x2(k)*(-(5*aa*density(k) - 3*NL(k) +...
    5*ab*density(k)*x2(k) + 5*ac*density(k)*x2(k)^2 + 5*ad*density(k)*x2(k)^3 + 5*ae*density(k)*x2(k)^4)/(5*density(k)*(aa + ac*x2(k)^2 +...
    ad*x2(k)^3 + ae*x2(k)^4 + ab*x2(k))))^(8/13)))/(13*density(k)*(Eog*rho(k) - 62*pi*x2(k)^2*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x2(k) +...
    5*ac*density(k)*x2(k)^2 + 5*ad*density(k)*x2(k)^3 + 5*ae*density(k)*x2(k)^4)/(5*density(k)*(aa + ac*x2(k)^2 + ad*x2(k)^3 + ae*x2(k)^4 +...
    ab*x2(k))))^(5/13))^2*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x2(k) + 5*ac*density(k)*x2(k)^2 + 5*ad*density(k)*x2(k)^3 +...
    5*ae*density(k)*x2(k)^4)/(5*density(k)*(aa + ac*x2(k)^2 + ad*x2(k)^3 + ae*x2(k)^4 + ab*x2(k))))^(8/13)*((Eog*rho(k)*x2(k)^2)/(Eog*rho(k) -...
    62*pi*x2(k)^2*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x2(k) + 5*ac*density(k)*x2(k)^2 + 5*ad*density(k)*x2(k)^3 +...
    5*ae*density(k)*x2(k)^4)/(5*density(k)*(aa + ac*x2(k)^2 + ad*x2(k)^3 + ae*x2(k)^4 + ab*x2(k))))^(5/13)))*(aa + ac*x2(k)^2 +...
    ad*x2(k)^3 + ae*x2(k)^4 + ab*x2(k))^2);

dCIdW(k)=kI*NT(k)*NPH(k);
dTdW(k)=density(k)*((a(6)+a(7)*x2(k)+a(8)*x2(k)^2+a(9)*x2(k)^3+a(10)*x2(k)^4)+...
      2*(a(11)+a(12)*x2(k)+a(13)*x2(k)^2+a(14)*x2(k)^3+a(15)*x2(k)^4)*(x2(k+n))^1+...
      3*(a(16)+a(17)*x2(k)+a(18)*x2(k)^2+a(19)*x2(k)^3+a(20)*x2(k)^4)*(x2(k+n))^2+...
      4*(a(21)+a(22)*x2(k)+a(23)*x2(k)^2+a(24)*x2(k)^3+a(25)*x2(k)^4)*(x2(k+n))^3);
  
A4(k) = dTdW(k)*.6*NL(k)/Tx(k)^2;
A5(k) = 31*(1/2.6)*(0.6*NL(k)/Tx(k)-1)^(inv(2.6)-1)*A4(k);
A6(k)=-.5*x2(k)*((1-2*pi*Irayo(k)*x2(k)^2/(rho(k)*Eog))^-1.5)*A5(k)*(2*pi*x2(k)^2)/(rho(k)*Eog);
dCGdW(k)=-NT(k)*b(2)*b(1)*((rho(k)/Rb(k))^(b(2)-1))*rho(k)*inv(Rb(k)^2)*A6(k);
end
    Ctot(u)=sum(CG)+sum(CI);
    Cgtot(u)=+sum(CG);
    syslambda(u)=-lmd;
    countl(u)=u;
end
%%
figure
plot(Ctot)
figure
plot(syslambda)
%% Smooting curves 
% CT = smooth(Ctot,0.1,'rloess');
Tk = linspace(0.1,1.2,110);
% figure
% plot(Tk',CT/10^6,'','LineWidth',1)
% xlabel('T_k')
% ylabel('C_{GK}+C_{Ik} [M$]')
%%
xlabel('T_k')
ylabel('\partial C_{Gk}/\partial T_k   [M$]')
%% 
f = polyfit(1./(Tk),smooth(smooth(smooth(Ctot)))',1);
f1 = polyval(f,1./Tk);
figure, hold on
plot(Tk,f1/10^6,'','LineWidth',1)
plot(Tk,Ctot/10^6,'','LineWidth',1)
xlabel('T_k')
ylabel('\partial C_{Gk}/\partial T_k   [M$]')
%%
% lambda = polyfit(-1./(Tk.^2),syslambda,1);
% lambda1 = polyval(lambda,-1./(Tk.^2));
% figure, hold on 
% plot(Tk,lambda1)
% plot(Tk,syslambda)
%%
lbda = diff(f1/10^6);
figure, hold on 
plot(linspace(0.1,1.2,109),lbda*10)