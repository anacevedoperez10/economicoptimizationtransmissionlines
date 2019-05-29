function [BFR] = cigre_method(Vnom,X,Y,R,span,W,TR,Rg,T)
    %% CIGRE method
    %% Setup Data
    %% Coordinates 
    % [m]  Conductor Horizontal Distance
    Xg = X(1:2);    % Shield Wire
    X = X(3:8);     % [m] A B C C' B' A'
    HT = Y(1);  % Tower Height
    Y = Y(3:8);  % [m] A B C C' B' A'
    Eo = 1500;  % kV/m
    C = 300;    %light speed
    NE = 2; % Shield Wire
    %% Input Data
    % TR Tower Base Radius
    % Rg Footing Resistance
    % T Keraunic Level
    % W Insulator Length
    % span [m]
    Rg0 = 0;
    tf = 2;
    tf0 = 0;
    nc = 6; % Conductors 
    b = Xg(2)-Xg(1); % Horizontal separation between two shield wires
    %% Tower Modeling 
    % Span travel time
    TS = span/(0.9*C);

    % Tower-top voltage
    TV = 1.8*820*W;

    % Shield wire corona adjustment
    RC = fsolve(@(RC) RC*log((2*HT)/RC)-TV/Eo,[1],optimoptions('fsolve','Display','off')); 

    % Self surge impedance of each shield wire
    GZ = 60*sqrt(log(2*HT/R(1))*log(2*HT/RC));

    % Combined surge impedance
    if NE == 1
        GC = GZ;
    else 
        BM = abs(Xg(2)-Xg(1));
        AM = sqrt((2*HT)^2+BM^2);
        GM = 60*log(AM/BM);
        GC = (GZ+GM)/2;
        Zg = GC;
    end

    % Mutual impedances between conductors and shield wires
    for i = 1 : NE
        for j = 1 : nc
            AM = sqrt((HT+Y(j))^2+(Xg(i)-X(j))^2);
            BM = sqrt((HT-Y(j))^2+(Xg(i)-X(j))^2);
            CZ(j, i) = 60 * log(AM / BM);
        end
    end

    % Coupling factors
    for i = 1:nc
        if NE == 1
            CF(i) = CZ(i, 1)/GC;
        else   
            CF(i) = (CZ(i, 1) + CZ(i, 2))/ 2/ GC;
        end
    end

    % Tower surge response 
    TU = HT/(0.85*C);
    ZT = 30*log(2*(HT^2+TR*TR/4)/(TR*TR/4));
    %% Main Loop
    while abs(tf-tf0)>0.1
        while abs(Rg-Rg0)>0.1
            alphaR = Zg/(Zg+2*Rg);
            alphaT = (ZT-Rg)/(ZT+Rg);
            Re = Zg*Rg/(Zg+2*Rg);

            Kspan = 1-alphaR*(1 - alphaT) * ...
                ((1 - 2*TS/tf) * double(2*TS/tf < 1) + alphaR*alphaT*(1 - 4*TS/tf) * double(4*TS/tf < 1) ...
                 + (alphaR*alphaT)^2*(1 - 6*TS/tf) * double(6*TS/tf < 1) + (alphaR*alphaT)^3*(1 - 8*TS/tf) ...
                 * double(8*TS/tf < 1));
            tao = TS*Zg/Rg;

            for i = 1:nc
                TA(i) = (HT-Y(i))/(0.85*C);
                KTA(i) = Re+alphaT*ZT*TA(i)/tf;
            end

            KTT = Re+alphaT*ZT*TU/tf;
            VTT = Kspan*KTT;
            VTA = Kspan*KTA;
            VF = Kspan*Re;
            CFO = 1800;
            Vo = Vnom*sqrt(2)/sqrt(3);

            for i = 1:nc
                VI(i) = VTA(i) - CF(i)*VTT;
                VIF(i) = (1 - CF(i)) * VF;
                deltaV(i) = VI(i) - VIF(i);
                CFO_NS(i) = CFO * (0.977 + 2.82/tao) * (1 + deltaV(i)/VIF(i)) *...
                (1 - 0.2*(1 + deltaV(i)/VIF(i))*Vo/CFO) * (1 - 0.09*(1 + 10/tao)*deltaV(i)/VIF(i)) *...
                exp(-deltaV(i)/VIF(i) * tf/13);
            end

            k = find(CF == min(CF));
            IC1 = (CFO_NS(k) - 0.4*sqrt(2/3)*Vo) / (Kspan*(KTA(k) - CF(k)*KTT));
            IR = Re*IC1/Rg;
            Rg0 = Rg;
            Rg = 50/sqrt(1 + IR/25.5);
            tf0 = tf;
            tf = 0.207*IC1^0.53;
        end 
    end

    ph = [0 -2*pi/3 2*pi/3 2*pi/3 -2*pi/3 0];
    for i = 1:nc
        for f = 0:360
            IC(i,f+1) = (CFO_NS(i) - Vo*sin(f*pi/180+ph(i)))/ VI(i);
            % Probability 
            if IC(i) < 20
                P(i,f+1) = 1 - logncdf(IC(i,f+1),log(61.1),1.33);
            else
                P(i,f+1) = 1 - logncdf(IC(i,f+1),log(33.3),0.605);
            end
        end
    end
    %% Lightning 
    PBF = max(P);
    Ng = 0.04*T^(1.25); % Ground Flash Density 
    %Ng = 3.6;
    NL = Ng*(28*HT^(0.6)+b)/10;
    BFR = 0.6*NL*sum(PBF)/length(PBF);
%     disp('Total backflashover per 100 km per year - CIGRE')
%     disp(BFR)
end

