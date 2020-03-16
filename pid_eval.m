function piderr = pid_eval(params)
du_max = 0.1;
u_min = 0.6;
u_max = 1.6;
u_pp = 1.1;
y_pp = 2.5;
piderr = 0;
% Nastawy
T=1;
K = params(1);  %8 to stale niegasnace oscylacje
Ti = params(2);%14 gud (dla K = 4 (po³owa oscylacji) i dla Td = 0 dobieramy Ti)
Td = params(3); %8 gud (dla wy¿ej dobranych K i Ti dobieramy Td)

r0 = K * (1 + (T/(2*Ti)) + (Td/T));
r1 = K * ((T/(2*Ti)) - (2*Td/T) - 1);
r2 = K*Td/T;


t_sim4 = 800;
y_zad = ones(t_sim4, 1) * 2.7;
y_zad(100:250) = 2.9;
y_zad(250:400) = 2.7;
y_zad(400:600) = 2.4;
y = ones(t_sim4, 1) * y_pp;
u = ones(t_sim4, 1) * u_pp;
e = zeros(t_sim4, 1);

% y_min = [2.02742187500021]
% y_max = [2.97257812416614]
for k = 3:t_sim4

    if k-11 <= 0
        y(k) = symulacja_obiektu6Y(u_pp,u_pp,y(k-1),y(k-2));
    else
        y(k) = symulacja_obiektu6Y(u(k-10),u(k-11),y(k-1),y(k-2));
    end
  
    e(k) = y_zad(k) - y(k);
    piderr = piderr + (e(k))^2;
    du = r2*e(k-2) + r1*e(k-1) + r0*e(k);
    
    % Na³o¿enie ograniczeñ
    if du>du_max
        du = du_max;
    elseif du<-du_max
        du = -du_max;
    end
    uk = u(k-1) + du;
    if uk>u_max
        uk = u_max;
    elseif uk<u_min
        uk = u_min;
    end
    u(k) = uk;
end

end