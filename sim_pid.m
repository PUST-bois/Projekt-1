% Ograniczenia
du_max = 0.1;
u_min = 0.6;
u_max = 1.6;
% Punkt pracy
u_pp = 1.1;
y_pp = 2.5;
% WskaŸnik jakoœci - b³¹d œredniokwadratowy
piderr = 0;
% Okres regulacji
T=1;
% Nastawy PID
params = [4.4817,   14.4945,    6.8265];
K = params(1);
Ti = params(2);
Td = params(3);

r0 = K * (1 + (T/(2*Ti)) + (Td/T));
r1 = K * ((T/(2*Ti)) - (2*Td/T) - 1);
r2 = K*Td/T;

% Czas symulacji
t_sim = 800;
% Zadana trajektoria
y_zad = ones(t_sim, 1) * 2.7;
y_zad(100:250) = 2.9;
y_zad(250:400) = 2.7;
y_zad(400:600) = 2.4;
%Inicjalizacja wektorów y, u, e
y = ones(t_sim, 1) * y_pp;
u = ones(t_sim, 1) * u_pp;
e = zeros(t_sim, 1);

% Pêtla, w której odbywa siê symulacja
for k = 3:t_sim
    % Symulacja obiektu
    if k-11 <= 0
        y(k) = symulacja_obiektu6Y(u_pp,u_pp,y(k-1),y(k-2));
    else
        y(k) = symulacja_obiektu6Y(u(k-10),u(k-11),y(k-1),y(k-2));
    end
    % Obliczenie uchybu
    e(k) = y_zad(k) - y(k);
    % Obliczenie wskaŸnika jakoœci
    piderr = piderr + e(k)^2;
    % Obliczenie przyrostu sterowania
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

% Wyœwietlenie wyników symulacji
figure(1)
stairs(0:t_sim-1, y, 'b')
hold on
stairs(0:t_sim-1, y_zad, '--r')

hold off