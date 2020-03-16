%% ZAD1
%wyznaczanie y_pp dla u_pp = 1.1
u_pp = 1.1;

t_sim = 300;
y = [0;0];

for k=2:t_sim
    y_temp = symulacja_obiektu6Y(u_pp,u_pp,y(k),y(k-1));
    y = [y;y_temp];
end

figure(1)
stairs(0:t_sim,y)

y_pp = y(end);

%% ZAD2
t_sim2 = 300;

%konstruowanie sygnaÂ³Ã³w sterujÂ¹cych
step_tim = 0;

u_base = ones(1,step_tim)*u_pp;
u_step_temp = ones(1,t_sim2 - step_tim)*1.3;
u_step2_temp = ones(1,t_sim2 - step_tim)*0.9;

u_step = [u_base, u_step_temp];
u_step2 = [u_base, u_step2_temp];

y = ones(t_sim2, 1)*y_pp;
y2 = ones(t_sim2, 1)*y_pp;

for k = 3:t_sim2
    if k-11 <= 0
        y(k) = symulacja_obiektu6Y(u_pp,u_pp,y(k-1),y(k-2));
        y2(k) = symulacja_obiektu6Y(u_pp,u_pp,y2(k-1),y2(k-2));
    else
        y(k) = symulacja_obiektu6Y(u_step(k-10),u_step(k-11),y(k-1),y(k-2));
        y2(k) = symulacja_obiektu6Y(u_step2(k-10),u_step2(k-11),y2(k-1),y2(k-2));
    end
end

figure(2)
stairs(0:t_sim2-1, y)
hold on;
stairs(0:t_sim2-1, y2)

%% ZAD 3  -- chyba ok ale ktos moze sprawdzic czy nie pomieszalem od ktorej chwili powinna byc odpowiedz skokowa

% TODO: naprawic wektory

t_sim3 = 600;

u = ones(t_sim3, 1) * u_pp;
y = zeros(t_sim3, 1);

% Stabilizuje obiekt w punkcie pracy po czym wykonuje skok w ramach
% ograniczen
for k = 3:t_sim3

    if k-11 <= 0
        y(k) = symulacja_obiektu6Y(0,0,y(k-1),y(k-2));
    else
        y(k) = symulacja_obiektu6Y(u(k-10),u(k-11),y(k-1),y(k-2));
    end
    
    if k == 200
        y_pp = y(k);
        u(200:t_sim3) = 1.4;
    end
end

figure(3)
stairs(0:t_sim3-1, y)

s = (y(201:t_sim3) - y_pp) / (abs(0.3));
figure(4)
stairs(0:size(s)-1, s)


%% ZAD 4 - PID

du_max = 0.1;
u_min = 0.6;
u_max = 1.6;
u_pp = 1.1;
y_pp = 2.5;

% Nastawy
T=1;
K = 2.5;
Ti = 60;
Td = 0;

r0 = K * (1 + (T/(2*Ti)) + (Td/T));
r1 = K * ((T/(2*Ti)) - (2*Td/T) - 1);
r2 = K*Td/T;


t_sim4 = 600;
y_zad = ones(t_sim4, 1) * 2.7;
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

figure(5)
stairs(0:t_sim4-1, y)





%% ZAD 5 - DMC moze dziala a moze nie xd (nie dobrane nastawy)


D = 200;
N = 15;
Nu = 15;
lambda = 1;

% Macierz M
M = zeros(N,Nu);
for i = 1:size(M,1)
    for j = 1:size(M,2)
        if i>=j
            M(i,j) = s(i-j+1);
        end
    end
end

% Macierz Mp
Mp = zeros(N,D-1);
for i = 1:size(Mp,1)
    for j = 1:size(Mp,2)
        if i+j<D
            Mp(i,j) = s(i+j) - s(j);
        else
            Mp(i,j) = s(D) - s(j);
        end
    end
end

K = ((M'*M + lambda*eye(Nu))^-1)*M';
ke = sum(K(1,:));

ku = zeros(D-1,1);
for i = 1:D-1
    ku(i) = K(1,:) * Mp(:,i);
end


t_sim4 = 600;
y_zad = ones(t_sim4, 1) * 2.7;
y = ones(t_sim4, 1) * y_pp;
u = ones(t_sim4, 1) * u_pp;
du = zeros(t_sim4, 1);
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
    
    current_sum = 0;
    for i = 1:D-1
        if k-j > 1
            current_sum = current_sum + ku(i) * du(k-j);
        end
    end
    duk = ke*e(k) - current_sum;
    
    % Na³o¿enie ograniczeñ
    if duk>du_max
        duk = du_max;
    elseif duk<-du_max
        duk = -du_max;
    end
    du(k) = duk;
    uk = u(k-1) + duk;
    if uk>u_max
        uk = u_max;
    elseif uk<u_min
        uk = u_min;
    end
    u(k) = uk;
end

figure(6)
stairs(0:t_sim4-1, y)

