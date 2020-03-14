%% ZAD1
%wyznaczanie y_pp dla u_pp = 1.1
u_pp = 1.1;

t_sim = 300;
y1 = [0;0];

for k=2:t_sim
    y_temp = symulacja_obiektu6Y(u_pp,u_pp,y1(k),y1(k-1));
    y1 = [y1;y_temp];
end

figure(1)
stairs(0:t_sim,y1)

y_pp = y1(end);

%% ZAD2
t_sim2 = 300;

%konstruowanie sygna³ów steruj¹cych
step_tim = 10;
u_base = ones(1,step_tim)*u_pp;
u_step1_temp = ones(1,t_sim2 - step_tim)*1.3;
u_step2_temp = ones(1,t_sim2 - step_tim)*0.9;

u_step1 = [u_base, u_step1_temp];
u_step2 = [u_base, u_step2_temp];

y1 = ones(t_sim, 1)*y_pp;
y2 = ones(t_sim, 1)*y_pp;

for k = 3:t_sim2
    if k-11 <= 0
        y1(k) = symulacja_obiektu6Y(u_pp,u_pp,y1(k-1),y1(k-2));
        y2(k) = symulacja_obiektu6Y(u_pp,u_pp,y2(k-1),y2(k-2));
    else
        y1(k) = symulacja_obiektu6Y(u_step1(k-10),u_step1(k-11),y1(k-1),y1(k-2));
        y2(k) = symulacja_obiektu6Y(u_step2(k-10),u_step2(k-11),y2(k-1),y2(k-2));
    end
end

figure(2)
stairs(0:t_sim2-1, y1)
hold on;
stairs(0:t_sim2-1, y2)


