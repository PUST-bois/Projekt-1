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

%konstruowanie sygna³ów steruj¹cych
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

%% ZAD 3  -- tu zle, trzeba zobaczyć jakie y ustali się dla u = 0
t_sim3 = 200;
step_tim = 0;

u_base = ones(1,step_tim)*0;
u_step_temp = ones(1,t_sim3 - step_tim);


u_step = [u_base, u_step_temp];


y = ones(t_sim3, 1)*0;

for k = 3:t_sim3
    if k-11 <= 0
        y(k) = symulacja_obiektu6Y(0,0,y(k-1),y(k-2));
    else
        y(k) = symulacja_obiektu6Y(u_step(k-10),u_step(k-11),y(k-1),y(k-2));
    end
end


figure(3)
stairs(0:t_sim3-1, y)

s = y;
