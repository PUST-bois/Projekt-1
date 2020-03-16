function dmcerr = dmc_eval(params)
    global s
    u_pp = 1.1;
    y_pp = 2.5;
    N = floor(params(1));
    Nu = floor(params(2));
    du_max = 0.1;
    u_max = 1.6;
    u_min = 0.6;
    lambda = 1;
    dmcerr = 0;
    D = 200;
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


    t_sim4 = 800;
    y_zad = ones(t_sim4, 1) * 2.7;
    y_zad(100:250) = 2.9;
    y_zad(250:400) = 2.7;
    y_zad(400:600) = 2.4;
    y = ones(t_sim4, 1) * y_pp;
    u = ones(t_sim4, 1) * u_pp;
    du = zeros(t_sim4, 1);
    e = zeros(t_sim4, 1);
    for k = 3:t_sim4

    if k-11 <= 0
        y(k) = symulacja_obiektu6Y(u_pp,u_pp,y(k-1),y(k-2));
    else
        y(k) = symulacja_obiektu6Y(u(k-10),u(k-11),y(k-1),y(k-2));
    end
  
    e(k) = y_zad(k) - y(k);
    dmcerr = dmcerr + e(k)^2;
    
    current_sum = 0;
    for i = 1:D-1
        if k-i > 1
            current_sum = current_sum + ku(i) * du(k-i);
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
end