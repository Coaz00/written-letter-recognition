clear;
close all;
clc;

N = 500;

%KRUG
R = 2*rand(1,N);
theta = 2*pi*rand(1,N);

K1 = [R.*cos(theta);R.*sin(theta)]';

%PRSTEN
R = 7 + 2*rand(1,N);
theta = 2*pi*rand(1,N);

K2 = [R.*cos(theta);R.*sin(theta)]';

figure(1)
hold all;
scatter(K1(:,1),K1(:,2))
scatter(K2(:,1),K2(:,2))
title('Raspored odabiraka')
xlabel('x1[]')
ylabel('x2[]')

%inicijalna kalsterizacija
omega = randi([1 2],1,2*N);
omega_novo = zeros(1,2*N);

K = [K1;K2];

while(1)
    flag = 0;    
    
    K1_ = K(omega == 1,:);
    K2_ = K(omega == 2,:);
    
    M1_ = mean(K1_);
    M2_ = mean(K2_);
    
    M_ = [M1_;M2_];
    
    S1_ = cov(K1_);
    S2_ = cov(K2_);
    
    S_ = [S1_;S2_];
   
    for i = 1:length(K)
        X = K(i,:);
        d_min = inf;
        for j = 1:2
            d = (X(1) - M_(j,1))^2 + (X(2) - M_(j,2))^2;
            if d < d_min
                d_min = d;
                idx = j;
            end
        end
        omega_novo(i) = idx;
        if idx ~= omega(i)
            flag = 1;
        end
    end
    
    if flag == 0
        break;
    end
    
    omega = omega_novo;
    
end

figure(2)
hold all;
scatter(K1_(:,1),K1_(:,2),'ro')
scatter(K2_(:,1),K2_(:,2),'go')
title('Inicijalna klasterizacija')
xlabel('x1[]')
ylabel('x2[]')

%% kvadratna dekompozicija

omega = randi([1 2],1,2*N);

cnt = 0;
while(1)
    cnt = cnt + 1;
    flag = 0;
    
    K1_ = K(omega == 1,:);
    K2_ = K(omega == 2,:);
    
    P1_ = length(K1_)/length(K);
    P2_ = length(K2_)/length(K);
    P_ = [P1_;P2_];
    
    M1_ = mean(K1_);
    M2_ = mean(K2_);
    M_ = [M1_;M2_];
    
    S1_ = cov(K1_);
    S2_ = cov(K2_);
    S_ = [S1_;S2_];
    
    for i = 1:length(K)
        X = K(i,:);
        d_min = inf;
        for j = 1:2
            d = (X'-M_(j,:)')'*(S_(2*(j-1)+1:2*j,:))^(-1)*(X'-M_(j,:)') + log(det(S_(2*(j-1)+1:2*j,:))) - log(P_(j));
            if d < d_min
                d_min = d;
                idx = j;
            end
        end
%         d1 = 0.5*(X' - M1_')'*S1_^(-1)*(X'-M1_') + 0.5*log(det(S1_)) - 0.5*log(P1_);
%         d2 = 0.5*(X' - M2_')'*S2_^(-1)*(X'-M2_') + 0.5*log(det(S2_)) - 0.5*log(P2_);
%         if d1 < d2
%             idx = 1;
%         else
%             idx = 2;
%         end
        omega_novo(i) = idx;
        if idx ~= omega(i)
            flag = 1;
        end
    end
    
    if flag == 0
        break;
    end
    
    omega = omega_novo;
 
end

figure(3)
hold all;
scatter(K1_(:,1),K1_(:,2),'ro')
scatter(K2_(:,1),K2_(:,2),'go')
title('Kvadratna dekompozicija')
xlabel('x1[]')
ylabel('x2[]')