clear;
close all;
clc;

N = 500;

M1 = [7.5;7.5]; S1 = eye(2);
M2 = [-7.5;7.5]; S2 = eye(2);
M3 = [-7.5;-7.5]; S3 = eye(2);
M4 = [7.5;-7.5]; S4 = eye(2);

K1 = mvnrnd(M1,S1,N);
K2 = mvnrnd(M2,S2,N);
K3 = mvnrnd(M3,S3,N);
K4 = mvnrnd(M4,S4,N);

figure(1)
hold all;
scatter(K1(:,1),K1(:,2),'ro')
scatter(K2(:,1),K2(:,2),'go')
scatter(K3(:,1),K3(:,2),'bo')
scatter(K4(:,1),K4(:,2),'yo')
legend('K1','K2','K3','K4')
title('Raspored odabiraka')
xlabel('x1[]')
ylabel('x2[]')

%% C-mean

%inicijalna kalsterizacija
omega = randi([1 4],1,4*N);
omega_novo = zeros(1,4*N);
K = [K1;K2;K3;K4];

K1_ = K(omega == 1,:);
K2_ = K(omega == 2,:);
K3_ = K(omega == 3,:);
K4_ = K(omega == 4,:);

figure(2)
hold all;
scatter(K1_(:,1),K1_(:,2),'ro')
scatter(K2_(:,1),K2_(:,2),'go')
scatter(K3_(:,1),K3_(:,2),'bo')
scatter(K4_(:,1),K4_(:,2),'yo')
legend('K1','K2','K3','K4')
title('Inicijalna klasterizacija')
xlabel('x1[]')
ylabel('x2[]')
cnt = 0;
while(1)
    cnt = cnt+1;
    flag = 0;    
    
    K1_ = K(omega == 1,:);
    K2_ = K(omega == 2,:);
    K3_ = K(omega == 3,:);
    K4_ = K(omega == 4,:);
    
    M1_ = mean(K1_);
    M2_ = mean(K2_);
    M3_ = mean(K3_);
    M4_ = mean(K4_);
    
    M_ = [M1_;M2_;M3_;M4_];
    
    for i = 1:length(K)
        X = K(i,:);
        d_min = inf;
        for j = 1:4
            %d = (X(1) - M_(j,1))^2 + (X(2) - M_(j,2))^2;
            d = (X'-M_(j,:)')'*(X'-M_(j,:)');
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


figure(3)
hold all;
scatter(K1_(:,1),K1_(:,2),'ro')
scatter(K2_(:,1),K2_(:,2),'go')
scatter(K3_(:,1),K3_(:,2),'bo')
scatter(K4_(:,1),K4_(:,2),'yo')
legend('K1','K2','K3','K4')
title('C-mean klasterizacija')
xlabel('x1[]')
ylabel('x2[]')

%% Maximum Likelihood
K = [K1;K2;K3;K4];
T = 0.01;

omega = [ones(1,0.75*N) randi([1,4],1,0.25*N) 2*ones(1,0.75*N) randi([1,4],1,0.25*N) 3*ones(1,0.75*N) randi([1,4],1,0.25*N) 4*ones(1,0.75*N) randi([1,4],1,0.25*N)];

K1_ = K(omega == 1,:);
K2_ = K(omega == 2,:);
K3_ = K(omega == 3,:);
K4_ = K(omega == 4,:);

figure(4)
hold all;
scatter(K1_(:,1),K1_(:,2),'ro')
scatter(K2_(:,1),K2_(:,2),'go')
scatter(K3_(:,1),K3_(:,2),'bo')
scatter(K4_(:,1),K4_(:,2),'yo')
legend('K1','K2','K3','K4')
title('Inicijalna klasterizacija')
xlabel('x1[]')
ylabel('x2[]')

M1_ = mean(K1_);
M2_ = mean(K2_);
M3_ = mean(K3_);
M4_ = mean(K4_);

M_ = [M1_;M2_;M3_;M4_];

P1_ = length(K1_)/length(K);
P2_ = length(K2_)/length(K);
P3_ = length(K3_)/length(K);
P4_ = length(K4_)/length(K);

P_ = [P1_;P2_;P3_;P4_];

S1_ = cov(K1_);
S2_ = cov(K2_);
S3_ = cov(K3_);
S4_ = cov(K4_);

S_ = [S1_;S2_;S3_;S4_];


f_ = zeros(1,4);
q = zeros(1,4);

Q_old = zeros(length(K),4);

for i = 1:length(K)
    X = K(i,:)';
    for j = 1:4
        f_(j) = 1/((2*pi)^2*(det(S_(j))))^0.5*exp(-0.5*(X-M_(j,:)')'*(S_(2*(j-1)+1:2*j,:))^(-1)*(X-M_(j,:)'));
    end
    for j = 1:4
        q(j) = P_(j)*f_(j)/(f_*P_);
    end
    Q_old(i,:) = q;
end
cnt = 0;
while(cnt < 100)
    cnt = cnt + 1;
    for i  = 1:4
        P_(i) = 1/length(K)*sum(Q_old(:,i));
        
        sum_M = zeros(2,1);
        for j = 1:length(K)
            sum_M = sum_M + Q_old(j,i)*K(j,:)';
        end
        M_(i,:) = (1/P_(i)/length(K)*sum_M)';
        
        sum_S = zeros(2,2);
        for j = 1:length(K)
            sum_S = sum_S + Q_old(j,i)*(K(j,:)'-M_(i,:)')*(K(j,:)'-M_(i,:)')';
        end
        S_(2*(i-1)+1:2*i,:) = 1/P_(i)/length(K)*sum_S;
    end
    
    Q_new = zeros(length(K),4);
    
    for i = 1:length(K)
        X = K(i,:)';
        for j = 1:4
            f_(j) = 1/((2*pi)^2*(det(S_(j))))^0.5*exp(-0.5*(X-M_(j,:)')'*(S_(2*(j-1)+1:2*j,:))^(-1)*(X-M_(j,:)'));
        end
        for j = 1:4
            q(j) = P_(j)*f_(j)/(f_*P_);
        end
        Q_new(i,:) = q;
    end
    
    omega_new = zeros(1,length(K));
    for i = 1:length(K)
        [~,omega_new(i)] = max(Q_new(i,:));
    end

    
    if max(max(abs(Q_new - Q_old))) < T
        break;
    end
    
    omega = omega_new;
    Q_old = Q_new;
    
end

K1_ = K(omega == 1,:);
K2_ = K(omega == 2,:);
K3_ = K(omega == 3,:);
K4_ = K(omega == 4,:);

figure(5)
hold all;
scatter(K1_(:,1),K1_(:,2),'ro')
scatter(K2_(:,1),K2_(:,2),'go')
scatter(K3_(:,1),K3_(:,2),'bo')
scatter(K4_(:,1),K4_(:,2),'yo')
legend('K1','K2','K3','K4')
title('ML klasterizacija')
xlabel('x1[]')
ylabel('x2[]')