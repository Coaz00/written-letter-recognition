clear;
close all;
clc;

N = 1000;

M1 = [0;7.5]; S1 = eye(2);
M2 = [-7.5;-7.5]; S2 = eye(2);
M3 = [7.5;-7.5]; S3 = eye(2);

K1 = mvnrnd(M1,S1,N);
K2 = mvnrnd(M2,S2,N);
K3 = mvnrnd(M3,S3,N);

K = [K1;K2;K3];

figure(1)
hold all;

scatter(K1(:,1), K1(:,2), 'ro');
scatter(K2(:,1), K2(:,2), 'bo');
scatter(K3(:,1), K3(:,2), 'go');

xlabel('x1[]')
ylabel('x2[]')
title('Raspored odabiraka')
legend('K1','K2','K3')

%% Resupstitucija

M1_ = mean(K1); S1_ = cov(K1);
M2_ = mean(K2); S2_ = cov(K2);
M3_ = mean(K3); S3_ = cov(K3);

N1 = length(K1);
N2 = length(K2);
N3 = length(K3);
N = length(K);

p1_ = N1/N;
p2_ = N2/N;
p3_ = N3/N;

%% Izmedju K1 i K2

v0_opt_s = [];
Neps_s = [];

s = 0:1e-3:1;

for i = 1:length(s)
     V = ((s(i)*S1_+(1-s(i))*S2_)^(-1))*(M2_'-M1_');
     
     Y1 = V'*K1';
     Y2 = V'*K2';
     Y = [Y1 Y2];
     Y = sort(Y);
     
     v0 = [];
     Neps = [];
     
     for j =1:length(Y)-1
         v0(j) =-(Y(j)+Y(j+1))/2;
         Neps(j) = 0;
         for k = 1:N1
             if Y1(k) > -v0(j)
                 Neps(j) = Neps(j) + 1;
             end
             if Y2(k) < -v0(j)
                 Neps(j) = Neps(j) + 1;
             end
         end
     end
     [Neps_s(i),idx] = min(Neps);
     v0_opt_s(i) = v0(idx);
end

[Neps_opt,idx] = min(Neps_s);
v0_opt = v0_opt_s(idx);
s_opt = s(idx);

V12 = (s_opt*S1_+(1-s_opt)*S2_)^(-1)*(M2_'-M1_');
v012 = v0_opt;

x1 = -10:0.1:0;
x2 = -(v0_opt+V(1)*x1)/V(2);

figure(1)
plot(x1,x2,'k');

%% Izmedju K2 i K3

v0_opt_s = [];
Neps_s = [];
s = 0:1e-3:1;

for i = 1:length(s)
     V = ((s(i)*S2_+(1-s(i))*S3_)^(-1))*(M3_'-M2_');
     
     Y2 = V'*K2';
     Y3 = V'*K3';
     Y = [Y2 Y3];
     Y = sort(Y);
     
     v0 = [];
     Neps = [];
     for j =1:length(Y)-1
         v0(j) =-(Y(j)+Y(j+1))/2;
         Neps(j) = 0;
         for k = 1:N2
             if Y2(k) > -v0(j)
                 Neps(j) = Neps(j) + 1;
             end
             if Y3(k) < -v0(j)
                 Neps(j) = Neps(j) + 1;
             end
         end
     end
     [Neps_s(i),idx] = min(Neps);
     v0_opt_s(i) = v0(idx);
end

[Neps_opt,idx] = min(Neps_s);
v0_opt = v0_opt_s(idx);
s_opt = s(idx);

V23 = (s_opt*S2_+(1-s_opt)*S3_)^(-1)*(M3_'-M2_');
v023 = v0_opt;

x1 = -0.2:0.01:0.1;
x2 = -(v0_opt+V(1)*x1)/V(2);

figure(1)
plot(x1,x2,'k');

%% Izmedju K3 i K1
v0_opt_s = [];
Neps_s = [];

s = 0:1e-3:1;

for i = 1:length(s)
     V = ((s(i)*S3_+(1-s(i))*S1_)^(-1))*(M1_'-M3_');
     
     Y3 = V'*K3';
     Y1 = V'*K1';
     Y = [Y3 Y1];
     Y = sort(Y);
     
     v0 = [];
     Neps = [];
     for j =1:length(Y)-1
         v0(j) =-(Y(j)+Y(j+1))/2;
         Neps(j) = 0;
         for k = 1:N3
             if Y3(k) > -v0(j)
                 Neps(j) = Neps(j) + 1;
             end
             if Y1(k) < -v0(j)
                 Neps(j) = Neps(j) + 1;
             end
         end
     end
     [Neps_s(i),idx] = min(Neps);
     v0_opt_s(i) = v0(idx);
end

[Neps_opt,idx] = min(Neps_s);
v0_opt = v0_opt_s(idx);
s_opt = s(idx);

V31 = (s_opt*S3_+(1-s_opt)*S1_)^(-1)*(M1_'-M3_');
v031 = v0_opt;
x1 = 0:0.1:10;
x2 = -(v0_opt+V(1)*x1)/V(2);

figure(1)
plot(x1,x2,'k');

%% Spajanje
pred = zeros(1,length(K));
for i = 1:length(K)
    glasovi = zeros(1,3);
    
    Y12 = V12'*K(i,:)';
    if Y12 > -v012
        glasovi(2) = glasovi(2) + 1;
    else
        glasovi(1) = glasovi(1) + 1;
    end
    
        Y23 = V23'*K(i,:)';
    if Y23 > -v023
        glasovi(3) = glasovi(3) + 1;
    else
        glasovi(2) = glasovi(2) + 1;
    end
    
        Y31 = V31'*K(i,:)';
    if Y31 > -v031
        glasovi(1) = glasovi(1) + 1;
    else
        glasovi(3) = glasovi(3) + 1;
    end
    
    [~,pred(i)] = max(glasovi);
end

true = [ones(1,N1) ones(1,N2)*2 ones(1,N3)*3];
C = confusionmat(true,pred);

%%
figure(2)
confusionchart(C,{'K1','K2','K3'})

