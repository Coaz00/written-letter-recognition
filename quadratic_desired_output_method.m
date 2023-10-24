clear;
close all;
clc;

N = 1000;

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
xlabel('x1[]')
ylabel('x2[]')
title('Raspored odabiraka')
K = [K1;K2];

%% Kvadratni klasifikator na bazi zeljenog izlaza
K1=K1';
K2=K2';

G = [ones(N,1); ones(N,1)];
U = [-ones(1,N) ones(1,N); -K1 K2; -(K1(1,:)).^2 (K2(1,:)).^2; -(K1(2,:)).^2 (K2(2,:)).^2; -2*K1(1,:).*K1(2,:) 2*K2(1,:).*K2(2,:)];

W = (U*U')^(-1)*U*G;

v0 = W(1);
V(1) = W(2);
V(2) = W(3);
Q(1,1) = W(4);
Q(2,2) = W(5);
Q(1,2) = W(6);


V = V';

x1 = -6:0.1:6;
x2 = -6:0.1:6;
h = zeros(length(x1), length(x2));

for i = 1:length(x1)
    for j = 1:length(x2)
        h(i,j)=v0+V(1)*x1(i)+V(2)*x2(j)+Q(1,1)*(x1(i))^2+Q(2,2)*(x2(j))^2+Q(1,2)*x1(i)*x2(j);
    end
end

figure(1)
hold all;
contour(x1,x2,h,[0 0]);
legend('K1','K2','klas. kriva')

%% Klasifikacija
pred = zeros(1,length(K));
true = [ones(1,N) 2*ones(1,N)];
for i = 1:length(K)
    X = K(i,:);
    h = v0+V(1)*X(1)+V(2)*X(2)+Q(1,1)*(X(1))^2+Q(2,2)*(X(2))^2+Q(1,2)*X(1)*X(2);
    if h < 0
        pred(i) = 1;
    else
        pred(i) = 2;
    end
end

C = confusionmat(true,pred);

%%
confusionchart(C,{'K1','K2'})
