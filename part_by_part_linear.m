clear;
close all;
clc;

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

%% Izmedju K1 i K2
U = [-ones(1,N) ones(1,N); -K1' K2'];
G = [ones(N,1);2*ones(N,1)];
W = ((U*U')^(-1))*U*G;
v0 = W(1);
V(1) = W(2);
V(2) = W(3);

x1 = -10:0.1:0;
x2 = -(v0+V(1)*x1)/V(2);
figure(1)
plot(x1,x2,'k');

v012 = v0;
V12 = V';

%% Izmedju K2 i K3

U = [-ones(1,N) ones(1,N); -K2' K3'];
G = [ones(N,1);2*ones(N,1)];
W = ((U*U')^(-1))*U*G;
v0 = W(1);
V(1) = W(2);
V(2) = W(3);

x1 = -0.2:0.01:0.1;
x2 = -(v0+V(1)*x1)/V(2);
figure(1)
plot(x1,x2,'k');

v023 = v0;
V23 = V';


%% Izmedju K3 i K1

U = [-ones(1,N) ones(1,N); -K3' K1'];
G = [ones(N,1);2*ones(N,1)];
W = ((U*U')^(-1))*U*G;
v0 = W(1);
V(1) = W(2);
V(2) = W(3);

x1 = 0:0.1:10;
x2 = -(v0+V(1)*x1)/V(2);
figure(1)
plot(x1,x2,'k');

v031 = v0;
V31 = V';

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
true = [ones(1,N) ones(1,N)*2 ones(1,N)*3];
C = confusionmat(true,pred);

%% 
figure(2)
confusionchart(C,{'K1','K2','K3'})

