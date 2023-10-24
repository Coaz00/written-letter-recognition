clear;
close all;
clc;

N = 500;
P11 = 0.7;
P12 = 1 - P11;
P21 = 0.6;
P22 = 1 - P21;

M11 = [-3;-4]; S11 = [3 1; 1 3];
M12 = [-3;4]; S12 = [1.5 0.8; 0.8 1.5];
M21 = [3;4]; S21 = [2 -0.5; -0.5 2];
M22 = [3;-4]; S22 = [1.5 0.6; 0.6 1.5];

%% Teorijska fgv
x = -10:0.1:10;
y = -10:0.1:10;

f11 = zeros(length(x), length(y));
f12 = zeros(length(x), length(y));
f21 = zeros(length(x), length(y));
f22 = zeros(length(x), length(y));
f1 = zeros(length(x),length(y));
f2 = zeros(length(x),length(y));

for i = 1:length(x)
    for j = 1:length(y)
        X = [y(j); x(i)];
        f11(i,j) = 1/((2*pi)^2*(det(S11)))^0.5*exp(-0.5*(X-M11)'*inv(S11)*(X-M11));
        f12(i,j) = 1/((2*pi)^2*(det(S12)))^0.5*exp(-0.5*(X-M12)'*inv(S12)*(X-M12));
        f21(i,j) = 1/((2*pi)^2*(det(S21)))^0.5*exp(-0.5*(X-M21)'*inv(S21)*(X-M21));
        f22(i,j) = 1/((2*pi)^2*(det(S22)))^0.5*exp(-0.5*(X-M22)'*inv(S22)*(X-M22));
        f1(i,j) = P11*f11(i,j) + P12*f12(i,j);
        f2(i,j) = P21*f21(i,j) + P22*f22(i,j);
    end
end



%% Generisanje odabiraka

pom = rand(N,1);
K11 = mvnrnd(M11, S11, N);
K12 = mvnrnd(M12, S12, N);
K1 = (pom<P11).*K11 + (pom>=P11).*K12;

pom = rand(N,1);
K21 = mvnrnd(M21, S21, N);
K22 = mvnrnd(M22, S22, N);
K2 = (pom<P21).*K21 + (pom>=P21).*K22;

%% Iscrtavnaje
figure(1)
hold all;
scatter(K1(:,1), K1(:,2), 'ro');
scatter(K2(:,1), K2(:,2), 'bo');
title('Raspored odabiraka')
xlabel('x1[]')
ylabel('x2[]')
legend('K1','K2')

figure(2)
subplot(1,2,1)
mesh(x,y,f1);
title('Teorijska fgv')
xlabel('x1[]')
ylabel('x2[]')
xlim([-10 10])
subplot(1,2,2)
histogram2(K1(:,1),K1(:,2),'Normalization','probability')
title('Histogram odabiraka')
xlabel('x1[]')
ylabel('x2[]')
xlim([-10 10])

figure(3)
subplot(1,2,1)
mesh(x,y,f2);
title('Teorijska fgv')
xlabel('x1[]')
ylabel('x2[]')
xlim([-10 10])
subplot(1,2,2)
histogram2(K2(:,1),K2(:,2),'Normalization','probability')
title('Histogram odabiraka')
xlabel('x1[]')
ylabel('x2[]')
xlim([-10 10])

%% Bayesov test minimalne greske

h = -log(f1./f2);
T = 0; % p1 = p2 -> T = ln(p2/p1) = 0

figure(4)
contour(x,y,h - T,[0 0],'g','LineWidth',1);
hold all;
scatter(K1(:,1), K1(:,2), 'ro');
scatter(K2(:,1), K2(:,2), 'bo');
title('Bayes-ov test')
xlabel('x1[]')
ylabel('x2[]')
legend('Diskriminaciona kriva','K1','K2')


Xs = [K1',K2'];
X_true = [ones(1,N), ones(1,N)*2]; %stvarne klase
X_pred = zeros(size(X_true));
for i = 1:length(Xs)
    X = Xs(:,i);
    f11_ = 1/((2*pi)^2*(det(S11)))^0.5*exp(-0.5*(X-M11)'*inv(S11)*(X-M11));
    f12_ = 1/((2*pi)^2*(det(S12)))^0.5*exp(-0.5*(X-M12)'*inv(S12)*(X-M12));
    f21_ = 1/((2*pi)^2*(det(S21)))^0.5*exp(-0.5*(X-M21)'*inv(S21)*(X-M21));
    f22_ = 1/((2*pi)^2*(det(S22)))^0.5*exp(-0.5*(X-M22)'*inv(S22)*(X-M22));
    f1_ = P11*f11_ + P12*f12_;
    f2_ = P21*f21_ + P22*f22_;
    if(f1_>f2_)
        X_pred(i) = 1;
    else
        X_pred(i) = 2;
    end
end

C_bayes = confusionmat(X_true,X_pred);
figure(10)
confusionchart(C_bayes,{'K1','K2'})
e1_bayes = C_bayes(1,2)/sum(C_bayes(1,:))*100; 
e2_bayes = C_bayes(2,1)/sum(C_bayes(2,:))*100;
e_bayes = 0.5*(e1_bayes + e2_bayes);

%% Test minimalne cene

c11 = 0;
c22 = 0;
c21 = 4;
c12 = 1;

T = (c12 - c22)/(c21-c11);

h = -log(f1./f2)+log(T); 

figure(5)
contour(x,y,h,[0 0],'g');
hold all;
scatter(K1(:,1), K1(:,2), 'ro');
scatter(K2(:,1), K2(:,2), 'bo');
title('Test minimalne cene')
xlabel('x1[]')
ylabel('x2[]')
legend('Diskriminaciona kriva','K1','K2')


Xs = [K1',K2'];
X_true = [ones(1,N), ones(1,N)*2]; %stvarne klase
X_pred = zeros(size(X_true));
for i = 1:length(Xs)
    X = Xs(:,i);
    f11_ = 1/((2*pi)^2*(det(S11)))^0.5*exp(-0.5*(X-M11)'*inv(S11)*(X-M11));
    f12_ = 1/((2*pi)^2*(det(S12)))^0.5*exp(-0.5*(X-M12)'*inv(S12)*(X-M12));
    f21_ = 1/((2*pi)^2*(det(S21)))^0.5*exp(-0.5*(X-M21)'*inv(S21)*(X-M21));
    f22_ = 1/((2*pi)^2*(det(S22)))^0.5*exp(-0.5*(X-M22)'*inv(S22)*(X-M22));
    f1_ = P11*f11_ + P12*f12_;
    f2_ = P21*f21_ + P22*f22_;
    if f1_/f2_ > T
        X_pred(i) = 1;
    else
        X_pred(i) = 2;
    end
end

C_cena = confusionmat(X_true,X_pred);
figure(11)
confusionchart(C_cena,{'K1','K2'})
e1_cena = C_cena(1,2)/sum(C_cena(1,:))*100; 
e2_cena = C_cena(2,1)/sum(C_cena(2,:))*100;
e_cena = 0.5*(e1_cena + e2_cena);

%% Teorijska greska
e1 = 0;
e2 = 0;
for i = 1:length(x)-1
    for j = 1:length(y)-1
        X = [x(i) y(j)]';
        h = -log(f1(i,j)/f2(i,j));
        if(h<0)
            e2 = e2 + 0.1*0.1*((f2(i,j)+f2(i+1,j)+f2(i,j+1)+f2(i+1,j+1))/4);
        else 
            e1 = e1+0.1*0.1*((f2(i,j)+f2(i+1,j)+f2(i,j+1)+f2(i+1,j+1))/4);
        end
    end
end

e = 1/2*(e1+e2);

%% Neyman-Pearson-ov test

h = -log(f1./f2);

br = 0;
mi_skup = 0.01:0.01:10;
for mi = mi_skup
    br = br+1;
    Eps0(br) = 0;
    for i = 1:length(x)-1
        for j = 1:length(y)-1
            if (h(i,j)<-log(mi))
                Eps0(br) = Eps0(br) + 0.1*0.1*((f2(i,j)+f2(i+1,j)+f2(i,j+1)+f2(i+1,j+1))/4); 
            end
        end
    end
end

figure(6);
plot(mi_skup,Eps0);
xlabel('\mu');
ylabel('\epsilon_0');
title('Zavisnost $\epsilon_0$ od $\mu$','Interpreter','Latex')

tol = 0.1;
E0 = 3;
mi = mi_skup(find(abs(Eps0*100 - E0) < tol));
mi = mi(1);

figure(7)
contour(x,y,h+log(mi),[0 0],'g');
hold all;
scatter(K1(:,1), K1(:,2), 'ro');
scatter(K2(:,1), K2(:,2), 'bo');
title('Neyman-Pearson-ov test')
xlabel('x1[]')
ylabel('x2[]')
legend('Diskriminaciona kriva','K1','K2')


Xs = [K1',K2'];
X_true = [ones(1,N), ones(1,N)*2];
X_pred = zeros(size(X_true));
for i = 1:length(Xs)
    X = Xs(:,i);
    f11_ = 1/(2*pi*(det(S11)))^0.5*exp(-0.5*(X-M11)'*inv(S11)*(X-M11));
    f12_ = 1/(2*pi*(det(S12)))^0.5*exp(-0.5*(X-M12)'*inv(S12)*(X-M12));
    f21_ = 1/(2*pi*(det(S21)))^0.5*exp(-0.5*(X-M21)'*inv(S21)*(X-M21));
    f22_ = 1/(2*pi*(det(S22)))^0.5*exp(-0.5*(X-M22)'*inv(S22)*(X-M22));
    f1_ = P11*f11_ + P12*f12_;
    f2_ = P21*f21_ + P22*f22_;
    h_ = -log(f1_./f2_);
    if h_ < -log(mi)
        X_pred(i) = 1;
    else
        X_pred(i) = 2;
    end
end

C_np = confusionmat(X_true,X_pred);
figure(12)
confusionchart(C_np,{'K1','K2'})
e1_np = C_np(1,2)/sum(C_np(1,:))*100; 
e2_np = C_np(2,1)/sum(C_np(2,:))*100;
e_np = 1/2*(e1_np+e2_np);



%% Waldov test
br_koraka = zeros(length(1e-10:1e-11:1e-8));
cnt = 0;
for eps2 =1e-10:1e-11:1e-8
    cnt = cnt + 1;
    eps1 = 1e-9;
    
    A = (1-eps1)/eps2;
    B = eps1/(1-eps2);
    
    a = -log(A);
    b = -log(B);
    
    Sm = 0;
    Seq = [K1(1:0.7*length(K1),:); K2(1:0.3*length(K2),:)];
    Seq = Seq(randperm(length(Seq)),:);
    
    for i = 1:length(Seq)
        X = Seq(i,:);
        X = X';
        f11_ = 1/(2*pi*(det(S11)))^0.5*exp(-0.5*(X-M11)'*inv(S11)*(X-M11));
        f12_ = 1/(2*pi*(det(S12)))^0.5*exp(-0.5*(X-M12)'*inv(S12)*(X-M12));
        f21_ = 1/(2*pi*(det(S21)))^0.5*exp(-0.5*(X-M21)'*inv(S21)*(X-M21));
        f22_ = 1/(2*pi*(det(S22)))^0.5*exp(-0.5*(X-M22)'*inv(S22)*(X-M22));
        f1_ = P11*f11_ + P12*f12_;
        f2_ = P21*f21_ + P22*f22_;
        
        Sm = Sm - log(f1_./f2_);
        
        if Sm < a
            pred = 1;
            break;
        elseif Sm > b
            pred = 2;
            break;
        end
    end
    
    br_koraka(cnt) = i;
end

figure(20)
plot(1e-10:1e-11:1e-8,br_koraka)
title('Zavisnost broja koraka od greske druge vrste')
xlabel('\epsilon_2')



