clear
close all
clc

slova = load('PO_slova.mat');
imena = fieldnames(slova);
izabrana_slova = ['a','y','h','i','s','w','b','z','v','c'];
min_primeraka = inf;

for i = 1:length(izabrana_slova)
    tmp = slova.(izabrana_slova(i));
    
    if min_primeraka > length(tmp)
        min_primeraka = length(tmp);
    end
    
    tmp = cell2mat(tmp(1));
    
    tmp_vx = tmp(1,:);
    tmp_vy = tmp(2,:);
    
    tmp_x = cumsum(tmp_vx); % x(i) = x(i-1) + vx(i)*Ts -> Ts = 1
    tmp_y = cumsum(tmp_vy); % y(i) = y(i-1) + vy(i)*Ts -> Ts = 1
    
    figure(i)
    
    subplot(2,2,2)
    plot(tmp_vx)
    title('Brzina po x osi', 'Interpreter','Latex')
    xlabel('odabirak[]','Interpreter','Latex')
    ylabel('$v_x[\frac{m}{s}]$','Interpreter','Latex')
    
    subplot(2,2,4)
    plot(tmp_vy)
    xlabel('odabirak[]','Interpreter','Latex')
    title('Brzina po y osi', 'Interpreter','Latex')
    ylabel('$v_y[\frac{m}{s}]$','Interpreter','Latex')
    
    subplot(2,2,[1,3])
    plot(tmp_x,tmp_y)
    title(strcat("Izgled slova: ",upper(izabrana_slova(i))),'Interpreter','Latex')
    xlabel('$x[]$','Interpreter','Latex')
    ylabel('$y[]$','Interpreter','Latex')
end

slova_obelezja = zeros(6*length(izabrana_slova),min_primeraka);
    
 for i = 1:length(izabrana_slova)
     tmp = slova.(izabrana_slova(i));
     
     f1tmp = zeros(1,length(tmp));
     f2tmp = zeros(1,length(tmp));
     f3tmp = zeros(1,length(tmp));
     f4tmp = zeros(1,length(tmp));
     f5tmp = zeros(1,length(tmp));
     f6tmp = zeros(1,length(tmp));
     
     for j = 1:length(tmp)
        slovo = cell2mat(tmp(j));
        
        f1tmp(j) = sum(slovo(1,:));
        f2tmp(j) = max(diff(slovo(1,:))) - min(diff(slovo(1,:)));
        f3tmp(j) = max(diff(slovo(2,:))) - min(diff(slovo(2,:)));
        
        for k = 1:length(slovo(1,:))-1
            if slovo(1,k)*slovo(1,k+1) < 0
                f4tmp(j) = f4tmp(j) + 1;
            end
            if slovo(2,k)*slovo(2,k+1) < 0
                f5tmp(j) = f5tmp(j) + 1;
            end
        end
        
        f6tmp(j) = sum(slovo(2,:));
     end
     
     slova_obelezja(6*(i-1)+1,:) = f1tmp(1:min_primeraka);
     slova_obelezja(6*(i-1)+2,:) = f2tmp(1:min_primeraka);
     slova_obelezja(6*(i-1)+3,:) = f3tmp(1:min_primeraka);
     slova_obelezja(6*(i-1)+4,:) = f4tmp(1:min_primeraka);
     slova_obelezja(6*(i-1)+5,:) = f5tmp(1:min_primeraka);
     slova_obelezja(6*(i-1)+6,:) = f6tmp(1:min_primeraka);
 end
%% LDA
S = zeros(60,6);
M = zeros(60,1);

for i = 1:length(izabrana_slova)
    M(6*(i-1)+1:6*i) = mean(slova_obelezja(6*(i-1)+1:6*i,:),2);
    S(6*(i-1)+1:6*i,:) = cov(slova_obelezja(6*(i-1)+1:6*i,:)');
end

M0 = zeros(6,1);
SW = zeros(6,6);
SB = zeros(6,6);
for i = 1:length(izabrana_slova)
    M0 = M0 + M(6*(i-1)+1:6*i);
    SW = SW + S(6*(i-1)+1:6*i,:);
    SB = SB + (M(6*(i-1)+1:6*i) - M0)*(M(6*(i-1)+1:6*i) - M0)';
end
M0 = 1/length(izabrana_slova)*M0;
SW = 1/length(izabrana_slova)*SW;
SB = 1/length(izabrana_slova)*SB;

S = SW^(-1)*SB;
[F, L] = eig(S);
A = F(:,1:3);

slova_obelezja_lda = zeros(3*length(izabrana_slova),min_primeraka);
for i = 1:length(izabrana_slova)
    slova_obelezja_lda(3*(i-1)+1:3*i,:) = A'*slova_obelezja(6*(i-1)+1:6*i,:) + 2;
end

f1min = inf;
f2min = inf;
f3min = inf;
f1max = -inf;
f2max = -inf;
f3max = -inf;

for i = 1:length(izabrana_slova)
    f1 =  slova_obelezja_lda(3*(i-1)+1,:);
    f2 =  slova_obelezja_lda(3*(i-1)+2,:);
    f3 =  slova_obelezja_lda(3*(i-1)+3,:);
    
    if f1min > min(f1)
        f1min = min(f1);
    end
    if f2min > min(f2)
        f2min = min(f2);
    end
    if f3min > min(f3)
        f3min = min(f3);
    end
    if f1max < max(f1)
        f1max = max(f1);
    end
    if f2max < max(f2)
        f2max = max(f2);
    end
    if f3max < max(f3)
        f3max = max(f3);
    end
    
    figure(length(izabrana_slova)+1)
    hold all;
    scatter3(f1,f2,f3);
    xlabel('lda1');
    ylabel('lda2');
    zlabel('lda3');
    title('Prostorni raspored odabiraka u redukovanom prostoru')
end

figure(length(izabrana_slova)+1)
legend('a','y','h','i','s','w','b','z','v','c')

 %% KNN procena fgv
N = length(izabrana_slova)*min_primeraka;
K = 70; % alfa = 0.6

x = f1min:1e-1:f1max;
y = f2min:1e-1:f2max;
z = f3min:1e-1:f3max;

fgvs = zeros(length(izabrana_slova),length(x),length(y),length(z));

f1 = zeros(1,min_primeraka);
f2 = zeros(1,min_primeraka);
f3 = zeros(1,min_primeraka);

for n = 1:length(izabrana_slova)
    n
    fgv = zeros(length(x),length(y),length(z));
    for i = 1:length(x)
        for j = 1:length(y)
            for k = 1:length(z)
                X = [x(i);y(j);z(k)];
                
                f1 = slova_obelezja_lda(3*(n-1)+1,:);
                f2 = slova_obelezja_lda(3*(n-1)+2,:);
                f3 = slova_obelezja_lda(3*(n-1)+3,:);
                
                d = zeros(1,length(f1));
                for m = 1:length(f1)
                    d(m) = (X - [f1(m);f2(m);f3(m)])'*(X - [f1(m);f2(m);f3(m)]);
                end
                d = sort(d);
                d = d(1:K);
                
                R = d(K);
                v = 4/3*R^3*pi;
                
                fgv(i,j,k) = (K-1)/N/v;
            end
        end
    end
    fgvs(n,:,:,:) = fgv;
end

%% Bayes-ov test vise hipoteza
pred = zeros(1,N);
true = [];

for i = 1:length(izabrana_slova)
    true = [true i*ones(1,min_primeraka)];
end

for a = 1:length(izabrana_slova)
    for b = 1:min_primeraka
        X = slova_obelezja_lda(3*(a-1)+1:3*a,b);
        
        i = ceil((X(1) - min(x))/(max(x)-min(x))*length(x))+1;
        j = ceil((X(2) - min(y))/(max(y)-min(y))*length(y))+1;
        k = ceil((X(3) - min(z))/(max(z)-min(z))*length(z))+1;
        
        if i > length(x)
            i = length(x);
        end
        if j > length(y)
            j = length(y);
        end
        if k > length(z)
            k = length(z);
        end
        
        [~,pred(1,(a-1)*min_primeraka+b)] = max(fgvs(:,i,j,k));
    end

end
%% Rezultati
C = confusionmat(true,pred);

s = 0;
for i = 1:length(izabrana_slova)
    s = s + C(i,i);
end

tacnost = s/N
 
figure(12)
confusionchart(C,{'a','y','h','i','s','w','b','z','v','c'});


