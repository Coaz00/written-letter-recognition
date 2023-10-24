clear;
close all;
clc;

slova = load('PO_slova.mat');

izabrana_slova = ['d' 's'];

min_primeraka = inf;
for i = 1:length(izabrana_slova)
    tmp = slova.(izabrana_slova(i));
    
    if min_primeraka > length(tmp)
        min_primeraka = length(tmp);
    end
end

slova_obelezja = zeros(2*length(izabrana_slova),min_primeraka);
for i = 1:length(izabrana_slova)
     tmp = slova.(izabrana_slova(i));
     
     f1 = zeros(1,length(tmp));
     f2 = zeros(1,length(tmp));
     
     for j = 1:length(tmp)
        slovo = cell2mat(tmp(j));
        
        f1(j) = max(diff(slovo(1,:))) - min(diff(slovo(1,:)));
        f2(j) = max(diff(slovo(2,:))) - min(diff(slovo(2,:)));
     end    
     
     figure(1)
     scatter(f1,f2)
     hold all;
     
     slova_obelezja(2*i-1,:) = f1(1:min_primeraka);
     slova_obelezja(2*i,:) = f2(1:min_primeraka);
end
 
 figure(1)
 legend('d','s')
 
 D_f = slova_obelezja(1:2,:);
 S_f = slova_obelezja(3:4,:);
 
 %% Klasifikator distance
 Xs = [D_f S_f];
 
 true = [ones(1,min_primeraka) 2*ones(1,min_primeraka)];
 pred = zeros(1,2*min_primeraka);
 
 M_D = mean(D_f,2);
 M_S = mean(S_f,2);
 
 for i = 1:length(pred)
     X = Xs(:,i);
     
     e_D = (X-M_D)'*(X-M_D);
     e_S = (X-M_S)'*(X-M_S);
     
     if e_D < e_S
         pred(i) = 1;
     else
         pred(i) = 2;
     end
 end
 
 C = confusionmat(true,pred);
 
 k = - (M_S(1) - M_D(1))/(M_S(2)-M_D(2));
 centar = 0.5*(M_S + M_D);
 n = centar(1) - k*centar(2);

 x = 0.25:0.01:1.5;
 figure(1)
 plot(x,k*x+n);
 xlabel('$max(a_x)-min(a_x)$','Interpreter','Latex')
 ylabel('$max(a_y)-min(a_y)$','Interpreter','Latex')
 title('Klasifikacija odabiraka','Interpreter','Latex')
 
 figure(2)
 confusionchart(C,{'d','s'});