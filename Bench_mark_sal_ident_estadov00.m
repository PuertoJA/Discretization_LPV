%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Estudo de LPV          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Topico: LFT up para conversao de freq
%%%% Autor:  Jorge Andres Puerto Acosta
%%%% Data :  marco 2016
%%%% Versao: 0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%% Definicao %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1: fu(G,1/s)= M22 +M21(I-M11(1/s))^(-1)M12 = D + C(1/s)(I-A(1/s))^(-1)B 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sistema massa mola amortecedor com incerteza  LFT em K,c,m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms cb mb kb  delta_m  ym s

M=[ 0    ,     1,        0,           0;
   -kb/mb,-cb/mb, -inv(mb),-0.1*inv(mb);
    1    ,     0,        0,           0;
    0    ,      1,       0,           0;
      -kb,   -cb,        1,        -0.1];

delta=[delta_m];

M11=M(1:4,1:3);
M12=M(1:4,4:4);
M21=M(5:5,1:3);
M22=M(5,4);
I1=eye(size(M22));
fl= M11+ M12*delta*inv(I1 - M22*delta)*M21; 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sistema massa mola amortecedor com incerteza  LFT em K,c,m continuo (1/s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta1=inv(s);
fl11=fl(1:2,1:2);
fl12=fl(1:2,3);
fl21=fl(3:4,1:2);
fl22=fl(3:4,3);
Ifl=eye(size(fl11));
fu= fl22 + fl21*delta1*inv(Ifl - fl11*delta1)*fl12;        
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Resposta da LFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m1=1.4*10^(-5);
m2=1.0*10^(-5);
m3=4*1.34*10^(-3);
mb=m1+m2+m3;
cb=6.25*10^(-3);
kb= 1;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% teste delta vetor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tmu= 1.05e-2;
t=0:Tmu:0.105;
u=(1 + 0.25.*randn(1,length(t)));%N(1,0.25)
% u=(1 + 0.5.*randn(1,length(t)));%N(1,0.5)
delta_masa=(0 + 0.01.*randn(1,length(t)));%N(0,0.01)
mag=[];
fase=[];
%freq=[];
freq=logspace(0.1,1.5051,100);
% % freq=logspace(0.1,3,100);

wtf=freq;
for l=1:length(delta_masa)
      delta_m=delta_masa(l);
% % % fuev=eval(fu);
% % % [symNum,symDen] = numden(fuev); %obtem num e den da tf simbolica
% % % TFnum1 = sym2poly(symNum(1,1));    %comverte simbolico a polinomio
% % % TFden1 = sym2poly(symDen(1,1));    %comverte simbolico a polinomio
% % % TFnum2 = sym2poly(symNum(2,1));    %comverte simbolico a polinomio
% % % TFden2 = sym2poly(symDen(2,1));    %comverte simbolico a polinomio
% % % nums = {TFnum1;TFnum2};
% % % dens={TFden1;TFden2};
% % % HZ =tf(nums,dens);
% % % Hss=ss(HZ);
% bode(HZ,'r',Hss,'+k');
% [magtf,phasetf,wtf]= bode(HZ);
HZ=ss(eval(fl11),eval(fl12),eval(fl21),eval(fl22));
[magtf,phasetf]= bode(HZ,wtf);
magtf1=magtf(1,:,:);
magtf2=magtf(2,:,:);
phasetf1=phasetf(1,:,:);
phasetf2=phasetf(2,:,:);
eval([ 'mag1' num2str(l) ' = magtf1;' ]);
eval([ 'fase1' num2str(l) ' = phasetf1;' ]);
eval([ 'mag2' num2str(l) ' = magtf2;' ]);
eval([ 'fase2' num2str(l) ' = phasetf2;' ]);
% eval([ 'freq' num2str(l) ' = wtf;' ]);
mag1(:,l)=magtf1;
mag2(:,l)=magtf2;
fase1(:,l)=phasetf1;
fase2(:,l)=phasetf2;

oo=length(wtf);
aa=delta_m.*ones(1,oo);
% % %     
% % %     figure(1)
% % %     hold on
% % %     subplot(2,1,1),plot3(wtf,aa,20*log10(abs(magtf1(:))));
% % %     xlabel('Frequencia [rad]'),ylabel('delta m'),zlabel('Magnitude')
% % %     grid
% % %     hold off
% % %     hold on
% % %     subplot(2,1,2),plot3(wtf,aa,phasetf1(:));
% % %     xlabel('Frequencia [rad]'),ylabel('delta m'),zlabel('Fase')
% % %     grid
% % %     hold off
% % %     
% % %     
% % %     figure(2)
% % %     hold on
% % %     subplot(2,1,1),plot3(wtf,aa,20*log10(abs(magtf2(:))));
% % %     xlabel('Frequencia [rad]'),ylabel('delta m'),zlabel('Magnitude')
% % %     grid
% % %     hold off
% % %     hold on
% % %     subplot(2,1,2),plot3(wtf,aa,phasetf2(:));
% % %     xlabel('Frequencia [rad]'),ylabel('delta m'),zlabel('Fase')
% % %     grid
% % %     hold off
end

figure;
[X,Y] = meshgrid(wtf,delta_masa);
subplot(2,1,1),mesh(mag1');
set(gca,'XScale','log'),zlabel('Magnitude'),xlabel('freq [rad]'),ylabel('delta masa')
subplot(2,1,2),mesh(fase1');
set(gca,'XScale','log'),zlabel('Fase'),xlabel('freq [rad]'),ylabel('delta masa')
title('Saida 1')
% print -depsc bode1_masa_mola_amortecedor.eps

figure;
[X,Y] = meshgrid(wtf,delta_masa);
subplot(2,1,1),mesh(mag2');
set(gca,'XScale','log'),zlabel('Magnitude'),xlabel('freq [rad]'),ylabel('delta masa')
subplot(2,1,2),mesh(fase2');
set(gca,'XScale','log'),zlabel('Fase'),xlabel('freq [rad]'),ylabel('delta masa')
title('Saida 2')
% print -depsc bode2_masa_mola_amortecedor.eps
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%sis discreto %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
II=eye(size(fl11));
At=[-inv(II-fl11)*(II+fl11)];
Bt=[-(sqrt(2))*inv(II-fl11)*fl12];
Ct=[sqrt(2)*fl21*inv(II-fl11)];
Dt=[fl21*inv(II-fl11)*fl12+fl22];
M_tilde=[At,Bt;Ct,Dt]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for l=1:length(delta_masa)
    delta_m=delta_masa(l);
% % fuevd=eval(fudis);
% % [symNumd,symDend] = numden(fuevd); %obtem num e den da tf simbolica
% % TFnum1d = sym2poly(symNumd(1,1));    %comverte simbolico a polinomio
% % TFden1d = sym2poly(symDend(1,1));    %comverte simbolico a polinomio
% % TFnum2d = sym2poly(symNumd(2,1));    %comverte simbolico a polinomio
% % TFden2d = sym2poly(symDend(2,1));    %comverte simbolico a polinomio
% % numsd = {TFnum1d;TFnum2d};
% % densd={TFden1d;TFden2d};
% % HZd =tf(numsd,densd,Tmu);
% % Hssd=ss(HZd);
% bode(HZ,'r',Hss,'+k');
% [magtf,phasetf,wtf]= bode(HZ);
HZd=ss(eval(At),eval(Bt),eval(Ct),eval(Dt),Tmu);
[magtfd,phasetfd]= bode(HZd,wtf);
magtf1d=magtfd(1,:,:);
magtf2d=magtfd(2,:,:);
phasetf1d=phasetfd(1,:,:);
phasetf2d=phasetfd(2,:,:);
eval([ 'mag1d' num2str(l) ' = magtf1d;' ]);
eval([ 'fase1d' num2str(l) ' = phasetf1d;' ]);
eval([ 'mag2d' num2str(l) ' = magtf2d;' ]);
eval([ 'fase2d' num2str(l) ' = phasetf2d;' ]);
% eval([ 'freq' num2str(l) ' = wtf;' ]);
mag1d(:,l)=magtf1d;
mag2d(:,l)=magtf2d;
fase1d(:,l)=phasetf1d;
fase2d(:,l)=phasetf2d;

oo=length(wtf);
aa=delta_m.*ones(1,oo);
% % %     
% % %     figure(5)
% % %     hold on
% % %     subplot(2,1,1),plot3(wtf,aa,20*log10(abs(magtf1d(:))));
% % %     xlabel('Frequencia [rad]'),ylabel('delta m'),zlabel('Magnitude')
% % %     grid
% % %     hold off
% % %     hold on
% % %     subplot(2,1,2),plot3(wtf,aa,phasetf1d(:));
% % %     xlabel('Frequencia [rad]'),ylabel('delta m'),zlabel('Fase')
% % %     grid
% % %     hold off
% % %     
% % %     
% % %     figure(6)
% % %     hold on
% % %     subplot(2,1,1),plot3(wtf,aa,20*log10(abs(magtf2d(:))));
% % %     xlabel('Frequencia [rad]'),ylabel('delta m'),zlabel('Magnitude')
% % %     grid
% % %     hold off
% % %     hold on
% % %     subplot(2,1,2),plot3(wtf,aa,phasetf2d(:));
% % %     xlabel('Frequencia [rad]'),ylabel('delta m'),zlabel('Fase')
% % %     grid
% % %     hold off
end

figure;
[X,Y] = meshgrid(wtf,delta_masa);
subplot(2,1,1),mesh(mag1d');
set(gca,'XScale','log'),zlabel('Magnitude'),xlabel('freq [rad]'),ylabel('delta masa')
subplot(2,1,2),mesh(fase1d');
set(gca,'XScale','log'),zlabel('Fase'),xlabel('freq [rad]'),ylabel('delta masa')
title('Saida 1 dis')
% print -depsc bode1d_masa_mola_amortecedor.eps

figure;
[X,Y] = meshgrid(wtf,delta_masa);
subplot(2,1,1),mesh(mag2d');
set(gca,'XScale','log'),zlabel('Magnitude'),xlabel('freq [rad]'),ylabel('delta masa')
subplot(2,1,2),mesh(fase2d');
set(gca,'XScale','log'),zlabel('Fase'),xlabel('freq [rad]'),ylabel('delta masa')
title('Saida 2 dis')
% print -depsc bode2d_masa_mola_amortecedor.eps
figure,bode(HZ,'r+',HZd),legend('CC','DT'),title('Comparaçao de bode modelo local')
print -depsc bode_2xmasa_mola_amortecedor.eps
% break
%% %%%%%%%%%%%%%%% Resposta ao impulso e ao degrau %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[y_cc_imp]=impulse(HZ,t);
[y_cc_degrau]=step(HZ,t);
%% %%%%%%%%%%%%%%% Resposta forcada %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y_cc=lsim(Hss,u,t);
%% %%%%%%%%%%%%%%% Sistema discreto %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%% Resposta ao impulso e ao degrau %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[y_dis_imp]=impulse(HZd,t);
[y_dis_degrau]=step(HZd,t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(t,y_cc_degrau,'bo',t,y_dis_degrau,'k')
legend('Continuo','Discreto')
grid

figure
plot(t,y_cc_imp,'b+',t,y_dis_imp,'k'),legend('Continuo','Discreto')
grid

