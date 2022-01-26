%収容数毎のパケットロス率をグラフで表示するプログラム
clear; %関数clearは、現在のワークスペースからすべての変数を削除して、システムメモリから解放する。
close all; %ハンドルが表示されているすべての Figure を閉じる。

z=1744;%[bit]              %G.711 1パケット当たり218[Byte]=1744[bit]
B=100000000;%[bps]           %10M  1×10^8
%b0=100000;%[bps]     https://www.school.ctc-g.co.jp/columns/inst/inst50.html

param.K=5;
%K_list=[0:1:8];
alphak=1/52800; %11行目alphakの値の分母の数値を変更する。 α1.5倍
%alphak=1/44000; % α1.25倍
%alphak=1/35200; % α1.0倍
%alphak=1/31680; % α0.9倍
%alphak=1/28160; % α0.8倍
betak=1/65000;
alphah=1/35200;
betah=1/65000;
alphag=1/35200;
betag=1/65000;
T=1600;
N_list=[1:1:20];
%param.mu=(B/z)/1000000;%%%%[s]→[μs]のため/1000000  1×10^6
param.mu=0.043;%0.057339;%0.0299;
%param.mu=0.1;

lamdas0k=(betak)/(T*(alphak+betak));%param
ca2jouk=(1-(1-alphak*T)^2)/((T^2)*(alphak+betak)^2); %ok 18.0950
ca3jouk=ca2jouk^(3/2);
ca4jouk=ca2jouk^2;
sk0k=(2*alphak*T*(((alphak*T)^2)-3*alphak*T+3))/((alphak*T*(2-alphak*T))^(3/2));
Dk=(3*lamdas0k*(ca2jouk-1))/(2*sk0k*ca3jouk-3*ca4jouk-1);%param
Fk=Dk*(3*ca4jouk-sk0k*ca3jouk-3*ca2jouk+2)/(3*(ca2jouk-1));%param
Ek=Dk*(ca2jouk-1)/(Fk^2);%param
lamdas0h=(betah)/(T*(alphah+betah));
ca2jouh=(1-(1-alphah*T)^2)/((T^2)*(alphah+betah)^2); %ok 18.0950
ca3jouh=ca2jouh^(3/2);
ca4jouh=ca2jouh^2;
sk0h=(2*alphah*T*(((alphah*T)^2)-3*alphah*T+3))/((alphah*T*(2-alphah*T))^(3/2));
Dh=(3*lamdas0h*(ca2jouh-1))/(2*sk0h*ca3jouh-3*ca4jouh-1);
Fh=Dh*(3*ca4jouh-sk0h*ca3jouh-3*ca2jouh+2)/(3*(ca2jouh-1));
Eh=Dh*(ca2jouh-1)/(Fh^2);
lamdas0g=(betag)/(T*(alphag+betag));
ca2joug=(1-(1-alphag*T)^2)/((T^2)*(alphag+betag)^2); %ok 18.0950
ca3joug=ca2joug^(3/2);
ca4joug=ca2joug^2;
sk0g=(2*alphag*T*(((alphag*T)^2)-3*alphag*T+3))/((alphag*T*(2-alphag*T))^(3/2));
Dg=(3*lamdas0g*(ca2joug-1))/(2*sk0g*ca3joug-3*ca4joug-1);
Fg=Dg*(3*ca4joug-sk0g*ca3joug-3*ca2joug+2)/(3*(ca2joug-1));
Eg=Dg*(ca2joug-1)/(Fg^2);

param.lamdas0k=lamdas0k;
param.lamdas0h=lamdas0h;
param.lamdas0g=lamdas0g;
param.Dk=Dk;
param.Ek=Ek;
param.Fk=Fk;
param.Dh=Dh;
param.Eh=Eh;
param.Fh=Fh;
param.Dg=Dg;
param.Eg=Eg;
param.Fg=Fg;
thin=17;
thout=15;
param.Nk=0;
param.Nh=0;
param.Ng=0;


tic;

tmp1=[];
tmp2=[];
tmp3=[];
tmp4=[];
tmp5=[];
tmp6=[];
tmp7=[];
tmp8=[];
tmp9=[];
tmp10=[];
tmp11=[];
tmp12=[];
tmp13=[];
tmp14=[];
tmp15=[];
tmp16=[];
tmp17=[];
tmp18=[];
tmp19=[];
tmp20=[];
ave=[];
res=0;
cou=0;

for i=0:20              %緊急
    for j=0:thin        %被災地
        for k=0:thout   %被災地外
            N=i+j+k;
            if(N<=20)
                param.Nk=i;
                param.Nh=j;
                param.Ng=k;

                if(N<=20)
                    res=pakettofunc(param);
                    switch N
                        case 1
                            tmp1=[tmp1 res];
                        case 2
                            tmp2=[tmp2 res];
                        case 3
                            tmp3=[tmp3 res];
                        case 4
                            tmp4=[tmp4 res];
                        case 5
                            tmp5=[tmp5 res];
                        case 6
                            tmp6=[tmp6 res];
                        case 7
                            tmp7=[tmp7 res];
                        case 8 
                            tmp8=[tmp8 res];
                        case 9
                            tmp9=[tmp9 res];
                        case 10
                            tmp10=[tmp10 res];
                        case 11
                            tmp11=[tmp11 res];
                        case 12
                            tmp12=[tmp12 res];
                        case 13
                            tmp13=[tmp13 res];
                        case 14
                            tmp14=[tmp14 res];
                        case 15
                            tmp15=[tmp15 res];
                        case 16
                            tmp16=[tmp16 res];
                        case 17
                            tmp17=[tmp17 res];
                        case 18
                            tmp18=[tmp18 res];
                        case 19
                            tmp19=[tmp19 res];
                        case 20
                            tmp20=[tmp20 res];
                    end
                end
            end
        end
    end
end
%mean()で配列の平均値
ave=[mean(tmp1) mean(tmp2) mean(tmp3) mean(tmp4) mean(tmp5) mean(tmp6) mean(tmp7) mean(tmp8) mean(tmp9) mean(tmp10) mean(tmp11) mean(tmp12) mean(tmp13) mean(tmp14) mean(tmp15) mean(tmp16) mean(tmp17) mean(tmp18) mean(tmp19) mean(tmp20)];

toc;

figure
semilogy(N_list, ave, 'bo-');
yline(0.001);
xlabel('Number of accommodated','FontSize',23);
ylabel('Packet loss rate','FontSize',23);



