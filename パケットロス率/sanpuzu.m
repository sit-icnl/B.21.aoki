%収容数を指定して、その際にパケットロス率が上限を超えるかを３次元プロットで示すプログラム

clear; %関数clearは、現在のワークスペースからすべての変数を削除して、システムメモリから解放する。
close all; %ハンドルが表示されているすべての Figure を閉じる。

z=1744;%[bit]              %G.711 1パケット当たり218[Byte]=1744[bit]
B=100000000;%[bps]           %10M  1×10^8
%b0=100000;%[bps]     https://www.school.ctc-g.co.jp/columns/inst/inst50.html

param.K=5;
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
param.mu=0.043;%0.0299;%0.057339;

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
restmp=[];
ave=[];
num=0;
cou=0;
itmp=[];
jtmp=[];
ktmp=[];
high=[];
low=[];
ftmp=[];
flg=0;

for i=0:20              %緊急
    for j=0:20          %被災地
        for k=0:20      %被災地外
            N=i+j+k;
            if(N<=20)
                param.Nk=i;
                param.Nh=j;
                param.Ng=k;
            
                if(N==17) %87行目のif文のなかの数値を変更する。
                    restmp=[restmp pakettofunc(param)];
                    itmp=[itmp i];
                    jtmp=[jtmp j];
                    ktmp=[ktmp k];
                    ftmp=[ftmp flg];
                    cou=cou+1;
                end
                flg=0;
            end
        end
    end
end

%
count=0;
for i=1:cou
    if(ftmp(i)==0)
        if(restmp(i) > 0.001)
            scatter3(itmp(i),jtmp(i),ktmp(i),36,[1 0 0]);%%%%%%%%%%　　赤で,'o','filled'
            hold on
            count=count+1;
        else
            scatter3(itmp(i),jtmp(i),ktmp(i),36,[0 0 1]);%%%%%%%%%%　　青で
            hold on
            count=count+1;
        end
    else
        if(restmp(i) > 0.001)
            scatter3(itmp(i),jtmp(i),ktmp(i),36,[0 0 0]);%%%%%%%%%%　　黒
            hold on
            count=count+1;
        else
            scatter3(itmp(i),jtmp(i),ktmp(i),36,[0 1 0]);%%%%%%%%%%　　緑
            hold on
            count=count+1;
        end
    end
end
hold off
%xlabel('緊急呼収容本数','FontSize',23);
%ylabel('被災地呼収容本数','FontSize',23);
%zlabel('被災地外呼収容本数','FontSize',23);
xlabel('emergency call','FontSize',23);
ylabel('general in call','FontSize',23);
zlabel('general out call','FontSize',23);
toc;