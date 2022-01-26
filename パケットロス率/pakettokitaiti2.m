%収容数等を引数としてパケットロス率の期待値を出力する関数
function Le = pakettokitaiti2(param)
z=1744;%[bit]              %G.711 1パケット当たり218[Byte]=1744[bit]
param.B=20;%[bps]           %10M  1×10^8
%b0=100000;%[bps]     https://www.school.ctc-g.co.jp/columns/inst/inst50.html
c1=0.05; %緊急呼の上限呼損率
%c2=[1]; %被災地呼の上限呼損率
c2=0.5;
param.K=5;
%K_list=[0:1:8];
%alphak=1/52800; %11行目alphakの値の分母の数値を変更する。 α1.5倍
%alphak=1/44000; % α1.25倍
alphak=1/35200; % α1.0倍
%alphak=1/31680; % α0.9倍
%alphak=1/28160; % α0.8倍
betak=1/65000;
alphah=1/35200;
betah=1/65000;
alphag=1/35200;
betag=1/65000;
T=1600;
%N_list=[1:1:20];
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
param.thin=thin;
param.thout=thout;
param.Nk=0;
param.Nh=0;
param.Ng=0;
rho2d_list = [0.4];%[0.05:0.05:1];  %被災地呼0.05から１まで0.05間隔でトラヒック密度をプロット
rho1d_list = [0.5];     %[0.3 0.9]　　%緊急呼トラヒック密度
rho3d_list = [0.8];  %被災地外呼0.05から１まで0.05間隔でトラヒック密度をプロット

%Th2_list    = 20;  %閾値2を0から全帯域まで調べる
%Th1_list    = 20;  %閾値1を0から全帯域まで調べる

 Th2_list    = [0:param.B];  %閾値2を0から全帯域まで調べる
 Th1_list    = [0:param.B];  %閾値1を0から全帯域まで調べる
 %Th2_list    = [0:20];  %閾値2を0から全帯域まで調べる
 %Th1_list    = [0:20];  %閾値1を0から全帯域まで調べる

tic;

Le=0;


%mean()で配列の平均値
%ave=[mean(tmp1) mean(tmp2) mean(tmp3) mean(tmp4) mean(tmp5) mean(tmp6) mean(tmp7) mean(tmp8) mean(tmp9) mean(tmp10) mean(tmp11) mean(tmp12) mean(tmp13) mean(tmp14) mean(tmp15) mean(tmp16) mean(tmp17) mean(tmp18) mean(tmp19) mean(tmp20)];
%ave=[sum(tmp1) sum(tmp2) sum(tmp3) sum(tmp4) sum(tmp5) sum(tmp6) sum(tmp7) sum(tmp8) sum(tmp9) sum(tmp10) sum(tmp11) sum(tmp12) sum(tmp13) sum(tmp14) sum(tmp15) sum(tmp16) sum(tmp17) sum(tmp18) sum(tmp19) sum(tmp20)];
start_i =1;

%もし，エラーで止まった場合には，ここのコメントアウトをはずす
%load temp.mat;
%start_i = i;
H=1;

for i=start_i:numel(rho1d_list);  %緊急呼のトラヒック密度    
    for j=1:numel(rho2d_list);   %被災地呼のトラヒック密度        
        for k=1:numel(rho3d_list);  %被災地外呼のトラヒック密度
            for x=H:numel(c2);
                for l=1:numel(Th2_list);   %閾値2を1から全帯域まで  n = numel(A)は配列Aの要素数nを返す。 
                    for m=1:numel(Th1_list);   %閾値1を1から(閾値2)-1まで 
                        for n=1:numel(Th1_list)
                            param.rho1dash = rho1d_list(i);
                            param.rho2dash = rho2d_list(j);
                            param.rho3dash = rho3d_list(k);
                            
                            param.Nk=l;
                            param.Nh=m;
                            param.Ng=n;
                            stateR = sanjigenn2(param);     
                            %各呼損率を各変数に格納
                            Le=Le+stateR;

                        end                                               
                    end
                end                       
            end %ここまで         
        end       
    end
    % iのループの最後で毎回ワークスペース中の全変数を保存
    save temp.mat
    
end


toc;