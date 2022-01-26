%三元シミュレーション シミュレーション用のプログラム。各呼のトラヒック密度の条件を設定すれば結果が出る。しかし、呼の到着部分に関しては、各呼が独立して到着するように修正をしなければならない。
deltatime=0.01;  %最小時間[秒/タイムスロット]

%%



tic;
% A=[];
% B=[];
% C=[];
% I=[];
% J=[];
% K=[];

%%

    %rho1 =[0.05 0.1];
    rho1 =[0.5];%緊急呼トラヒック密度
    rho2 =[0.05:0.05:1];%被災地呼トラヒック密度
    %rho2 =[0.05];%被災地外呼トラヒック密度
    rho3 =[0.8];
D=[950	850	800	750	700	700	700	700	700	700	650	650	650	650	650	650	650	650	650	650];%被災地呼閾値
 E=[600	550	600	600	600	600	650	650	650	650	600	600	600	600	600	600	600	600	600	600];%被災地外呼閾値

%     D=1000;
%     E=1000;
%     O=[];
%     P=[];
%     Q=[];
    for p=1:numel(rho2)%閾値とトラヒック密度回す用
    D_a=D(p);
    E_a=E(p);
  rho2_a=rho2(p);
    for q=1:10%試行回数
        
        

        %配列などリセット
        a=0;
        lambda1=0;
        lambda2=0;
        lambda3=0;
        e=0;
        gin=0;
        gout=0;
        kaisen=0;
        e_t=[];
        gin_t=[];
        gout_t=[];
        e_loss=0;
        gin_loss=0;
        gout_loss=0;
        e_ok=0;
        gin_ok=0;
        gout_ok=0;
        kaisen1=0;
        kaisen2=0;
        kaisen3=0;
        kaisen4=0;
        kaisen5=0;
        kaisen6=0;
        
        for o=1:1000000 %試行回数
            a=a+1;
            if(a==300000)%最初のほうのデータをリセット
                
                e=0;
                gin=0;
                gout=0;
                e_loss=0;
                gin_loss=0;
                gout_loss=0;
                e_ok=0;
                gin_ok=0;
                gout_ok=0;
            end
            lambda1=1000*rho1*exprnd(0.01);%到着率を導出
            e_l=lambda1;%到着率を記録
            lambda2=1000*rho2_a*exprnd(0.01);
            gin_l=lambda2;
            lambda3=1000*rho3*exprnd(0.01);
            gout_l=lambda3;
            
            if(poissrnd(0.01*lambda1)~=0)
                e=e+1;%到着した緊急呼を加算
                if(kaisen<1000)
                    e_ok=e_ok+1;%格納した緊急呼を加算
                    kaisen4=kaisen+1;
                    e_t(e_ok)=(1000*rho1)/lambda1;%処理時間を計算
                else
                    e_loss=e_loss+1;%呼損した緊急呼を加算
                end
            end
            if(poissrnd(0.01*lambda2)~=0)
                gin=gin+1;
                %% 
                if(kaisen4<D_a)
                    gin_ok=gin_ok+1;
                    kaisen5=kaisen4+1;
                    gin_t(gin_ok)=(1000*rho2_a)/lambda2;
                else
                    gin_loss=gin_loss+1;
                end
            end
            
            if(poissrnd(0.01*lambda3)~=0)
                gout=gout+1;
                if(kaisen5<E_a)
                    gout_ok=gout_ok+1;
                  %  kaisen6=kaisen5+1;
                    gout_t(gout_ok)=(1000*rho3)/lambda3;
                else
                    gout_loss=gout_loss+1;
                end
            end
            %kaisen=0;
            
            
            
            for t=1:e_ok
                e_t(t)=e_t(t)-0.01;%処理時間から経過時間を引く
                if(e_t(t)>0)
                    kaisen1=kaisen1+1;%処理が終わってない呼を数える
                    end
            end
            for t=1:gin_ok
                gin_t(t)=gin_t(t)-0.01;
                if(gin_t(t)>0)
                    kaisen2=kaisen2+1;
                end
            end
            for t=1:gout_ok
                gout_t(t)=gout_t(t)-0.01;
                if(gout_t(t)>0)
                    kaisen3=kaisen3+1;
                end
            end
            kaisen=0;
            kaisen=kaisen1+kaisen2+kaisen3;
            kaisen1=0;
            kaisen2=0;
            kaisen3=0;
            kaisen4=kaisen;
            kaisen5=kaisen;
    %        kaisen6=kaisen;
        end
        
        
        
        A(q)=e_loss/e;%呼損率計算
        B(q)=gin_loss/gin;
        C(q)=gout_loss/gout;
        
%         I(p)=e;
%         J(p)=gin;
%         K(p)=gout;
    end
    
%     L=sum(I)/1;
%     M=sum(J)/1;
%     N=sum(K)/1;
    
    F=sum(A)/10;%試行回数で割る
    G=sum(B)/10;
    H=sum(C)/10;
    
    O(p)=F;%各トラヒック条件の箱に格納
    P(p)=G;
    Q(p)=H;
    
    fprintf('%0.3f  %0.3f  %0.3f \n',  F,  G,  H);
    
    end



%%
figure(1);
plot(rho2,F,'ro-',rho2,G,'b^-',rho2,H,'g*-');
set(gca,'FontSize',14);
xlabel('Traffic intensity of emergency call,rho2') % x-axis label
ylabel('Call blocking probability') % y-axis label
legend('Emergency call','Disaster general call','Non-disaster general call')

figure(2);
plot(rho2,L,'ro-',rho2,M,'b^-',rho2,N,'g*-');
set(gca,'FontSize',14);
xlabel('Traffic intensity of emergency call,rho2') % x-axis label
ylabel('Call blocking probability') % y-axis label
legend('Emergency call','Disaster general call','Non-disaster general call')


% %通常コンテンツと簡易コンテンツの平均数のグラフを追加(呼損率の確認及び考察用)
% figure(2);
% plot(linkspeed,B,'ro-')%,linkspeed,B,'b^-');
% set(gca,'FontSize',14);
% xlabel('Linkspeed,L') % x-axis label
% ylabel('通常コンテンツの個数(平均)') % y-axis label
% legend('Conventional')%,'Proposed')
% 
% figure(3);
% plot(linkspeed,C,'ro-')%,linkspeed,B,'b^-');
% set(gca,'FontSize',14);
% xlabel('Linkspeed,L') % x-axis label
% ylabel('簡易コンテンツの個数(平均)') % y-axis label
% legend('Conventional')%,'Proposed')

toc;


