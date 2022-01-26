% 閾値を総当たりするプログラム。また、各呼のトラヒック密度や上限呼損率、退去率等を設定している。現在は何の縛りもなくすべての状態における各呼の呼損率を表示している。もし、ある一定の呼損率以下の場合のみ実数値を格納したい場合は、100～122行目のコメントマークを解除し、自分で各種条件を設定する。
% 2元受付制御(遠慮)トータル呼損率比較グラフおよび最適閾値グラフ
% 2017.07.24
% by Sumiko Miyata

%3元受付制御　グラフ出力
%2018.09.19
%川瀬

% 遠慮行動.による退去を実装するために
% deltaとmu3を追加

clear;
close all;

param.B    = 20;    %全帯域
param.b1   = 1;      %緊急呼占有帯域
param.b2   = 1;      %被災地呼占有帯域
param.b3   = 1;      %被災地外呼占有帯域

param.mu1  = 0.01; %緊急呼退去率
param.mu2  = 0.01; %被災地呼退去率
param.mu3  = 0.01; %被災地外呼退去率
ok=0;

c1=0.05; %緊急呼の上限呼損率
%c2=[1]; %被災地呼の上限呼損率
c2=0.5;

rho2d_list = [0.4];%[0.05:0.05:1];  %被災地呼0.05から１まで0.05間隔でトラヒック密度をプロット
rho1d_list = [0.5];     %[0.3 0.9]　　%緊急呼トラヒック密度
rho3d_list = [0.8];  %被災地外呼0.05から１まで0.05間隔でトラヒック密度をプロット

%Th2_list    = 20;  %閾値2を0から全帯域まで調べる
%Th1_list    = 20;  %閾値1を0から全帯域まで調べる

 Th2_list    = [0:param.B];  %閾値2を0から全帯域まで調べる
 Th1_list    = [0:param.B];  %閾値1を0から全帯域まで調べる
 %Th2_list    = [0:20];  %閾値2を0から全帯域まで調べる
 %Th1_list    = [0:20];  %閾値1を0から全帯域まで調べる
H=1;

tic; %時間測定

% for i=1:numel(rho1d_list);
%     for j=1:numel(rho2d_list);
%           for l=1:numel(rho3d_list);
%         param.rho1dash = rho1d_list(i);
%         param.rho2dash = rho2d_list(j);
%         param.rho3dash = rho3d_list(l);
%         res_S(i,j,l) = func_CAC_simple(param);
%           end
%     end
% end

%原因は不明だが，固有値計算の途中でエラーになることがあるため，もし，エラーで止まった場合には，下記をコメントアウトして，
start_i =1;

%もし，エラーで止まった場合には，ここのコメントアウトをはずす
%load temp.mat;
%start_i = i;

for i=start_i:numel(rho1d_list);  %緊急呼のトラヒック密度
    
    for j=1:numel(rho2d_list);   %被災地呼のトラヒック密度
        
        for l=1:numel(rho3d_list);  %被災地外呼のトラヒック密度
            
            for x=H:numel(c2);
                c3=c2(H);
                ok=0;
                th_gin=[];
                th_gout=[];
                e_call_late=[];
                gin_call_late=[];
                gout_call_late=[];
                a=1;
                b=1;
                for k=1:numel(Th2_list);   %閾値2を1から全帯域まで  n = numel(A)は配列Aの要素数nを返す。 
                     a=k;
                    for m=1:numel(Th1_list);   %閾値1を1から(閾値2)-1まで
                     b=m;   
                        param.rho1dash = rho1d_list(i);
                        param.rho2dash = rho2d_list(j);
                        param.rho3dash = rho3d_list(l);
                        param.B_th2     = Th2_list(k);
                        param.B_th1     = Th1_list(m);
                        [res_temp_e(k,m), res_temp_gin(k,m),res_temp_gout(k,m)] = sanjigenn(param);
                        %90行目：各呼に重みづけをした評価を行おうとしたがやめた。
                          %res_temp_sum(k,m)=0.5*res_temp_e(k,m)+0.4*res_temp_gin(k,m)+0.8*res_temp_gout(k,m);
                        
                        %各呼損率を各変数に格納
                        e_call_rate(k,m)=res_temp_e(k,m);
                        gin_call_rate(k,m)=res_temp_gin(k,m);
                        gout_call_rate(k,m)=res_temp_gout(k,m);
                        
                        
                        
                        %100～122行目：被災地呼の呼損率が被災地外呼より小さく、被災地呼の呼損率がc２以下であり、緊急呼の呼損率がc１以下の時は各数値を変数に格納。そうでないときは０にする。
                         if(res_temp_gin(k,m)<res_temp_gout(k,m)&&res_temp_e(k,m)<=c1&&res_temp_gin(k,m)<=c3)
                             
                             gout_call_rate(k,m)=res_temp_gout(k,m);
                             ok=ok+1;
                             e_call_late(a,b)=res_temp_e(k,m);
                             gin_call_late(a,b)=res_temp_gin(k,m);
                             gout_call_late(a,b)=res_temp_gout(k,m);
                             th_gin(a)=k-1;
                             th_gout(b)=m-1;
 
                             
                             
                             
                             
                         else
                             
                             gout_call_rate(k,m) =0;
                             gout_call_late(a,b)=1;                             
                         end
                        
                        
                        fprintf('e_call_rate=%0.3f  gin_call_rate=%0.3f  gout_call_rate=%0.3f  Th_gin=%0.3f  Th_gout=%0.3f ok=%0.3f  \n',  e_call_rate(k,m),  gin_call_rate(k,m),  gout_call_rate(k,m), k, m , ok);
                   
                        fprintf('%0.3f  %0.3f  %0.3f \n',  e_call_rate(k,m),  gin_call_rate(k,m),  gout_call_rate(k,m));
                    end
                end
                
                %128～140行目：条件を満たす閾値の数が一定の状況において、被災地外呼の呼損率が最小になる閾値及び各呼の呼損率を導出
                 if(ok>44)
                     Min=min(gout_call_late,[],'all');
                     I=a*b;                    
                     [C,I]=min(gout_call_late(:));
                     [ii,jj]=ind2sub(size(gout_call_late),I);
                     
                     fprintf('e_call_rate=%0.3f  gin_call_rate=%0.3f  gout_call_rate=%0.3f  Th_gin=%0.3f  Th_gout=%0.3f ok=%0.3f c2=%0.3f \n',  e_call_late(ii,jj),  gin_call_late(ii,jj),  gout_call_late(ii,jj), th_gin(ii), th_gout(jj) , ok, c3);
                     fprintf('%0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f \n',  e_call_late(ii,jj),  gin_call_late(ii,jj),  gout_call_late(ii,jj), th_gin(ii), th_gout(jj) , ok, c3);
                     break;                     
                 else
                     H=H+1;
                 end
                
                
                %ここから
                %     if(res_temp_e(k,m)<c1&&res_temp_gin(k,m)<c2)           %緊急呼がある一定の呼損率以下なら
                %    res_gout(k,m) =res_temp_gout(k,m);  %ある一定の呼損率以下の緊急呼の呼損率の同じ閾値での一般呼の呼損率
                %    else
                %      res_gout(k,m)=1;           %それ以外は反映されないように
                %  end
                %  if(res_e(k,m)<c)           %緊急呼がある一定の呼損率以下なら
                % gout_call_rate(k,m) =res_temp_gout(k,m);  %ある一定の呼損率以下の緊急呼の呼損率の同じ閾値での一般呼の呼損率
                % else
                % gout_call_rate(k,m) =1;           %それ以外は反映されないように
                % end
                
                % emergency_r2(k,m,j,l,i) =res_temp_e(k,m);
                % generallin_r2(k,m,j,l,i) =res_temp_gin(k,m);
                % generallout_r2(k,m,j,l,i) =res_temp_gout(k,m);
                
                
                %   fprintf('rho1=%0.3f  rho2=%0.3f  rho3=%0.3f   \n', param.rho1dash, param.rho2dash, param.rho3dash);
                %    [res_min, min_ind] = min(gout_call_rate);    %被災地外呼の呼損率の最小値
                
                %         res_T(i,l,j)  = res_min;
                % r_opt(i,l,j)= res_min;  %最適閾値での呼損率
                %     Th2_opt(:,j,l,i) = Th2_list(min_ind);   %最適閾値2
                %     Th1_opt(:,j,l,i) = Th1_list(min_ind);   %最適閾値1
                %      th_opt(j)= Th_list(min_ind);
                %  emergency_r(:,j,l,i) =res_temp_e(min_ind);  %最適閾値での緊急呼の呼損率
                %    r_1(j) = res_temp_n(min_ind);
                %    generallin_r(:,j,l,i) = res_temp_gin(min_ind);   %最適閾値での被災地呼の呼損率
                %    r_1(j) = res_temp_n(min_ind);
                %  generallout_r(:,j,l,i) = res_temp_gout(min_ind);   %最適閾値での被災地外呼の呼損率
            end %ここまで
            
            
        end
        
        
    end
    % iのループの最後で毎回ワークスペース中の全変数を保存
    save temp.mat
    
end

%%

%for i=1:numel(rho1d_list);
%     rho1filename = round(param.rho1dash*100);
%     roptfile = sprintf('ropt%d.mat', rho1filename);
%     load(roptfile);
%end

toc;


 
 %最適緊急呼呼損率特性
%ylim([0 0.5])
%xlabel('Traffic intencity \rho_2','FontSize',18)             %横軸
%ylabel('emergency call blocking rate','FontSize',18)          %縦軸
 %legend({'m=0','m=1','m=2'},'FontSize',15)      %右上の説明欄
 

 %閾値を変化させた時の緊急呼呼損率特性
figure
 x=Th1_list;
 y=Th2_list; 
 z=e_call_rate;
 stem3(x,y,z,'r','o','LineStyle','none') % 3次元離散データ列のプロット マーカーは赤い円,ラインスタイルはライン無し
 zlim([0 1.0])
 xlabel('th gout','FontSize',18)
 ylabel('th gin','FontSize',18)
 zlabel('emergency call blocking late,r_e','FontSize',18)
 
 %閾値を変化させた時の被災地呼呼損率特性
 figure
 x=Th1_list;
 y=Th2_list; 
 z=gin_call_rate;
 stem3(x,y,z,'r','o','LineStyle','none') % 3次元離散データ列のプロット マーカーは赤い円,ラインスタイルはライン無し
 zlim([0 1.0])
 xlabel('th gout','FontSize',18)
 ylabel('th gin','FontSize',18)
 zlabel('Generalin call blocking late,r_{gin}','FontSize',18)

 %閾値を変化させた時の緊急呼呼損率特性
figure
 x=Th1_list;
 y=Th2_list; 
 z=gout_call_rate;
 stem3(x,y,z,'r','o','LineStyle','none') % 3次元離散データ列のプロット マーカーは赤い円,ラインスタイルはライン無し
 zlim([0 1.0])
 xlabel('th gout','FontSize',18)
 ylabel('th gin','FontSize',18)
 zlabel('Generalout call blocking late,r_{gout}','FontSize',18)

 %閾値を変化させた時の全ての呼の呼損率特性
 %figure
 %x=Th1_list;
 %y=Th2_list;
 %stem3(x,y,gout_call_rate,'b','d','LineStyle','none')
 %hold on
 %stem3(x,y,gin_call_rate,'k','s','LineStyle','none')
 %hold on
 %stem3(x,y,e_call_rate,'r','o','LineStyle','none')
 %hold on
 %stem3(x,y,sum_call_rate,'r','o','LineStyle','none')
 %xlabel('threshold th1','FontSize',30)
 %ylabel('threshold th2','FontSize',30)
%zlabel('call blocking late','FontSize',30)
%legend({'被災地外呼','被災地呼','緊急呼'},'FontSize',15)


%figure(1);
%plot(rho1d_list, res_T(:,1), 'ko', rho1d_list, res_T(:,2), 'k^', rho1d_list, res_T_c(:,1), 'k+', rho1d_list, res_T_c(:,2), 'k*', rho1d_list, res_T_o(:,1), 'kx', rho1d_list, res_T_o(:,2), 'ks','MarkerSize',8);
%ylim([0 0.3])
%xlabel('threshold th2','FontSize',30)
%ylabel('call blocking late','FontSize',30)
%legend('r^{\delta=0}_{opt}，\rho_2%=0.3','r^{\delta=0}_{opt}，\rho_2’=0.9','r^{prop}_{opt}，\rho_2’=0.3','r^{prop}_{opt}，\rho_2’=0.9','r^{con}_{opt}，\rho_2’=0.3','r^{con}_{opt}，\rho_2’=0.9','Location','northwest')
%title('トータル呼損率特性')




 
 %最適被災地呼呼損率特性
%  figure(2);
%plot(rho2d_list, generall_r(:,1,1), 'rp', rho2d_list,  generall_r(:,2,1), 'go', rho2d_list,  generall_r(:,3,1), 'b*','MarkerSize',15);
%ylim([0 1.0])
%xlabel('Traffic intencity \rho_2','FontSize',18)  
%ylabel('generall call blocking rate','FontSize',18)
%legend({'m=0','m=1','m=2'},'FontSize',15)
 
%最適被災地外呼呼損率特性
%  figure(2);
%plot(rho2d_list, generall_r(:,1,1), 'rp', rho2d_list,  generall_r(:,2,1), 'go', rho2d_list,  generall_r(:,3,1), 'b*','MarkerSize',15);
%ylim([0 1.0])
%xlabel('Traffic intencity \rho_2','FontSize',18)  
%ylabel('generall call blocking rate','FontSize',18)
%legend({'m=0','m=1','m=2'},'FontSize',15)


 %最適閾値呼損率特性
 %figure(3);
 %plot(rho2d_list,Th_opt(:,1,1), 'rp',rho2d_list,Th_opt(:,2,1), 'go',rho2d_list,Th_opt(:,3,1), 'b*','MarkerSize',15);
 %xlabel('Traffic intencity \rho_2','FontSize',18)
%ylabel('optimal threshold th_{opt}','FontSize',18)
% legend({'m=0','m=1','m=2'},'FontSize',15)
 %plot(Th1_list, gout_call_rate(:,:) , 'rp', Th1_list, gout_call_rate(2,:) , 'go', Th1_list,  gout_call_rate(3,:) , 'b*',Th1_list, gout_call_rate(4,:) , 'rx',Th1_list, gout_call_rate(5,:) , 'g>',Th1_list, gout_call_rate(6,:) , 'b<','MarkerSize',11);

 
 %閾値変化させたとき
% figure(4);
%plot(Th1_list, gout_call_rate(13,:)', 'ks', Th1_list, gout_call_rate(14,:)', 'ko', Th1_list,  gout_call_rate(15,:)', 'kx','MarkerSize',11);
%ylim([0.7 1.0]) % y軸の範囲の設定
%xlabel('th gout','FontSize',18) % x軸：閾値th_gout
%ylabel('Generalout call blocking rate,r_{gout}','FontSize',18) % y軸：被災地外呼損率
%legend({'th gin=25','th gin=27','th gin=29'},[0.44 0.5 0.35 0.25],'FontSize',15)% 座標軸への凡例の追加

% figure(5);
%plot(Th1_list, gin_call_rate(13,:)', 'ks', Th1_list, gin_call_rate(14,:)', 'ko', Th1_list,  gin_call_rate(15,:)', 'kx','MarkerSize',11);
%ylim([0 1.0]) % y軸の範囲の設定
%xlabel('th gout','FontSize',18) % x軸：閾値th_gout
%ylabel('Generalin call blocking rate,r_{gin}','FontSize',18) % y軸：被災地呼損率
%legend({'th gin=25','th gin=27','th gin=29'},[0.44 0.5 0.35 0.25],'FontSize',15)% 座標軸への凡例の追加

% figure(6);
%plot(Th1_list, e_call_rate(13,:)', 'ks', Th1_list, e_call_rate(14,:)', 'ko', Th1_list,  e_call_rate(15,:)', 'kx','MarkerSize',11);
%ylim([0 0.15]) % y軸の範囲の設定
%xlabel('th gout','FontSize',18) % x軸：閾値th_gout
%ylabel('Emergency call blocking rate,r_e','FontSize',18) % y軸：緊急呼損率
%legend({'th gin=25','th gin=27','th gin=29'},[0.44 0.5 0.35 0.25],'FontSize',15)% 座標軸への凡例の追加

%figure(7);
%plot(Th1_list, sum_call_rate(17,:)', 'ks', Th1_list, sum_call_rate(18,:)', 'ko', Th1_list,  sum_call_rate(19,:)', 'kx','MarkerSize',11);
%ylim([0.8 1.1]) % y軸の範囲の設定
%xlabel('th gout','FontSize',18)  % x軸：閾値th_gout
%ylabel('sum call blocking rate','FontSize',18) % y軸：トータル呼損率
%legend({'th gin=16','th gin=17','th gin=18'},[0.44 0.5 0.35 0.25],'FontSize',15)% 座標軸への凡例の追加



 %閾値変化させたとき
% figure(5);
%plot(Th2_list, gout_call_rate(:,8)', 'p', Th1_list, gout_call_rate(:,9)', 'o', Th1_list,  gout_call_rate(:,10)', '*', 'MarkerSize',11);
%ylim([0 1.0])
%xlabel('th gout','FontSize',18)
%ylabel('generallout call blocking rate','FontSize',18)
%legend({'th gin=25','th gin=26','th gin=27','th gin=28','th gin=29','th gin=30'},[0.44 0.5 0.35 0.25],'FontSize',15)


% figure(6);
%plot(Th1_list, sum_call_rate(25,:)', 'rp', Th1_list, sum_call_rate(26,:)', 'go', Th1_list,  sum_call_rate(27,:)', 'b*',Th1_list, sum_call_rate(28,:)', 'rx',Th1_list, sum_call_rate(29,:)', 'g>',Th1_list, sum_call_rate(30,:)', 'b<', 'MarkerSize',11);
%ylim([0 1.0])
%xlabel('th gout','FontSize',18)
%ylabel('generallout call blocking rate','FontSize',18)
%legend({'th gin=25','th gin=26','th gin=27','th gin=28','th gin=29','th gin=30'},[0.44 0.5 0.35 0.25],'FontSize',15)


 
% figure(5);
%plot(Th_list, generall_r2(:,1,1,1) , 'kp', Th_list, generall_r2(:,1,2,1) , 'ko', Th_list,  generall_r2(:,1,3,1) , 'k*',Th_list, generall_r2(:,1,4,1) , 'kx',Th_list, generall_r2(:,1,5,1) , 'k>',Th_list, generall_r2(:,1,6,1) , 'k<','MarkerSize',8);
%xlabel('th')
%ylabel('generall call blocking rate')
% legend('yoyaku=0','yoyaku=1','yoyaku=2','yoyaku=3','yoyaku=5','yoyaku=10','location','best')