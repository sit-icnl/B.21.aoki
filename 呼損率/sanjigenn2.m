%function [emergencyr,generallrin,generallrout] = sanjigenn2(param)% paramを受け入れ引数とした収容本数における定常確率を出力として返すsanjigenを宣言
function sr = sanjigenn2(param)

rho1dash =param.rho1dash;  % 緊急呼の（正規化）トラヒック密度
rho2dash =param.rho2dash;  % 被災地呼の（正規化）トラヒック密度
rho3dash =param.rho3dash;  % 被災地外呼の（正規化）トラヒック密度
%rho1dash =0.4;%param.rho1dash;  % 緊急呼の（正規化）トラヒック密度
%rho2dash =0.5;%param.rho2dash;  % 被災地呼の（正規化）トラヒック密度
%rho3dash =0.8;%param.rho3dash;  % 被災地外呼の（正規化）トラヒック密度

B   =20;              % 全帯域(=サーバ数）
b1  =1;             % 緊急呼フローの占有帯域
b2  =1;             % 被災地呼フローの占有帯域
b3  =1;             % 被災地外呼フローの占有帯域

mu1 =0.1;            %  緊急呼の退去率(小さめの値にしたほうが安全）
mu2 =0.1;            %　被災地呼の退去率
mu3 =0.1;            %　被災地外呼の退去率

B_th1 =param.thout;        %被災地外呼の閾値1 
B_th2 =param.thin;        %被災地呼の閾値2 
%B_th1 =15;          %閾値1 
%B_th2 =17;          %閾値2

rho1 = rho1dash*B/b1;       % 緊急呼のトラヒック密度
rho2 = rho2dash*B/b2;       % 被災地呼のトラヒック密度
rho3 = rho3dash*B/b3;       % 被災地外呼のトラヒック密度

lambda1=mu1*rho1;           % 緊急呼の到着率
lambda2=mu2*rho2;           % 被災地呼の到着率
lambda3=mu3*rho3;           % 被災地外呼の到着率


n1max = floor(B/b1);     %緊急呼の本数の最大
n2max = floor(B/b2);     %被災地呼の本数の最大
n3max = floor(B/b3);     %被災地外呼の本数の最大

Nk=param.Nk;
Nh=param.Nh;
Ng=param.Ng;

%Nk=0;
%Nh=0;
%Ng=1;

%% 収容済み帯域→状態番号　状態番号→収容済み帯域　の写像をつくる
%for v=1:5
% 状態番号を格納する行列smatrixと，状態ごとの（各帯域の）収容済み呼数を格納するstateをつくる
Smatrix = ones(n1max+3, n2max+3,n3max+3)*NaN; % 状態数は最大呼数+1，さらに前後に番兵を設置するためもう+2する
statenum = 0;

for i=0:n1max
  for j=0:n2max
    for z=0:n3max
        Bnow = i*b1 + j*b2 + z*b3;       %いま考えている状態での仕様帯域
        if(Bnow <= B)                   %それが全帯域以下なら・・・
            statenum = statenum+1;  
            
            %まずは状態番号をインクリメント
            Smatrix(i+2,j+2,z+2)=statenum;  %Smatrixの適切な場所に状態番号を書き込み
                                        % 本数0の状態があるので+1，番兵の分さらに+1
            state(statenum).n1 = i;     %状態番号がstatenumのときの狭帯域の本数を保存
            state(statenum).n2 = j;     %同様
            state(statenum).n3 = z;
            state(statenum).Bnow = Bnow;   %一応仕様帯域も保存しておく．
            
        end
    end
  end
end
% Smatrixへのaccessor
statemat = @(i,j,z) Smatrix(i+2,j+2,z+2);     %これ以降，直接Smatrixは参照しないこと！

                                                                                                                                                                                                                                                                                                                
% 隣接行列Aをゼロで初期化
A=zeros(statenum,statenum);
A=sparse(A);

loss1=[];                               %緊急呼が呼損する状態番号のリストを空に初期化
loss2=[];                               %被災地呼が呼損する状態番号のリストを空に初期化
loss3=[];                               %被災地外呼が呼損する状態番号のリストを空に初期化
L=[];
% 隣接行列Aの中身を作成
for s=1:statenum                        %状態番号sをスタート地点とする
    n1=state(s).n1;                     %現在の状態番号における緊急呼の収容済み本数
    n2=state(s).n2;                     %被災地呼の収容済み本数
    n3=state(s).n3;                     %被災地外呼の収容済み本数
    Bnow = state(s).Bnow;

    if(n1==Nk && n2==Nh && n3==Ng) %収容数が引数としたものと同じ時
                L=[L s];%その状態番号を配列Lに控えておく
    end
    
    %---緊急呼到着
    t = statemat(n1+1,n2,n3);            %n1を一つ増やしてみたときの状態番号を行き先tとする．
    if(~isnan(t))%上のような状態が存在するかチェック．存在しないならt=NaNのはず
            A(t,s)=lambda1;                 %状態がある場合：スタート地点sから行き先tに遷移する確率はλ1
    else                                %網内マックス呼損
        loss1 = [loss1 s];              %状態が無い場合：いま到着した狭帯域は呼損するので，                                        
    end                                 %その番号を配列loss1に控えておく（後で呼損率の計算に使う）．
    
    %---被災地呼到着
    t = statemat(n1,n2+1,n3);
    if(~isnan(t))
        if(n1+n2+n3<B_th2)          %収容済回線数が閾値2を超えない時
            A(t,s)=lambda2;         %到着
            elseif(n1+n2+n3>=B_th2) %収容済回線数が閾値2を超える時
            loss2 =[loss2 s];       %呼損
        end
    else
            loss2 =[loss2 s]; 
    end

     %---被災地外呼到着
    t = statemat(n1,n2,n3+1);
    if(~isnan(t))
        if(n1+n2+n3<B_th1)          %収容済回線数が閾値1を超えない時
            A(t,s)=lambda3;         %到着
            elseif(n1+n2+n3>=B_th1) %収容済回線数が閾値1を超える時
            loss3 =[loss3 s];       %呼損
        end
    else
            loss3 =[loss3 s]; 
    end     

    %----緊急呼退去
    t = statemat(n1-1,n2,n3);
    if(~isnan(t))
        A(t,s)=n1*mu1;   
    end

    %----被災地呼退去
    t = statemat(n1,n2-1,n3);
    if(~isnan(t))
        A(t,s)=n2*mu2;
    end
    
    %----被災地外呼退去
    t = statemat(n1,n2,n3-1);
    if(~isnan(t))
        A(t,s)=n3*mu3;
    end
                
end
%% 各状態の定常確率の計算
P = A+diag(1-sum(A,1));                     %列和が1になるように対角要素をつくる．x=diag(A)は、Aの主対角要素の列ベクトルを返します。
[U, ~] = eigs(P,1);                
%最大固有値に対応する固有ベクトルを計算
Prob = U/sum(U);                            %確率なので，総和が1になるように定数倍

%if(Prob(L)>0)
    sr=Prob(L);
%else
 %   sr=0;

%end

%sr=Prob(L); %指定された収容数における定常確率



%% 呼損率の計算
                     %一般呼の呼損率を計算
%r1=sum(Prob(loss1));
%r2=sum(Prob(loss2));
%r3=sum(Prob(loss3));

%emergencyr=r1;
%generallrin=r2;
%generallrout=r3;

