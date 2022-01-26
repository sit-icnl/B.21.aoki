%収容数等を引数としてパケットロス率を出力する関数
%pakettofunc.m というファイル内に、入力を受け入れ、パケットロス率を計算し、単一の結果を返す関数を定義。
function totalr = pakettofunc(param)

%z=1744;%[bit]              %G.711 1パケット当たり218[Byte]=1744[bit]
%B=100000000;%[bps]
%b0=100000;%[bps]     https://www.school.ctc-g.co.jp/columns/inst/inst50.html
%b1=100;             %通話あたりの音声帯域(正確には87.2kbpsだが、実環境では多少余裕を見た値を想定することが一般的)
%b2=100;

K=param.K;  %%%% キューの容量
mu=param.mu;%%%% 処理率

Nk=param.Nk; %緊急呼の収容数
Nh=param.Nh; %被災地呼の収容数
Ng=param.Ng; %被災地外呼の収容数

lamdas0k=param.lamdas0k; % 緊急呼のパケットの到着率 λs=β/T(α+β)
lamdas0h=param.lamdas0h; % 被災地呼のパケットの到着率
lamdas0g=param.lamdas0g; % 被災地外呼のパケットの到着率
%平方変動係数Ca^2:18.095、ひずみ度Sk:9.837
Dk=param.Dk;
Ek=param.Ek;
Fk=param.Fk;
Dh=param.Dh;
Eh=param.Eh;
Fh=param.Fh;
Dg=param.Dg;
Eg=param.Eg;
Fg=param.Fg;
%継続時間がそれぞれ平均r0^-1,r1^-1の指数分布に従う
r0k=Dk*(1-(1/sqrt(1+Nk*lamdas0k*Ek)));
r1k=Dk*(1+(1/sqrt(1+Nk*lamdas0k*Ek)));
r0h=Dh*(1-(1/sqrt(1+Nh*lamdas0h*Eh)));
r1h=Dh*(1+(1/sqrt(1+Nh*lamdas0h*Eh)));
r0g=Dg*(1-(1/sqrt(1+Ng*lamdas0g*Eg)));
r1g=Dg*(1+(1/sqrt(1+Ng*lamdas0g*Eg)));
%到着率がlamda0,lamda1のポアソン分布となる
lamdak0=Nk*lamdas0k+Fk+Fk*sqrt(1+Nk*lamdas0k*Ek);
lamdak1=Nk*lamdas0k+Fk-Fk*sqrt(1+Nk*lamdas0k*Ek);
lamdah0=Nh*lamdas0h+Fh+Fh*sqrt(1+Nh*lamdas0h*Eh);
lamdah1=Nh*lamdas0h+Fh-Fh*sqrt(1+Nh*lamdas0h*Eh);
lamdag0=Ng*lamdas0g+Fg+Fg*sqrt(1+Ng*lamdas0g*Eg);
lamdag1=Ng*lamdas0g+Fg-Fg*sqrt(1+Ng*lamdas0g*Eg);

smatrix = ones(K+3,K+3,K+3,K+3,K+3,K+3,K+3,K+3,11)*NaN;%%%%%%%%%%%%%%%%%,K+3,K+3,5
statenum=0;

for i=0:K
    for j=0:K
        for k=0:K
            for l=0:K
                for m=0:K
                    for n=0:K
                        for o=0:K
                            for p=0:K
                                for q=0:7
                                    Know=i+j+k+l+m+n+o+p;
                                    if(Know<=K)
                                        statenum=statenum+1;
                                        smatrix(i+2,j+2,k+2,l+2,m+2,n+2,o+2,p+2,q+2)=statenum;%%%%%%%%%%%%%
                                        state(statenum).q000=i;
                                        state(statenum).q001=j;
                                        state(statenum).q010=k;
                                        state(statenum).q011=l;
                                        state(statenum).q100=m;
                                        state(statenum).q101=n;
                                        state(statenum).q110=o;
                                        state(statenum).q111=p;
                                        state(statenum).q=q;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

statemat = @(i,j,k,l,m,n,o,p,q) smatrix(i+2,j+2,k+2,l+2,m+2,n+2,o+2,p+2,8);%%%%%%%%%%%%%%%% %
A=zeros(statenum,statenum);
A=sparse(A);

loss000=[];
loss001=[];
loss010=[];
loss011=[];
loss100=[];
loss101=[];
loss110=[];
loss111=[];

for s=1:statenum

    q000=state(s).q000;             %000
    q001=state(s).q001;             %001
    q010=state(s).q010;             %010
    q011=state(s).q011;             %011
    q100=state(s).q100;             %100
    q101=state(s).q101;             %101
    q110=state(s).q110;             %110
    q111=state(s).q111;             %111
    q=state(s).q;
    
    t = statemat(q000+1,q001,q010,q011,q100,q101,q110,q111,q);%%%%%%%%%%% %iを一つ増やしてみたときの状態番号を行き先tとする．
    if(~isnan(t))                       %上のような状態が存在するかチェック．存在しないならt=NaNのはず %TF=isnan(A)は、Aの要素がNaNの位置に1(true)を含み、要素がそれ以外である位置に 0 (false) を含むlogical配列を返します。
        A(t,s)=lamdak0+lamdah0+lamdag0;                 %状態がある場合：スタート地点sから行き先tに遷移する確率はλ0
    else                                %網内マックス呼損
        loss000 = [loss000 s];              %状態が無い場合：いま到着した狭帯域は呼損するので，                                        
    end
    
     t = statemat(q000,q001+1,q010,q011,q100,q101,q110,q111,q);%%%%%%%%%%%%
    if(~isnan(t))
        A(t,s)=lamdak0+lamdah0+lamdag1;
    else
        loss001 = [loss001 s];
    end
   
    t = statemat(q000,q001,q010+1,q011,q100,q101,q110,q111,q);%%%%%%%%%%%%,g0,g1,qk
    if(~isnan(t))
        A(t,s)=lamdak0+lamdah1+lamdag0;
    else
        loss010 = [loss010 s];
    end
    
    t = statemat(q000,q001,q010,q011+1,q100,q101,q110,q111,q);%%%%%%%%%%%%,g0,g1,qk
    if(~isnan(t))
        A(t,s)=lamdak0+lamdah1+lamdag1;
    else
        loss011 = [loss011 s];
    end
    
    t = statemat(q000,q001,q010,q011,q100+1,q101,q110,q111,q);%%%%%%%%%%%%,g0,g1,qk
    if(~isnan(t))
        A(t,s)=lamdak1+lamdah0+lamdag0;
    else
        loss100 = [loss100 s];
    end
    
    t = statemat(q000,q001,q010,q011,q100,q101+1,q110,q111,q);%%%%%%%%%%%%,g0,g1,qk
    if(~isnan(t))
        A(t,s)=lamdak1+lamdah0+lamdag1;
    else
        loss101 = [loss101 s];
    end
    
    t = statemat(q000,q001,q010,q011,q100,q101,q110+1,q111,q);%%%%%%%%%%%%,g0+1,g1,qk
    if(~isnan(t))
        A(t,s)=lamdak1+lamdah1+lamdag0;
    else
        loss110 = [loss110 s];
    end
    
    t = statemat(q000,q001,q010,q011,q100,q101,q110,q111+1,q);%%%%%%%%%%%%,g0,g1+1,qk
    if(~isnan(t))
        A(t,s)=lamdak1+lamdah1+lamdag1;
    else
        loss111 = [loss111 s];
    end
    
    %----退去
    t = statemat(q000-1,q001,q010,q011,q100,q101,q110,q111,q);%%%%%%%,g0,g1,qk
    if(~isnan(t))
        A(t,s)=q000*mu;   
    end
    
    t = statemat(q000,q001-1,q010,q011,q100,q101,q110,q111,q);%%%%%%%%%%%%
    if(~isnan(t))
        A(t,s)=q001*mu;   
    end
    
    t = statemat(q000,q001,q010-1,q011,q100,q101,q110,q111,q);%%%%%%%,g0,g1,qk
    if(~isnan(t))
        A(t,s)=q010*mu;   
    end
    
    t = statemat(q000,q001,q010,q011-1,q100,q101,q110,q111,q);%%%%%%%,g0,g1,qk
    if(~isnan(t))
        A(t,s)=q011*mu;   
    end
    
    t = statemat(q000,q001,q010,q011,q100-1,q101,q110,q111,q);%%%%%%%,g0,g1,qk
    if(~isnan(t))
        A(t,s)=q100*mu;   
    end
    
    t = statemat(q000,q001,q010,q011,q100,q101-1,q110,q111,q);%%%%%%%,g0,g1,qk
    if(~isnan(t))
        A(t,s)=q101*mu;
    end
    
    t = statemat(q000,q001,q010,q011,q100,q101,q110-1,q111,q);%%%%%%%
    if(~isnan(t))
        A(t,s)=q110*mu;  
    end
    
    t = statemat(q000,q001,q010,q011,q100,q101,q110,q111-1,q);%%%%%%%
    if(~isnan(t))
        A(t,s)=q111*mu; 
    end
    switch q
        case 0
            t = statemat(q000,q001,q010,q011,q100,q101,q110,q111,1);%%%%%%%+
            if(~isnan(t))
                A(t,s)=r0g; 
            end
            t = statemat(q000,q001,q010,q011,q100,q101,q110,q111,2);%%%%%%%+
            if(~isnan(t))
                A(t,s)=r0h;
            end
            t = statemat(q000,q001,q010,q011,q100,q101,q110,q111,4);%%%%%%%+
            if(~isnan(t))
                A(t,s)=r0k; 
            end
        case 1
            t = statemat(q000,q001,q010,q011,q100,q101,q110,q111,0);%%%%%%%+
            if(~isnan(t))
                A(t,s)=r1g; 
            end
            t = statemat(q000,q001,q010,q011,q100,q101,q110,q111,3);%%%%%%%+
            if(~isnan(t))
                A(t,s)=r0h; 
            end
            t = statemat(q000,q001,q010,q011,q100,q101,q110,q111,5);%%%%%%%+
            if(~isnan(t))
                A(t,s)=r0k; 
            end
        case 2
            t = statemat(q000,q001,q010,q011,q100,q101,q110,q111,0);%%%%%%%+
            if(~isnan(t))
                A(t,s)=r1h; 
            end
            t = statemat(q000,q001,q010,q011,q100,q101,q110,q111,3);%%%%%%%+
            if(~isnan(t))
                A(t,s)=r0g; 
            end
            t = statemat(q000,q001,q010,q011,q100,q101,q110,q111,6);%%%%%%%+
            if(~isnan(t))
                A(t,s)=r0k; 
            end
        case 3
            t = statemat(q000,q001,q010,q011,q100,q101,q110,q111,1);%%%%%%%+
            if(~isnan(t))
                A(t,s)=r1h; 
            end
            t = statemat(q000,q001,q010,q011,q100,q101,q110,q111,2);%%%%%%%+
            if(~isnan(t))
                A(t,s)=r1g; 
            end
            t = statemat(q000,q001,q010,q011,q100,q101,q110,q111,7);%%%%%%%+
            if(~isnan(t))
                A(t,s)=r0k; 
            end
        case 4
            t = statemat(q000,q001,q010,q011,q100,q101,q110,q111,0);%%%%%%%+
            if(~isnan(t))
                A(t,s)=r1k; 
            end
            t = statemat(q000,q001,q010,q011,q100,q101,q110,q111,5);%%%%%%%+
            if(~isnan(t))
                A(t,s)=r0g; 
            end
            t = statemat(q000,q001,q010,q011,q100,q101,q110,q111,6);%%%%%%%+
            if(~isnan(t))
                A(t,s)=r0h; 
            end
        case 5
            t = statemat(q000,q001,q010,q011,q100,q101,q110,q111,1);%%%%%%%+
            if(~isnan(t))
                A(t,s)=r1k; 
            end
            t = statemat(q000,q001,q010,q011,q100,q101,q110,q111,4);%%%%%%%+
            if(~isnan(t))
                A(t,s)=r0g; 
            end
            t = statemat(q000,q001,q010,q011,q100,q101,q110,q111,7);%%%%%%%+
            if(~isnan(t))
                A(t,s)=r1h; 
            end
        case 6
            t = statemat(q000,q001,q010,q011,q100,q101,q110,q111,2);%%%%%%%+
            if(~isnan(t))
                A(t,s)=r1k; 
            end
            t = statemat(q000,q001,q010,q011,q100,q101,q110,q111,4);%%%%%%%+
            if(~isnan(t))
                A(t,s)=r1h; 
            end
            t = statemat(q000,q001,q010,q011,q100,q101,q110,q111,7);%%%%%%%+
            if(~isnan(t))
                A(t,s)=r0g; 
            end
        case 7
            t = statemat(q000,q001,q010,q011,q100,q101,q110,q111,3);%%%%%%%+
            if(~isnan(t))
                A(t,s)=r1k; 
            end
            t = statemat(q000,q001,q010,q011,q100,q101,q110,q111,5);%%%%%%%+
            if(~isnan(t))
                A(t,s)=r1h; 
            end
            t = statemat(q000,q001,q010,q011,q100,q101,q110,q111,6);%%%%%%%+
            if(~isnan(t))
                A(t,s)=r1g; 
            end
        otherwise
            fprintf("error\n");
    end

end

P = A+diag(1-sum(A,1));     %列和が1になるように対角要素をつくる。D=diag(v,k)は、ベクトルvの要素をk番目の対角に配置します。k=0は主対角、k>0は主対角より上の対角、k<0は主対角より下の対角を表します。
[U V] = eigs(P,1);          %最大固有値に対応する固有ベクトルを計算 eigs(A,k)は、絶対値が最も大きいk個の固有値を返します。
Prob = U/sum(U);

r000 = sum(Prob(loss000));
r001 = sum(Prob(loss001));
r010 = sum(Prob(loss010));
r011 = sum(Prob(loss011));
r100 = sum(Prob(loss100));
r101 = sum(Prob(loss101));
r110 = sum(Prob(loss110));
r111 = sum(Prob(loss111));

totalr = ((r000+r001+r010+r011+r100+r101+r110+r111)/8);
end