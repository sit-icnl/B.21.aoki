%�O���V�~�����[�V���� �V�~�����[�V�����p�̃v���O�����B�e�Ẵg���q�b�N���x�̏�����ݒ肷��Ό��ʂ��o��B�������A�Ă̓��������Ɋւ��ẮA�e�Ă��Ɨ����ē�������悤�ɏC�������Ȃ���΂Ȃ�Ȃ��B
deltatime=0.01;  %�ŏ�����[�b/�^�C���X���b�g]

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
    rho1 =[0.5];%�ً}�ăg���q�b�N���x
    rho2 =[0.05:0.05:1];%��Вn�ăg���q�b�N���x
    %rho2 =[0.05];%��Вn�O�ăg���q�b�N���x
    rho3 =[0.8];
D=[950	850	800	750	700	700	700	700	700	700	650	650	650	650	650	650	650	650	650	650];%��Вn��臒l
 E=[600	550	600	600	600	600	650	650	650	650	600	600	600	600	600	600	600	600	600	600];%��Вn�O��臒l

%     D=1000;
%     E=1000;
%     O=[];
%     P=[];
%     Q=[];
    for p=1:numel(rho2)%臒l�ƃg���q�b�N���x�񂷗p
    D_a=D(p);
    E_a=E(p);
  rho2_a=rho2(p);
    for q=1:10%���s��
        
        

        %�z��Ȃǃ��Z�b�g
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
        
        for o=1:1000000 %���s��
            a=a+1;
            if(a==300000)%�ŏ��̂ق��̃f�[�^�����Z�b�g
                
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
            lambda1=1000*rho1*exprnd(0.01);%�������𓱏o
            e_l=lambda1;%���������L�^
            lambda2=1000*rho2_a*exprnd(0.01);
            gin_l=lambda2;
            lambda3=1000*rho3*exprnd(0.01);
            gout_l=lambda3;
            
            if(poissrnd(0.01*lambda1)~=0)
                e=e+1;%���������ً}�Ă����Z
                if(kaisen<1000)
                    e_ok=e_ok+1;%�i�[�����ً}�Ă����Z
                    kaisen4=kaisen+1;
                    e_t(e_ok)=(1000*rho1)/lambda1;%�������Ԃ��v�Z
                else
                    e_loss=e_loss+1;%�đ������ً}�Ă����Z
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
                e_t(t)=e_t(t)-0.01;%�������Ԃ���o�ߎ��Ԃ�����
                if(e_t(t)>0)
                    kaisen1=kaisen1+1;%�������I����ĂȂ��Ă𐔂���
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
        
        
        
        A(q)=e_loss/e;%�đ����v�Z
        B(q)=gin_loss/gin;
        C(q)=gout_loss/gout;
        
%         I(p)=e;
%         J(p)=gin;
%         K(p)=gout;
    end
    
%     L=sum(I)/1;
%     M=sum(J)/1;
%     N=sum(K)/1;
    
    F=sum(A)/10;%���s�񐔂Ŋ���
    G=sum(B)/10;
    H=sum(C)/10;
    
    O(p)=F;%�e�g���q�b�N�����̔��Ɋi�[
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


% %�ʏ�R���e���c�ƊȈՃR���e���c�̕��ϐ��̃O���t��ǉ�(�đ����̊m�F�y�эl�@�p)
% figure(2);
% plot(linkspeed,B,'ro-')%,linkspeed,B,'b^-');
% set(gca,'FontSize',14);
% xlabel('Linkspeed,L') % x-axis label
% ylabel('�ʏ�R���e���c�̌�(����)') % y-axis label
% legend('Conventional')%,'Proposed')
% 
% figure(3);
% plot(linkspeed,C,'ro-')%,linkspeed,B,'b^-');
% set(gca,'FontSize',14);
% xlabel('Linkspeed,L') % x-axis label
% ylabel('�ȈՃR���e���c�̌�(����)') % y-axis label
% legend('Conventional')%,'Proposed')

toc;


