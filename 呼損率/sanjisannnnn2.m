% 臒l�𑍓����肷��v���O�����B�܂��A�e�Ẵg���q�b�N���x�����đ����A�ދ�������ݒ肵�Ă���B���݂͉��̔�����Ȃ����ׂĂ̏�Ԃɂ�����e�Ă̌đ�����\�����Ă���B�����A������̌đ����ȉ��̏ꍇ�̂ݎ����l���i�[�������ꍇ�́A100�`122�s�ڂ̃R�����g�}�[�N���������A�����Ŋe�������ݒ肷��B
% 2����t����(����)�g�[�^���đ�����r�O���t����эœK臒l�O���t
% 2017.07.24
% by Sumiko Miyata

%3����t����@�O���t�o��
%2018.09.19
%�쐣

% �����s��.�ɂ��ދ����������邽�߂�
% delta��mu3��ǉ�

clear;
close all;

param.B    = 20;    %�S�ш�
param.b1   = 1;      %�ً}�Đ�L�ш�
param.b2   = 1;      %��Вn�Đ�L�ш�
param.b3   = 1;      %��Вn�O�Đ�L�ш�

param.mu1  = 0.01; %�ً}�đދ���
param.mu2  = 0.01; %��Вn�đދ���
param.mu3  = 0.01; %��Вn�O�đދ���
ok=0;

c1=0.05; %�ً}�Ă̏���đ���
%c2=[1]; %��Вn�Ă̏���đ���
c2=0.5;

rho2d_list = [0.4];%[0.05:0.05:1];  %��Вn��0.05����P�܂�0.05�Ԋu�Ńg���q�b�N���x���v���b�g
rho1d_list = [0.5];     %[0.3 0.9]�@�@%�ً}�ăg���q�b�N���x
rho3d_list = [0.8];  %��Вn�O��0.05����P�܂�0.05�Ԋu�Ńg���q�b�N���x���v���b�g

%Th2_list    = 20;  %臒l2��0����S�ш�܂Œ��ׂ�
%Th1_list    = 20;  %臒l1��0����S�ш�܂Œ��ׂ�

 Th2_list    = [0:param.B];  %臒l2��0����S�ш�܂Œ��ׂ�
 Th1_list    = [0:param.B];  %臒l1��0����S�ш�܂Œ��ׂ�
 %Th2_list    = [0:20];  %臒l2��0����S�ш�܂Œ��ׂ�
 %Th1_list    = [0:20];  %臒l1��0����S�ш�܂Œ��ׂ�
H=1;

tic; %���ԑ���

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

%�����͕s�������C�ŗL�l�v�Z�̓r���ŃG���[�ɂȂ邱�Ƃ����邽�߁C�����C�G���[�Ŏ~�܂����ꍇ�ɂ́C���L���R�����g�A�E�g���āC
start_i =1;

%�����C�G���[�Ŏ~�܂����ꍇ�ɂ́C�����̃R�����g�A�E�g���͂���
%load temp.mat;
%start_i = i;

for i=start_i:numel(rho1d_list);  %�ً}�Ẵg���q�b�N���x
    
    for j=1:numel(rho2d_list);   %��Вn�Ẵg���q�b�N���x
        
        for l=1:numel(rho3d_list);  %��Вn�O�Ẵg���q�b�N���x
            
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
                for k=1:numel(Th2_list);   %臒l2��1����S�ш�܂�  n = numel(A)�͔z��A�̗v�f��n��Ԃ��B 
                     a=k;
                    for m=1:numel(Th1_list);   %臒l1��1����(臒l2)-1�܂�
                     b=m;   
                        param.rho1dash = rho1d_list(i);
                        param.rho2dash = rho2d_list(j);
                        param.rho3dash = rho3d_list(l);
                        param.B_th2     = Th2_list(k);
                        param.B_th1     = Th1_list(m);
                        [res_temp_e(k,m), res_temp_gin(k,m),res_temp_gout(k,m)] = sanjigenn(param);
                        %90�s�ځF�e�Ăɏd�݂Â��������]�����s�����Ƃ�������߂��B
                          %res_temp_sum(k,m)=0.5*res_temp_e(k,m)+0.4*res_temp_gin(k,m)+0.8*res_temp_gout(k,m);
                        
                        %�e�đ������e�ϐ��Ɋi�[
                        e_call_rate(k,m)=res_temp_e(k,m);
                        gin_call_rate(k,m)=res_temp_gin(k,m);
                        gout_call_rate(k,m)=res_temp_gout(k,m);
                        
                        
                        
                        %100�`122�s�ځF��Вn�Ă̌đ�������Вn�O�Ă�菬�����A��Вn�Ă̌đ�����c�Q�ȉ��ł���A�ً}�Ă̌đ�����c�P�ȉ��̎��͊e���l��ϐ��Ɋi�[�B�����łȂ��Ƃ��͂O�ɂ���B
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
                
                %128�`140�s�ځF�����𖞂���臒l�̐������̏󋵂ɂ����āA��Вn�O�Ă̌đ������ŏ��ɂȂ�臒l�y�ъe�Ă̌đ����𓱏o
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
                
                
                %��������
                %     if(res_temp_e(k,m)<c1&&res_temp_gin(k,m)<c2)           %�ً}�Ă�������̌đ����ȉ��Ȃ�
                %    res_gout(k,m) =res_temp_gout(k,m);  %������̌đ����ȉ��ً̋}�Ă̌đ����̓���臒l�ł̈�ʌĂ̌đ���
                %    else
                %      res_gout(k,m)=1;           %����ȊO�͔��f����Ȃ��悤��
                %  end
                %  if(res_e(k,m)<c)           %�ً}�Ă�������̌đ����ȉ��Ȃ�
                % gout_call_rate(k,m) =res_temp_gout(k,m);  %������̌đ����ȉ��ً̋}�Ă̌đ����̓���臒l�ł̈�ʌĂ̌đ���
                % else
                % gout_call_rate(k,m) =1;           %����ȊO�͔��f����Ȃ��悤��
                % end
                
                % emergency_r2(k,m,j,l,i) =res_temp_e(k,m);
                % generallin_r2(k,m,j,l,i) =res_temp_gin(k,m);
                % generallout_r2(k,m,j,l,i) =res_temp_gout(k,m);
                
                
                %   fprintf('rho1=%0.3f  rho2=%0.3f  rho3=%0.3f   \n', param.rho1dash, param.rho2dash, param.rho3dash);
                %    [res_min, min_ind] = min(gout_call_rate);    %��Вn�O�Ă̌đ����̍ŏ��l
                
                %         res_T(i,l,j)  = res_min;
                % r_opt(i,l,j)= res_min;  %�œK臒l�ł̌đ���
                %     Th2_opt(:,j,l,i) = Th2_list(min_ind);   %�œK臒l2
                %     Th1_opt(:,j,l,i) = Th1_list(min_ind);   %�œK臒l1
                %      th_opt(j)= Th_list(min_ind);
                %  emergency_r(:,j,l,i) =res_temp_e(min_ind);  %�œK臒l�łً̋}�Ă̌đ���
                %    r_1(j) = res_temp_n(min_ind);
                %    generallin_r(:,j,l,i) = res_temp_gin(min_ind);   %�œK臒l�ł̔�Вn�Ă̌đ���
                %    r_1(j) = res_temp_n(min_ind);
                %  generallout_r(:,j,l,i) = res_temp_gout(min_ind);   %�œK臒l�ł̔�Вn�O�Ă̌đ���
            end %�����܂�
            
            
        end
        
        
    end
    % i�̃��[�v�̍Ō�Ŗ��񃏁[�N�X�y�[�X���̑S�ϐ���ۑ�
    save temp.mat
    
end

%%

%for i=1:numel(rho1d_list);
%     rho1filename = round(param.rho1dash*100);
%     roptfile = sprintf('ropt%d.mat', rho1filename);
%     load(roptfile);
%end

toc;


 
 %�œK�ً}�Čđ�������
%ylim([0 0.5])
%xlabel('Traffic intencity \rho_2','FontSize',18)             %����
%ylabel('emergency call blocking rate','FontSize',18)          %�c��
 %legend({'m=0','m=1','m=2'},'FontSize',15)      %�E��̐�����
 

 %臒l��ω����������ً̋}�Čđ�������
figure
 x=Th1_list;
 y=Th2_list; 
 z=e_call_rate;
 stem3(x,y,z,'r','o','LineStyle','none') % 3�������U�f�[�^��̃v���b�g �}�[�J�[�͐Ԃ��~,���C���X�^�C���̓��C������
 zlim([0 1.0])
 xlabel('th gout','FontSize',18)
 ylabel('th gin','FontSize',18)
 zlabel('emergency call blocking late,r_e','FontSize',18)
 
 %臒l��ω����������̔�Вn�Čđ�������
 figure
 x=Th1_list;
 y=Th2_list; 
 z=gin_call_rate;
 stem3(x,y,z,'r','o','LineStyle','none') % 3�������U�f�[�^��̃v���b�g �}�[�J�[�͐Ԃ��~,���C���X�^�C���̓��C������
 zlim([0 1.0])
 xlabel('th gout','FontSize',18)
 ylabel('th gin','FontSize',18)
 zlabel('Generalin call blocking late,r_{gin}','FontSize',18)

 %臒l��ω����������ً̋}�Čđ�������
figure
 x=Th1_list;
 y=Th2_list; 
 z=gout_call_rate;
 stem3(x,y,z,'r','o','LineStyle','none') % 3�������U�f�[�^��̃v���b�g �}�[�J�[�͐Ԃ��~,���C���X�^�C���̓��C������
 zlim([0 1.0])
 xlabel('th gout','FontSize',18)
 ylabel('th gin','FontSize',18)
 zlabel('Generalout call blocking late,r_{gout}','FontSize',18)

 %臒l��ω����������̑S�Ă̌Ă̌đ�������
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
%legend({'��Вn�O��','��Вn��','�ً}��'},'FontSize',15)


%figure(1);
%plot(rho1d_list, res_T(:,1), 'ko', rho1d_list, res_T(:,2), 'k^', rho1d_list, res_T_c(:,1), 'k+', rho1d_list, res_T_c(:,2), 'k*', rho1d_list, res_T_o(:,1), 'kx', rho1d_list, res_T_o(:,2), 'ks','MarkerSize',8);
%ylim([0 0.3])
%xlabel('threshold th2','FontSize',30)
%ylabel('call blocking late','FontSize',30)
%legend('r^{\delta=0}_{opt}�C\rho_2%=0.3','r^{\delta=0}_{opt}�C\rho_2�f=0.9','r^{prop}_{opt}�C\rho_2�f=0.3','r^{prop}_{opt}�C\rho_2�f=0.9','r^{con}_{opt}�C\rho_2�f=0.3','r^{con}_{opt}�C\rho_2�f=0.9','Location','northwest')
%title('�g�[�^���đ�������')




 
 %�œK��Вn�Čđ�������
%  figure(2);
%plot(rho2d_list, generall_r(:,1,1), 'rp', rho2d_list,  generall_r(:,2,1), 'go', rho2d_list,  generall_r(:,3,1), 'b*','MarkerSize',15);
%ylim([0 1.0])
%xlabel('Traffic intencity \rho_2','FontSize',18)  
%ylabel('generall call blocking rate','FontSize',18)
%legend({'m=0','m=1','m=2'},'FontSize',15)
 
%�œK��Вn�O�Čđ�������
%  figure(2);
%plot(rho2d_list, generall_r(:,1,1), 'rp', rho2d_list,  generall_r(:,2,1), 'go', rho2d_list,  generall_r(:,3,1), 'b*','MarkerSize',15);
%ylim([0 1.0])
%xlabel('Traffic intencity \rho_2','FontSize',18)  
%ylabel('generall call blocking rate','FontSize',18)
%legend({'m=0','m=1','m=2'},'FontSize',15)


 %�œK臒l�đ�������
 %figure(3);
 %plot(rho2d_list,Th_opt(:,1,1), 'rp',rho2d_list,Th_opt(:,2,1), 'go',rho2d_list,Th_opt(:,3,1), 'b*','MarkerSize',15);
 %xlabel('Traffic intencity \rho_2','FontSize',18)
%ylabel('optimal threshold th_{opt}','FontSize',18)
% legend({'m=0','m=1','m=2'},'FontSize',15)
 %plot(Th1_list, gout_call_rate(:,:) , 'rp', Th1_list, gout_call_rate(2,:) , 'go', Th1_list,  gout_call_rate(3,:) , 'b*',Th1_list, gout_call_rate(4,:) , 'rx',Th1_list, gout_call_rate(5,:) , 'g>',Th1_list, gout_call_rate(6,:) , 'b<','MarkerSize',11);

 
 %臒l�ω��������Ƃ�
% figure(4);
%plot(Th1_list, gout_call_rate(13,:)', 'ks', Th1_list, gout_call_rate(14,:)', 'ko', Th1_list,  gout_call_rate(15,:)', 'kx','MarkerSize',11);
%ylim([0.7 1.0]) % y���͈̔͂̐ݒ�
%xlabel('th gout','FontSize',18) % x���F臒lth_gout
%ylabel('Generalout call blocking rate,r_{gout}','FontSize',18) % y���F��Вn�O�đ���
%legend({'th gin=25','th gin=27','th gin=29'},[0.44 0.5 0.35 0.25],'FontSize',15)% ���W���ւ̖}��̒ǉ�

% figure(5);
%plot(Th1_list, gin_call_rate(13,:)', 'ks', Th1_list, gin_call_rate(14,:)', 'ko', Th1_list,  gin_call_rate(15,:)', 'kx','MarkerSize',11);
%ylim([0 1.0]) % y���͈̔͂̐ݒ�
%xlabel('th gout','FontSize',18) % x���F臒lth_gout
%ylabel('Generalin call blocking rate,r_{gin}','FontSize',18) % y���F��Вn�đ���
%legend({'th gin=25','th gin=27','th gin=29'},[0.44 0.5 0.35 0.25],'FontSize',15)% ���W���ւ̖}��̒ǉ�

% figure(6);
%plot(Th1_list, e_call_rate(13,:)', 'ks', Th1_list, e_call_rate(14,:)', 'ko', Th1_list,  e_call_rate(15,:)', 'kx','MarkerSize',11);
%ylim([0 0.15]) % y���͈̔͂̐ݒ�
%xlabel('th gout','FontSize',18) % x���F臒lth_gout
%ylabel('Emergency call blocking rate,r_e','FontSize',18) % y���F�ً}�đ���
%legend({'th gin=25','th gin=27','th gin=29'},[0.44 0.5 0.35 0.25],'FontSize',15)% ���W���ւ̖}��̒ǉ�

%figure(7);
%plot(Th1_list, sum_call_rate(17,:)', 'ks', Th1_list, sum_call_rate(18,:)', 'ko', Th1_list,  sum_call_rate(19,:)', 'kx','MarkerSize',11);
%ylim([0.8 1.1]) % y���͈̔͂̐ݒ�
%xlabel('th gout','FontSize',18)  % x���F臒lth_gout
%ylabel('sum call blocking rate','FontSize',18) % y���F�g�[�^���đ���
%legend({'th gin=16','th gin=17','th gin=18'},[0.44 0.5 0.35 0.25],'FontSize',15)% ���W���ւ̖}��̒ǉ�



 %臒l�ω��������Ƃ�
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