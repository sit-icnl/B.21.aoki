function [emergencyr,generallrin,generallrout] = sanjigenn(param)% param���󂯓���e�đ������o�͂Ƃ��ĕԂ�sanjigen��錾
% 3���������l�^��t����̃g�[�^���đ����̌v�Z �Ă̎�ނ𑝌������Ȃ�����A������K�v�Ȃ��Bsanjisannnnn2�ƃZ�b�g�ŉ^�p����B
% 2018.09.06
% by kawase

rho1dash =param.rho1dash;  % �ً}�Ắi���K���j�g���q�b�N���x
rho2dash =param.rho2dash;  % ��Вn�Ắi���K���j�g���q�b�N���x
rho3dash =param.rho3dash;  % ��Вn�O�Ắi���K���j�g���q�b�N���x

B   =param.B;              % �S�ш�(=�T�[�o���j
b1  =param.b1;             % �ً}�ăt���[�̐�L�ш�
b2  =param.b2;             % ��Вn�ăt���[�̐�L�ш�
b3  =param.b3;             % ��Вn�O�ăt���[�̐�L�ш�

mu1 =param.mu1;            %  �ً}�Ă̑ދ���(�����߂̒l�ɂ����ق������S�j
mu2 =param.mu2;            %�@��Вn�Ă̑ދ���
mu3 =param.mu3;            %�@��Вn�O�Ă̑ދ���
B_th1 =param.B_th1;        %臒l1 
B_th2 =param.B_th2;        %臒l2 

%clear;
%close all;
%rho1dash =0.5;%param.rho1dash;  % �ً}�Ắi���K���j�g���q�b�N���x
%rho2dash =0.5;%param.rho2dash;  % ��Вn�Ắi���K���j�g���q�b�N���x
%rho3dash =0.5;%param.rho3dash;  % ��Вn�O�Ắi���K���j�g���q�b�N���x

%B   =3;              % �S�ш�(=�T�[�o���j
%b1  =1;             % �ً}�ăt���[�̐�L�ш�
%b2  =1;             % ��Вn�ăt���[�̐�L�ш�
%b3  =1;             % ��Вn�O�ăt���[�̐�L�ш�

%mu1 =0.01;            %  �ً}�Ă̑ދ���(�����߂̒l�ɂ����ق������S�j
%mu2 =0.02;            %�@��Вn�Ă̑ދ���
%mu3 =0.03;            %�@��Вn�O�Ă̑ދ���
%B_th1 =1;          %臒l1 
%B_th2 =2;          %臒l2




rho1 = rho1dash*B/b1;       % �ً}�Ẵg���q�b�N���x
rho2 = rho2dash*B/b2;       % ��Вn�Ẵg���q�b�N���x
rho3 = rho3dash*B/b3;       % ��Вn�O�Ẵg���q�b�N���x

lambda1=mu1*rho1;           % �ً}�Ă̓�����
lambda2=mu2*rho2;           % ��Вn�Ă̓�����
lambda3=mu3*rho3;           % ��Вn�O�Ă̓�����


n1max = floor(B/b1);     %�ً}�Ă̖{���̍ő�
n2max = floor(B/b2);     %��Вn�Ă̖{���̍ő�
n3max = floor(B/b3);     %��Вn�O�Ă̖{���̍ő�


%% ���e�ςݑш恨��Ԕԍ��@��Ԕԍ������e�ςݑш�@�̎ʑ�������
%for v=1:5
% ��Ԕԍ����i�[����s��smatrix�ƁC��Ԃ��Ƃ́i�e�ш�́j���e�ς݌Đ����i�[����state������
Smatrix = ones(n1max+3, n2max+3,n3max+3)*NaN; % ��Ԑ��͍ő�Đ�+1�C����ɑO��ɔԕ���ݒu���邽�߂���+2����
statenum = 0;
% for j=0:B_th1
%     for i=0:B-j %���܍l���Ă����Ԃł̎d�l�ш�
%         Bnow=i*b1+j*b2+z*b3;
%        if(j+i+z<B_th1&&j+i+z<B_th2)  %���ꂪ�S�ш��臒l1�ȉ��Ȃ�E�E�E
%            statenum = statenum+1;      %�܂��͏�Ԕԍ����C���N�������g
%             Smatrix(i+2,j+2,z+2)=statenum;  %Smatrix�̓K�؂ȏꏊ�ɏ�Ԕԍ�����������
%                                        % �{��0�̏�Ԃ�����̂�+1�C�ԕ��̕������+1
%             state(statenum).n1 = i;     %��Ԕԍ���statenum�̂Ƃ��̋��ш�̖{����ۑ�
%             state(statenum).n2 = j;     %���l
%             state(statenum).n3 = z;
%             state(statenum).Bnow = Bnow;   %�ꉞ�d�l�ш���ۑ����Ă����D
%        elseif(j+i+z>=B_th1&&j+i+z<B_th2)
%                for z=0:n3max
%                     statenum=statenum+1;
%                          Smatrix(i+2,j+2,z+2)=statenum;  %Smatrix�̓K�؂ȏꏊ�ɏ�Ԕԍ�����������
%                                         % �{��0�̏�Ԃ�����̂�+1�C�ԕ��̕������+1
%             state(statenum).n1 = i;     %��Ԕԍ���statenum�̂Ƃ��̋��ш�̖{����ۑ�
%             state(statenum).n2 = j;     %���l
%             state(statenum).n3 = z;
%             state(statenum).Bnow = Bnow;   %�ꉞ�d�l�ш���ۑ����Ă����D
%                 end
%         end
%     end
% end


for i=0:n1max
  
   %if(i<=6)
   %    B=17;
   %elseif(i>6)
   %    B=16;
   %end
        for j=0:n2max
            for z=0:n3max
                Bnow = i*b1 + j*b2 + z*b3;       %���܍l���Ă����Ԃł̎d�l�ш�
                if(Bnow <= B)      %���ꂪ�S�ш�ȉ��Ȃ�E�E�E
                    statenum = statenum+1;      %�܂��͏�Ԕԍ����C���N�������g
                    Smatrix(i+2,j+2,z+2)=statenum;  %Smatrix�̓K�؂ȏꏊ�ɏ�Ԕԍ�����������
                                        % �{��0�̏�Ԃ�����̂�+1�C�ԕ��̕������+1
                    state(statenum).n1 = i;     %��Ԕԍ���statenum�̂Ƃ��̋��ш�̖{����ۑ�
                    state(statenum).n2 = j;     %���l
                    state(statenum).n3 = z;
                    state(statenum).Bnow = Bnow;   %�ꉞ�d�l�ш���ۑ����Ă����D
                end
            end
        end

end
% Smatrix�ւ�accessor
statemat = @(i,j,z) Smatrix(i+2,j+2,z+2);     %����ȍ~�C����Smatrix�͎Q�Ƃ��Ȃ����ƁI


                                                                                                                                                                                                                                                                                                                
% �אڍs��A���[���ŏ�����
A=zeros(statenum,statenum);
A=sparse(A);

loss1=[];                               %�ً}�Ă��đ������Ԕԍ��̃��X�g����ɏ�����
loss2=[];                               %��Вn�Ă��đ������Ԕԍ��̃��X�g����ɏ�����
loss3=[];                               %��Вn�O�Ă��đ������Ԕԍ��̃��X�g����ɏ�����

% �אڍs��A�̒��g���쐬
for s=1:statenum                        %��Ԕԍ�s���X�^�[�g�n�_�Ƃ���
    n1=state(s).n1;                     %���݂̏�Ԕԍ��ɂ�����ً}�Ă̎��e�ςݖ{��
    n2=state(s).n2;                     %��Вn�Ă̎��e�ςݖ{��
    n3=state(s).n3;                     %��Вn�O�Ă̎��e�ςݖ{��
    Bnow = state(s).Bnow;

    %---�ً}�ē���
    t = statemat(n1+1,n2,n3);            %n1������₵�Ă݂��Ƃ��̏�Ԕԍ����s����t�Ƃ���D
    if(~isnan(t))%��̂悤�ȏ�Ԃ����݂��邩�`�F�b�N�D���݂��Ȃ��Ȃ�t=NaN�̂͂�
   %     if(n1+n2+n3<18)
            A(t,s)=lambda1;                 %��Ԃ�����ꍇ�F�X�^�[�g�n�_s����s����t�ɑJ�ڂ���m���̓�1
   %         elseif(n1+n2+n3>=18)
   %         loss1 =[loss1 s]
   %     end
    else                                %�ԓ��}�b�N�X�đ�
        loss1 = [loss1 s];              %��Ԃ������ꍇ�F���ܓ����������ш�͌đ�����̂ŁC                                        
    end                                 %���̔ԍ���z��loss1�ɍT���Ă����i��Ōđ����̌v�Z�Ɏg���j�D
    
    %---��Вn�ē���
    t = statemat(n1,n2+1,n3);
    if(~isnan(t))
        if(n1+n2+n3<B_th2)          %���e�ω������臒l2�𒴂��Ȃ���
            A(t,s)=lambda2;         %����
            elseif(n1+n2+n3>=B_th2) %���e�ω������臒l2�𒴂��鎞
            loss2 =[loss2 s];       %�đ�
        end
    else
            loss2 =[loss2 s]; 
    end


     %---��Вn�O�ē���
    t = statemat(n1,n2,n3+1);
    if(~isnan(t))
        if(n1+n2+n3<B_th1)          %���e�ω������臒l1�𒴂��Ȃ���
            A(t,s)=lambda3;         %����
            elseif(n1+n2+n3>=B_th1) %���e�ω������臒l1�𒴂��鎞
            loss3 =[loss3 s];       %�đ�
        end
    else
            loss3 =[loss3 s]; 
    end
     
    


    %----�ً}�đދ�
    t = statemat(n1-1,n2,n3);
    if(~isnan(t))
        A(t,s)=n1*mu1;   
    end

    %----��Вn�đދ�
    t = statemat(n1,n2-1,n3);
    if(~isnan(t))
        A(t,s)=n2*mu2;
    end
    
    %----��Вn�O�đދ�
    t = statemat(n1,n2,n3-1);
    if(~isnan(t))
        A(t,s)=n3*mu3;
    end
            
end
%% �e��Ԃ̒��m���̌v�Z
P = A+diag(1-sum(A,1));                     %��a��1�ɂȂ�悤�ɑΊp�v�f������D
[U, V] = eigs(P,1);                
%�ő�ŗL�l�ɑΉ�����ŗL�x�N�g�����v�Z
Prob = U/sum(U);                            %�m���Ȃ̂ŁC���a��1�ɂȂ�悤�ɒ萔�{

%% �đ����̌v�Z
                     %��ʌĂ̌đ������v�Z
r1=sum(Prob(loss1));
r2=sum(Prob(loss2));
r3=sum(Prob(loss3));

emergencyr=r1;
generallrin=r2;
generallrout=r3;

