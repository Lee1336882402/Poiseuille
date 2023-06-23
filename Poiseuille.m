clear
clc
Mx=21;
My=101;
dt=1e-4;
%X_length=0.1;
Y_length=1;
%dx=X_length/(Mx-1);
dy=Y_length/(My-1);
dx=dy;
Re=100;
U=[zeros(1,Mx);ones(My-2,Mx);zeros(1,Mx)];
V=zeros(My,Mx);
uin_all=1;
omega=([U(2:end,:);zeros(1,Mx)]-[zeros(1,Mx);U(1:end-1,:)])/dy-([V(:,2:end),V(:,2)]-[V(:,end-1),V(:,1:end-1)])/dx;
psi=repmat(linspace(0,uin_all,My)',1,Mx);
omega(1,:)=2*psi(2,:)/dy^2;
omega(My,:)=(2*psi(My-1,:)-2*psi(My,:))/dy^2;
step=0;
eps=3e-5;
while true
    step=step+1;
    U0=U;
    V0=V;
    omega0=omega;
    new_omega=omega;
    omega(1,:)=2*psi(2,:)/dy^2;
    omega(My,:)=(2*psi(My-1,:)-2*psi(My,:))/dy^2;
%     omega(2:My-1,1)=(psi(2:My-1,2)-2*psi(2:My-1,1)+psi(2:My-1,Mx-1))/dx^2+(psi(3:My,1)+psi(1:My-2,1)-2*psi(2:My-1,1))/dy^2;
%     omega(1,1)=(psi(1,2)-2*psi(1,1)+psi(1,Mx-1))/dx^2+2*(psi(2,1)-psi(1,1))/dy^2;
%     omega(My,1)=(psi(My,2)-2*psi(My,1)+psi(My,Mx-1))/dx^2+2*(psi(My-1,1)-psi(My,1))/dy^2;
%     omega(:,Mx)=omega(:,1);
    for i = 2:My-1
        for j = 1:Mx
            if j == 1
                new_omega(i,j)=(1/Re*(omega(i,j+1)+omega(i,Mx-1)-2*omega(i,j))/dx^2+1/Re*(omega(i+1,j)+omega(i-1,j)-2*omega(i,j))/dy^2 ...
                            -U(i,j)*(omega(i,j+1)-omega(i,Mx-1))/2/dx-V(i,j)*(omega(i+1,j)-omega(i-1,j))/2/dy)*dt+omega(i,j);
                continue
            end
            if j ==Mx
                new_omega(i,j)=(1/Re*(omega(i,2)+omega(i,j-1)-2*omega(i,j))/dx^2+1/Re*(omega(i+1,j)+omega(i-1,j)-2*omega(i,j))/dy^2 ...
                            -U(i,j)*(omega(i,2)-omega(i,j-1))/2/dx-V(i,j)*(omega(i+1,j)-omega(i-1,j))/2/dy)*dt+omega(i,j);
                continue
            end
            new_omega(i,j)=(1/Re*(omega(i,j+1)+omega(i,j-1)-2*omega(i,j))/dx^2+1/Re*(omega(i+1,j)+omega(i-1,j)-2*omega(i,j))/dy^2 ...
                            -U(i,j)*(omega(i,j+1)-omega(i,j-1))/2/dx-V(i,j)*(omega(i+1,j)-omega(i-1,j))/2/dy)*dt+omega(i,j);
        end
    end
    omega(2:My-1,1:Mx)=new_omega(2:My-1,1:Mx);
    instep=0;
    ineps=1e-3;
    while true
        instep=instep+1;
        psi0=psi;
        for i = 2:My-1
            for j = 1:Mx
                if j == 1
                    psi(i,j)=(psi(i,j+1)+psi(i,end-1)+psi(i+1,j)+psi(i-1,j))/4-dx^2/4*omega(i,j);
                    continue
                end
                if j ==Mx
                    psi(i,j)=(psi(i,2)+psi(i,j-1)+psi(i+1,j)+psi(i-1,j))/4-dx^2/4*omega(i,j);
                    continue
                end
                psi(i,j)=(psi(i,j+1)+psi(i,j-1)+psi(i+1,j)+psi(i-1,j))/4-dx^2/4*omega(i,j);
            end
        end
        if max(max(abs(psi-psi0)))<ineps
            break
        end
        if instep>100000
            error('内迭代错误！');
        end
    end
    V=-([psi(:,2:end),psi(:,2)]-[psi(:,end-1),psi(:,1:end-1)])/2/dx;
    U(2:end-1,:)=(psi(3:end,:)-psi(1:end-2,:))/2/dy;
    if(max(max(max(abs(U-U0))),max(max(abs(V-V0)))))<eps && max(max(abs(omega-omega0)))<eps
        break
    end
    if step>100000
        error('整体迭代步数过大！');
    end

    
end
%求压力场
P_source=-((([U(:,2:end),U(:,2)]-[U(:,end-1),U(:,1:end-1)])/2).^2+2*([U(2:end,:);2*U(end,:)-U(end-1,:)]-[2*U(1,:)-U(2,:);U(1:end-1,:)]).*([V(:,2:end),V(:,2)]-[V(:,end-1),V(:,1:end-1)])/4/dx/dy+(([V(2:end,:);2*V(end,:)-V(end-1,:)]-[2*V(1,:)-V(2,:);V(1:end-1,:)])/2/dy).^2);
P=ones(My,Mx);
instep=0;
ineps=1e-4;
while true
    P0=P;
    instep=instep+1;
    for i = 1:My
        for j = 1:Mx
            if j == 1
                if i ==1
                    P(i,j)=(P(i+1,j)+P(i,j)+P(i,Mx-1)+P(i,j+1))/4-dx^2/4*P_source(i,j);
                    continue
                end
                if i == My
                    P(i,j)=(P(2,j)+P(i-1,j)+P(i,Mx-1)+P(i,j+1))/4-dx^2/4*P_source(i,j);
                    continue
                end
                P(i,j)=(P(i+1,j)+P(i-1,j)+P(i,Mx-1)+P(i,j+1))/4-dx^2/4*P_source(i,j);
                continue
            end
            if j ==Mx
                if i ==1
                    P(i,j)=(P(i+1,j)+P(i,j)+P(i,j-1)+P(i,2))/4-dx^2/4*P_source(i,j);
                    continue
                end
                if i == My
                    P(i,j)=(P(2,j)+P(i-1,j)+P(i,j-1)+P(i,2))/4-dx^2/4*P_source(i,j);
                    continue
                end
                P(i,j)=(P(i+1,j)+P(i-1,j)+P(i,j-1)+P(i,2))/4-dx^2/4*P_source(i,j);
                continue
            end
            if i== 1
                P(i,j)=(P(i+1,j)+P(i,j)+P(i,j-1)+P(i,j+1))/4-dx^2/4*P_source(i,j);
                continue
            end
            if i == My
                P(i,j)=(P(i,j)+P(i-1,j)+P(i,j-1)+P(i,j+1))/4-dx^2/4*P_source(i,j);
                continue
            end
                
            P(i,j)=(P(i+1,j)+P(i-1,j)+P(i,j-1)+P(i,j+1))/4-dx^2/4*P_source(i,j);
        end
    end
    if max(max(abs(P-P0)))<ineps
        break
    end
    if instep>10000
        error('压力计算错误！');
    end
end

figure(1)
plot(U(:,1));
%figure(2)
%plot(U(:,25));
%[X ,Y]=meshgrid(linspace(0,1,Mx),linspace(0,1,My));
% figure(2)
% streamslice(U,V)

