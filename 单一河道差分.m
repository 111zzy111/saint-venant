%模拟工况：单个河道，共有115个断面，长102493.5m，不考虑旁侧入流
%边界条件：上游水位为-5m，下游水位为-15m，初始水位-10m，初始流量为0。
%--------------------------------------------------

%--------------------------------------------------
%定义变量及初始化

clc,clear;

%导入断面原始测量数据，共115个断面，且每个断面给出了203个起点距和高程对
load('data.mat'); 

%设置大断面水位间距dh,T为计算时间段(s)，dt为时间步长（s）,a,c:连续方程离散后的系数，
%Zu、Zd分别为上下游水位，q0为河道初始流量，n0为糙率
dh=0.5;T=24*3600;dt=300;a=1;c=1;Zu=-5;Z0=-10;Zd=-15;q0=0;g=9.81;n0=0.025;

%为矩阵分配存储空间data数据表
%data(1,1)为断面总数115，n为河段总数，N为断面数*2即Zi,Qi变量总数115*2=230，t为时间分段数，x为Zi,Qi解向量 zeros表示1行N列，N个元素，初始化为0， ceil取整 
n=data(1,1)-1;N=2*data(1,1);t=ceil(T/dt);x=zeros(1,N);

%B：五对角矩阵方程的常数项，I,J,K,M,O:五对角矩阵A的五列  Ax=B
B=zeros(1,N);I=zeros(1,N);J=zeros(1,N-1);K=zeros(1,N-2);M=zeros(1,N-1);O=zeros(1,N-2);
 
%Z，Q分别为模拟时段内的水位和流量过程，X：水位、流量组合矩阵
Z=zeros(t+1,n+1);Q=zeros(t+1,n+1);X=zeros(N,t+1);

%x解向量初始化
for i=1:2:N
    x(i)=Z0;x(i+1)=q0;
end
%上下边界水位初始化
x(1)=Zu;x(N-1)=Zd;
%将初始时刻的解向量放入组合矩阵中存储
X(:,1)=x';

%初始化B,I,J,M向量
B(1)=Zu;B(N)=Zd;
I(1)=1;I(N)=0;
J(1)=0;M(N-1)=1;
%--------------------------------------------------

%--------------------------------------------------
%处理大断面实测数据data
%对应的水面宽度、水面面积、湿周、水力半径

%先转置原矩阵，再将矩阵各列合并成一个长的列向量 （'）转置，（：）改为列向量
AA=data';BB=AA(:);

%将列向量BB分成两列的矩阵Combinedata
%使若干个起点距和高程对能够对齐，将高程与起点距分开，BB（1，1）第一个断面起点距。BB(2，1）第一个断面高程，BB（3，1）第二个断面起点距 
for i=2:2:size(BB)
    Combinedata(i/2,1)=BB(i-1,1); 
    Combinedata(i/2,2)=BB(i,1);
end 
%处理断面间距，Combinedata(17,1)数据表示第一个断面距离河道起点的距离。17是在数据表中位于第17对 以此类推，第二个断面就是第247对
%nsection表示断面数量
nsection=Combinedata(1,1);
for i=1:1:nsection
    dx(i)=Combinedata(17+230*(i-1),1); 
end 
%差分求相邻断面间距,存储在dx向量中 diff是差分，计算相邻元素差值 
dx=diff(dx);

%生成每个断面的数据表
for i=1:1:nsection 
    
    %把每一个断面的起点距和高程数赋值给临时矩阵
    %第一个断面从第31对开始，第233对结束
    temp=Combinedata(31+230*(i-1):233+230*(i-1),:);
    
    %第二列的最小值及行号，即为大断面最低点
    [min_value,min_row]=min(temp(:,2));
    
    %设定工况：由于起点处有堤坝，而河道对岸没有堤坝，  所以大断面最后一个测量点的高程总是小于第一个  （不通用）
    %因此需要找到次最大值，即大断面最后一个测量点的高程，以便确定断面数据表的行数
    max_value=max(temp(min_row:end,2));
     
    %nrow表示断面数据表格的行数
    %Sct为断面数据总表
    %前五列分别为水位、水面宽、过水断面面积、湿周和水力半径，
    %往后就是第二个断面的水位、水面宽、过水断面面积、湿周和水力半径，以此类推
    nrow=ceil((max_value-min_value)/dh)+1;
    Sct(1:nrow,1+5*(i-1):5+5*(i-1))=zeros(nrow,5);
    
    for j=1:1:nrow 
        
        %水位
        z=min_value+dh*(j-1);
        
        %Plus_minus记录断面各点高程与所给水位的大小关系，用-1和0表示？？
        Plus_minus=diff(temp(:,2)>z);
        
        %记录下由正转负的零点位置
        bigtosmall=find(Plus_minus==-1);
        
        %记录下由负转正的零点位置
        smalltobig=find(Plus_minus==1);
        
        %判断正零点和负零点的个数是否相同，为了防止次最大值并不是对岸最后一点，如果最大值在最后一点右方某个凸起，应该在最右岸添加一个提防
        if numel(bigtosmall)>numel(smalltobig) 
            smalltobig(numel(bigtosmall),1)=size(temp,1);
            temp(size(temp,1)+1,1)=temp(end,1);temp(end,2)=z;
        end
        
        %生成断面数据表
        for k=1:1:size(bigtosmall)
            
            %求左边交点起点距
             dl=temp(bigtosmall(k)+1,1)-(z-temp(bigtosmall(k)+1,2))*(temp(bigtosmall(k)+1,1) ...
                -temp(bigtosmall(k),1))/(temp(bigtosmall(k),2)-temp(bigtosmall(k)+1,2));

            %求右边交点起点距
            dr=temp(smalltobig(k),1)+(z-temp(smalltobig(k),2))*(temp(smalltobig(k)+1,1) ...
                -temp(smalltobig(k),1))/(temp(smalltobig(k)+1,2)-temp(smalltobig(k),2));
            
            %Width表示相应水深的水面宽
            Width=dr-dl;
            
            %计算过水断面面积
            darea=(z-temp(bigtosmall(k)+1,2))*(temp(bigtosmall(k)+2,1)-dl)/2 ...
                +(z-temp(smalltobig(k),2))*(dr-temp(smalltobig(k)-1,1))/2;
            %累加梯形面积
            for m=bigtosmall(k)+2:1:smalltobig(k)-1 
                darea=darea+(z-temp(m,2))*(temp(m+1,1)-temp(m-1,1))/2;
            end
            
            %计算湿周
            dwetcycle=sqrt((temp(bigtosmall(k)+1,1)-dl)^2+(temp(bigtosmall(k)+1,2)-z)^2) ...
                +sqrt((temp(smalltobig(k),1)-dr)^2+(temp(smalltobig(k),2)-z)^2);
            %累加湿周
            for n=bigtosmall(k)+2:1:smalltobig(k) 
                dwetcycle=dwetcycle+sqrt((temp(n,1)-temp(n-1,1))^2+(temp(n,2)-temp(n-1,2))^2);
            end
            
            %分别为水位、水面宽、过水断面面积、湿周和水力半径
            Sct(j,1+5*(i-1))=z;
            Sct(j,2+5*(i-1))=Sct(j,2+5*(i-1))+Width;
            Sct(j,3+5*(i-1))=Sct(j,3+5*(i-1))+darea;
            Sct(j,4+5*(i-1))=Sct(j,4+5*(i-1))+dwetcycle;
            Sct(j,5+5*(i-1))=Sct(j,3+5*(i-1))/Sct(j,4+5*(i-1));
        end 
    end
end

%--------------------------------------------------
%分为t个时间层循环，求解五对角矩阵   为什么循环到N-2
for k=1:1:t  
    for i=1:2:N-2
        
        %利用差分确定前后断面与相应水位的大小关系
        p_m1=diff(Sct(:,1+5*(i-1)/2)>x(i));
        p_m2=diff(Sct(:,6+5*(i-1)/2)>x(i+2));
        %记录下前一个断面由负转正的零点位置
        stob1=find(p_m1==1);
        %记录下后一个断面由负转正的零点位置
        stob2=find(p_m2==1);
        
        %Zu1表示内插后上断面水位较小值，Zu2表示上断面水位较大值
        Zu1=Sct(stob1,1+5*(i-1)/2);
        Zu2=Sct(stob1+1,1+5*(i-1)/2);
        %x(i)表示上断面实际水位，alpha表示上断面的(Z-Z1)/(Z2-Z1)
        alpha=(x(i)-Zu1)/(Zu2-Zu1);
        
        %Zd1表示内插后下断面水位较小值，Zd2表示下断面水位较大值
        Zd1=Sct(stob2,6+5*(i-1)/2);
        Zd2=Sct(stob2+1,6+5*(i-1)/2);
        %x(i+2)表示下断面实际水位，beta表示下断面的(Z-Z1)/(Z2-Z1)
        beta=(x(i+2)-Zd1)/(Zd2-Zd1);
        
        %Bu1表示内插后上断面水面宽较小值，Bu2表示上断面水面宽较大值，Bu表示上断面实际水面宽
        Bu1=Sct(stob1,2+5*(i-1)/2);
        Bu2=Sct(stob1+1,2+5*(i-1)/2);
        Bu=Bu1+alpha*(Bu2-Bu1);
        %Bd1表示内插后下断面水面宽较小值，Bd2表示下断面水面宽较大值，Bd表示下断面实际水面宽
        Bd1=Sct(stob2,7+5*(i-1)/2);
        Bd2=Sct(stob2+1,7+5*(i-1)/2);
        Bd=Bd1+beta*(Bd2-Bd1);
        %Baverage表示上下断面水面宽平均值
        Baverage=(Bu+Bd)/2;
        
        %Au1表示内插后上断面过水面积较小值，Au2表示上断面过水面积较大值，Au表示上断面实际过水面积
        Au1=Sct(stob1,3+5*(i-1)/2);
        Au2=Sct(stob1+1,3+5*(i-1)/2);
        Au=Au1+alpha*(Au2-Au1);
        %Ad1表示内插后下断面过水面积较小值，Ad2表示下断面过水面积较大值，Ad表示下断面实际过水面积
        Ad1=Sct(stob2,8+5*(i-1)/2);
        Ad2=Sct(stob2+1,8+5*(i-1)/2);
        Ad=Ad1+beta*(Ad2-Ad1);
        %Aaverage表示上下断面过水面积平均值
        Aaverage=(Au+Ad)/2;
        
        %uu表示上断面流速，ud表示下断面流速，x(i+1)、x(i+3)分别表示Qi、Qi+1
        uu=x(i+1)/Au;ud=x(i+3)/Ad;
        %uaverage表示上下断面平均流速
        uaverage=(uu+ud)/2;
        
        %Ru1表示内插后上断面水力半径较小值，Ru2表示上断面水力半径较大值，Ru表示上断面实际水力半径
        Ru1=Sct(stob1,5+5*(i-1)/2);
        Ru2=Sct(stob1+1,5+5*(i-1)/2);
        Ru=Ru1+alpha*(Ru2-Ru1);
        %Rd1表示内插后下断面水力半径较小值，Rd2表示下断面水力半径较大值，Rd表示下断面实际水力半径
        Rd1=Sct(stob2,10+5*(i-1)/2);
        Rd2=Sct(stob2+1,10+5*(i-1)/2);
        Rd=Rd1+beta*(Rd2-Rd1);
        %Raverage表示上下断面水力半径平均值
        Raverage=(Ru+Rd)/2;
        
        %计算B向量各元素，B(i+1),B(i+2)分别为河段上下断面水位和流量之和即C1，C2
        B(i+1)=x(i)+x(i+2);
        B(i+2)=x(i+1)+x(i+3);
        
        %lambda表示网格比
        lambda=2*dt/dx((i+1)/2);
        %计算五对角阵对角线向量I,分别为B,G
        I(i+1)=-lambda/Baverage;I(i+2)=g*Aaverage*lambda;
        
        %计算对角线上一行向量J,分别为C,H
        J(i+1)=c;
        J(i+2)=1+uaverage*lambda+g*n0*n0*abs(uaverage)*dt/Raverage^(4/3);
        
        %计算对角线上两行向量K
        K(i)=0;K(i+1)=-I(i+1);
        
        %计算对角线下一行向量
        M(i)=a;M(i+1)=J(i+2)-2*uaverage*lambda;
        
        %O为对角线下两行向量.O(i)为方程组系数E,E=-G，O(i+1)为0已经初始化
        O(i)=-I(i+2);
    end
    
    %构造五对角矩阵A
    A=diag(I)+diag(J,1)+diag(K,2)+diag(M,-1)+diag(O,-2);
    %解五对角矩阵差分方程组并赋给X向量存储
    X(:,k+1)=A\B';
    %更新解向量
    x=X(:,k+1)';
end

%--------------------------------------------------
%从组合矩阵X中提取出水位Z和流量Q
%X的每一列表示一个时段的计算结果
%Z,Q的每一行表示一个时段的计算结果
for i=1:2:N
j=(i+1)/2;
Z(:,j)=X(i,:)';Q(:,j)=X(i+1,:)';
end  

%绘制出水位和流量三维图
figure(1);
mesh(Z);
xlabel('河段n');
ylabel('时间k/300s');
zlabel('水位/m');

figure(2);
mesh(Q);
xlabel('河段n');
ylabel('时间k/300s');
zlabel('流量(m^3/s)');

%-------------------------------------------------- 
