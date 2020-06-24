%%% This is the network code for two incommensurate fractional order L-H model (See fig.4 in the manuscript)%%%%%%%%
%%%% Here I refers to V^{shift}_{K2} and u refers to the membrane voltage(v) of the manuscript.
%%To run this code put the axis_setting.m and this file in the same folder.
%%% Title- Emergence of bursting in a network of memory dependent excitable and spiking Leech-Heart neurons
%%%  Authors-Sanjeev Kumar Sharma; Argha Monda; Arnab Mondal; Ranjit Kumar Upadhyay and Chittaranjan Hens


clear All
 alpha=1;
 beta=0.65;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    I=0.001;
    D=0;   % %% Coupling strength
    N = 100;   %%% Total number of neuron
    M=60;   %%%% Number of oscillatory neurons.

%    E=[];
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 dt =0.001; tspan = 0:dt:30;

 uu=zeros(length(tspan),N); VV=zeros(length(tspan),N); ww=zeros(length(tspan),N);
%   uu1=zeros(length(tspan)-9999,N);
%   UU=zeros(length(tspan)-9997,N);
  uu1=zeros(length(tspan),N);
  UU=zeros(length(tspan),N);

 %%%%%%%% Initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%
 

for m=1:N
    u(m)=rand/10; 
    V(m)=rand/10;
    w(m)=rand/10;
 end 
 

  for p1=1:N
     
  uu(1,p1)=u(p1);
 VV(1,p1)=V(p1);
  ww(1,p1)=w(p1);
  
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% COUPLING MATRIX%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=zeros(N,N);
for i1=1:N
    for j1=1:N
        if (i1~=j1)
           a(i1,j1)=1;
        end
     end
end

 for i1=1:N
     a(i1,i1)=-sum(a(i1,:));
 end

    
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% T1=tspan(end)/10;
% T2=T1+60;
NN=length(tspan);
nn=1:NN-1;

WCoet=(NN+1-nn).^(1-alpha)-(NN-nn).^(1-alpha);
WCoet1=(NN+1-nn).^(1-beta)-(NN-nn).^(1-beta);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
for j=1:length(tspan)-1    

 for m1=1:N 
 for k=1:N
            coup(k)=u(1,:)*a(:,k);
 end

    if j<2
   u(m1)=u(m1)+dt*(-2*(30*((w(m1)).^2)*(u(m1)+0.07)+8*(u(m1)+0.046)+200*(1/(1+exp(-150*(0.0305+u(m1)))).^3)*V(m1)*(u(m1)-0.045))+(D/N)*coup(m1));
   V(m1)=V(m1)+dt*(24.69*(1/(1+exp(500*(0.0333+u(m1))))-V(m1)));
  w(m1)=w(m1)+dt*(4*(1/(1+exp(-83*(0.018+I+u(m1))))-w(m1)));
  
     
  
    else
  %%%%% Fractional derivative starts  here
       WCoe=WCoet(end-j+2:end);  % The weight   of the fractional drivative  at each  tiime t 
       kr = dt^alpha*gamma(2-alpha);     %  the kernel   from the fractional derivative and  weighted  the markovian term
       WCoe1=WCoet1(end-j+2:end);
       kr1 = dt^beta*gamma(2-beta);
        %%%%% Fractional derivative for u1 starts  here 
        if m1 <= M 
             
        TeDiu=uu(2:j,m1)-uu(1:j-1,m1); % Delta uu (using all past values of u)  of  the  voltage memory trace of the fractional drivative  at each  tiime t 
        fraccalcu(m1)=WCoe*TeDiu-uu(j,m1);  %  The fraction derivative 
         %%%%% Fractional derivative ends   here 
        
        %%%%% ==V1==  Fractional derivative for V1 starts  here 
       
        TeDiV=VV(2:j,m1)-VV(1:j-1,m1); % Delta V (using all past values of VV) 
        fraccalcuV(m1)=WCoe*TeDiV-VV(j,m1);  %  The fraction derivative 
        %%%%% Fractional derivative ends   here 
         %%%% ==w1==  Fractional derivative for w starts  here 
       
        TeDiw=ww(2:j,m1)-ww(1:j-1,m1); % Delta w (using all past values of ww) 
        fraccalcuw(m1)=WCoe*TeDiw-ww(j,m1);  %  The fraction derivative 
         
       
        %%%%% Fractional derivative ends   here 
   u(m1)=kr*((-2*(30*((w(m1)).^2)*(u(m1)+0.07)+8*(u(m1)+0.046)+200*(1/(1+exp(-150*(0.0305+u(m1)))).^3)*V(m1)*(u(m1)-0.045)))+(D/N)*coup(m1))-fraccalcu(m1);
   V(m1)=kr*((24.69*(1/(1+exp(500*(0.0333+u(m1))))-V(m1))))-fraccalcuV(m1);
  w(m1)=kr*((4*(1/(1+exp(-83*(0.018+I+u(m1))))-w(m1))))-fraccalcuw(m1);
   
        else
        TeDiu=uu(2:j,m1)-uu(1:j-1,m1); % Delta uu (using all past values of u)  of  the  voltage memory trace of the fractional drivative  at each  tiime t 
        fraccalcu(m1)=WCoe1*TeDiu-uu(j,m1);  %  The fraction derivative 
         %%%%% Fractional derivative ends   here 
        
        %%%%% ==V1==  Fractional derivative for V starts  here 
        
        TeDiV=VV(2:j,m1)-VV(1:j-1,m1); % Delta V (using all past values of VV) 
        fraccalcuV(m1)=WCoe1*TeDiV-VV(j,m1);  %  The fraction derivative 
        %%%%% Fractional derivative ends   here 
         %%%% ==w1==  Fractional derivative for w starts  here 
        
        TeDiw=ww(2:j,m1)-ww(1:j-1,m1); % Delta w (using all past values of ww) 
        fraccalcuw(m1)=WCoe1*TeDiw-ww(j,m1);  %  Th
        
        u(m1)=kr1*((-2*(30*((w(m1)).^2)*(u(m1)+0.07)+8*(u(m1)+0.046)+200*(1/(1+exp(-150*(0.0305+u(m1)))).^3)*V(m1)*(u(m1)-0.045)))+(D/N)*coup(m1))-fraccalcu(m1);
        V(m1)=kr1*((24.69*(1/(1+exp(500*(0.0333+u(m1))))-V(m1))))-fraccalcuV(m1);
        w(m1)=kr1*((4*(1/(1+exp(-83*(0.018+I+u(m1))))-w(m1))))-fraccalcuw(m1);
        end
        
%         Memo2um(j,:)=fraccalcum+uum(j,:); % to save memory over time
%         Memo2Vm(j,:)=fraccalcuVm+VVm(j,:); % to save memory over time
%         Memo2wm(j,:)=fraccalcuwm+wwm(j,:); 
        
        
%        
%   
   end
  
   
    uu(j+1,m1)=u(m1);
     VV(j+1,m1)=V(m1);
     ww(j+1,m1)=w(m1);
     
     
end

%  
% A=max(uu1(:,1));
% B=min(uu1(:,1));
% A1=max(uu7(:,1));
% B1=min(uu7(:,1));
% E1=[E1;D,A,B,A1,B1];
% E=[E;D,A,B,A1,B1];
    %S=sqrt(mean((u1-u2).^2)./sqrt(mean(u1.*u1)*mean(u2.*u2)));
   % A=[A;D,S];
%    end
end
x1=smooth(uu(1:end,1));
     x7=smooth(uu(1:end,M+1));


  plot(tspan,(x1),'b')
  hold on
  plot(tspan,(x7),'r')
 set(gca,'FontName','Times New Roman','FontSize',24)
xlabel(' \it \bf t ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')
ylabel(' \it \bf v_1, v_{61} ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')

   figure;
   for m1=1:N
uu1(:,m1)=uu(1:end,m1);
UU=[tspan',uu1(:,1:N)];
imagesc(UU(1:end,2:end))
   end

 
     set(gca,'YDir','normal')
     axis_setting