      %Fractional code for time series (See fig.2 in the manuscript) of Leech's heart interneuron model.

%%%% Here I refers to V^{shift}_{K2}  and u refers to the membrane voltage(v) of the manuscript.

%%% Title- Emergence of bursting in a network of memory dependent excitable and spiking Leech-Heart neurons
%%% Authors-Sanjeev Kumar Sharma; Argha Monda; Arnab Mondal; Ranjit Kumar Upadhyay and Chittaranjan Hens


clear;
 alpha=1;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set I
 I=-0.021;
%Set II
%  I=-0.015;
 %Set III
%  I=0.001;
 %Set IV
%  I=0.003;

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 dt = 0.001; tspan = 0:dt:30;
  %%%%%%%% Initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%
 u=rand/10;
 V=rand/10;
 w=rand/10;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uu=zeros(length(tspan),1); VV=zeros(length(tspan),1); ww=zeros(length(tspan),1);
uu(1,1)=u;
VV(1,1)=V;
ww(1,1)=w;

T1=tspan(end)/10;
T2=T1+60;
NN=length(tspan);
nn=1:NN-1;

WCoet=(NN+1-nn).^(1-alpha)-(NN-nn).^(1-alpha);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:length(tspan)-1    

    if j<2
  u=u+dt*(-2*(30*(w.^2)*(u+0.07)+8*(u+0.046)+200*(1/(1+exp(-150*(0.0305+u))).^3)*V*(u-0.045)));
   V=V+dt*(24.69*(1/(1+exp(500*(0.0333+u)))-V));
  w=w+dt*(4*(1/(1+exp(-83*(0.018+I+u)))-w));
  
    else
  %%%%% Fractional derivative starts  here%%%%%%%%%%%%%%%%%
  
       WCoe=WCoet(end-j+2:end);  % The weight   of the fractional drivative  at each  tiime t 
       kr = dt^alpha*gamma(2-alpha);     %  the kernel   from the fractional derivative and  weighted  the markovian term
        
  %%%%% Fractional derivative for u starts  here 
        
        d2dM=uu(1:j,:); % to call all past values of voltage u  for fractioanl integration
       
        TeDi=d2dM(2:j,:)-d2dM(1:j-1,:); % Delta uu (using all past values of u)  of  the  voltage memory trace of the fractional drivative  at each  tiime t 
        fraccalcu=WCoe*TeDi-d2dM(j,:);  %  The fraction derivative 
   
   %%%%% Fractional derivative of u ends   here 
        
    %%%%%   Fractional derivative for V starts  here 
    
        d2dMV=VV(1:j,1); % to call all past values  VV  for fractioanl integration
        TeDiV=d2dMV(2:j,1)-d2dMV(1:j-1,1); % Delta V (using all past values of VV) 
        fraccalcuV=WCoe*TeDiV-d2dMV(j,1);  %  The fraction derivative 
        
    %%%%% Fractional derivative ends V  here 
    %%%% ==w==  Fractional derivative for w starts  here 
      
        d2dMw=ww(1:j,1); % to call all past values  ww  for fractioanl integration
        TeDiw=d2dMw(2:j,1)-d2dMw(1:j-1,1); % Delta w (using all past values of ww) 
        fraccalcuw=WCoe*TeDiw-d2dMw(j,1);  %  The fraction derivative 
        %%%%% Fractional derivative ends w  here 
       
         u =kr*(-2*(30*(w.^2)*(u+0.07)+8*(u+0.046)+200*(1/(1+exp(-150*(0.0305+u))).^3)*V*(u-0.045)))- fraccalcu; 
         V =kr*(24.69*(1/(1+exp(500*(0.0333+u)))-V))-fraccalcuV;
         w=kr*(4*(1/(1+exp(-83*(0.018+I+u)))-w))-fraccalcuw;
        
        Memo2(j,:)=fraccalcu+uu(j,:); % to save memory over time
        Memo2V(j,:)=fraccalcuV+VV(j,:); % to save memory over time
        Memo2w(j,:)=fraccalcuw+ww(j,:);
    end
    
     uu(j+1,1)=u;
     VV(j+1,1)=V;
     ww(j+1,1)=w;
end
% figure

 plot(tspan,uu,'b');
 set(gca,'FontName','Times New Roman','FontSize',24)
 xlabel('{\bf{\it t}}','FontName', 'Times New Roman','FontSize',30,'Color','k', 'Interpreter', 'tex')
ylabel('{\bf{\it v}}','FontName', 'Times New Roman','FontSize',30,'Color','k', 'Interpreter', 'tex')
% figure;
%  plot(tspan,ww,'b');
% xlabel('t')
% ylabel('w')
%   figure
%   set(gca,'FontName','Times New Roman','FontSize',22)
%   plot(uu,VV);
%   xlabel('{\bf{\it u}}','Interpreter','tex','FontName','Times New Roman','FontSize',28)
% ylabel('{\bf{\it v}}','Interpreter','tex','FontName','Times New Roman','FontSize',28)
%   figure
%   set(gca,'FontName','Times New Roman','FontSize',22)
%   plot3(uu,VV,ww);
% xlabel('{\bf{\it u}}','Interpreter','tex','FontName','Times New Roman','FontSize',28)
% ylabel('{\bf{\it v}}','Interpreter','tex','FontName','Times New Roman','FontSize',28)
% zlabel('{\bf{\it w}}','Interpreter','tex','FontName','Times New Roman','FontSize',28)

% figure
% plot(tspan,uu,'b');
% xlabel('t')
% ylabel('u')
        
        
        
        