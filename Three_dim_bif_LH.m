%%%% Three dimensional bifurcation code of Leech's heart interneuron model (See fig.1 in the manuscript).
%%%% Here I refers to V^{shift}_{K2} and u refers to the membrane voltage(v) of the manuscript.
%%To run this code put the peakdet.m and this file in the same folder.
%%% Title- Emergence of bursting in a network of memory dependent excitable and spiking Leech-Heart neurons
%%%  Authors-Sanjeev Kumar Sharma; Argha Monda; Arnab Mondal; Ranjit Kumar Upadhyay and Chittaranjan Hens


for alpha=0.7:0.01:1
    alpha
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 dt = 0.001; tspan = 0:dt:5;
 %%%%%%%% Initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%
 u=rand/10;
 V=rand/10;
 w=rand/10;
%   u=(-0.0271517+rand); 
%   V=(0.0442+rand);
%   w=(0.0513+rand);
%   u=-0.0288274+rand;
%   V=0.0965+rand;
%   w=0.3067+rand;

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

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for I=-0.04:0.001:0.04
     I
for j=1:length(tspan)-1    
%     if (tspan(j)>T1) 
%      I=42;
%      %I=-99;
%     else 
%         I=0;
%     end;
    if j<2
  u=u+dt*(-2*(30*(w.^2)*(u+0.07)+8*(u+0.046)+200*(1/(1+exp(-150*(0.0305+u))).^3)*V*(u-0.045)));
   V=V+dt*(24.69*(1/(1+exp(500*(0.0333+u)))-V));
  w=w+dt*(4*(1/(1+exp(-83*(0.018+I+u)))-w));
  
    else
  %%%%% Fractional derivative starts  here
       WCoe=WCoet(end-j+2:end);  % The weight   of the fractional drivative  at each  tiime t 
       kr = dt^alpha*gamma(2-alpha);     %  the kernel   from the fractional derivative and  weighted  the markovian term
        %%%%% Fractional derivative for u starts  here 
        
        d2dM=uu(1:j,:); % to call all past values of voltage u  for fractioanl integration
       
        TeDi=d2dM(2:j,:)-d2dM(1:j-1,:); % Delta uu (using all past values of u)  of  the  voltage memory trace of the fractional drivative  at each  tiime t 
        fraccalcu=WCoe*TeDi-d2dM(j,:);  %  The fraction derivative 
         %%%%% Fractional derivative ends   here 
        
        %%%%% Fractional derivative for V starts  here 
        d2dMV=VV(1:j,1); % to call all past values  VV  for fractioanl integration
        TeDiV=d2dMV(2:j,1)-d2dMV(1:j-1,1); % Delta V (using all past values of VV) 
        fraccalcuV=WCoe*TeDiV-d2dMV(j,1);  %  The fraction derivative 
        %%%%% Fractional derivative ends   here 
             %%%% ==w==  Fractional derivative for w starts  here 
        d2dMw=ww(1:j,1); % to call all past values ww  for fractioanl integration
        TeDiw=d2dMw(2:j,1)-d2dMw(1:j-1,1); % Delta w (using all past values of ww) 
        fraccalcuw=WCoe*TeDiw-d2dMw(j,1);  %  The fraction derivative 
        %%%%% Fractional derivative ends   here 
       
         u =kr*(-2*(30*(w.^2)*(u+0.07)+8*(u+0.046)+200*(1/(1+exp(-150*(0.0305+u))).^3)*V*(u-0.045)))- fraccalcu; 
         V =kr*(24.69*(1/(1+exp(500*(0.0333+u)))-V))-fraccalcuV;
         w=kr*(4*(1/(1+exp(-83*(0.018+I+u)))-w))-fraccalcuw;
        
%         Memo2(j,:)=fraccalcu+uu(j,:); % to save memory over time
%         Memo2V(j,:)=fraccalcuV+VV(j,:); % to save memory over time
%         Memo2w(j,:)=fraccalcuw+ww(j,:);
    end
  
    uu(j+1,1)=u;
     VV(j+1,1)=V;
     ww(j+1,1)=w;
end
y(:,1)=uu(:,1);
y(:,2)=VV(:,1);
y(:,3)=ww(:,1);
%%%%%%%%%%%%%%%%%%%%%
a=size(y(:,1))
b=round(3*a/4)
x=y(b:end,1)
[maxtab,mintab]=peakdet(x,0.001)
%%%%%%%%%%%%%%%%%%%%%%%%%% Plot command for different slices %%%%%%%%%%%%% 
if(alpha==0.82)
    
if((size(maxtab)~=[0 0]) & (size(mintab)~=[0 0]))
    set(gca,'FontName','Times New Roman','FontSize',24)
plot3(alpha,I,maxtab(:,2),'g .','markersize',2)
   xlabel(' \it \bf \alpha ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')
ylabel(' \it \bf V^{shift}_{K2} ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')
   zlabel(' \it \bf v ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')
hold on
plot3(alpha,I,mintab(:,2) ,'g .','markersize',2)
 xlabel(' \it \bf \alpha ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')
ylabel(' \it \bf V^{shift}_{K2} ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')
   zlabel(' \it \bf v ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')
else
     y1=y(size(y,1),:);
     y2=y1(:,1);
plot3(alpha,I,y2, 'g .','markersize',2)
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
if(alpha==0.87)
     
if((size(maxtab)~=[0 0]) & (size(mintab)~=[0 0]))
    set(gca,'FontName','Times New Roman','FontSize',24)
plot3(alpha,I,maxtab(:,2),'r .','markersize',2)
   xlabel(' \it \bf \alpha ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')
ylabel(' \it \bf V^{shift}_{K2} ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')
   zlabel(' \it \bf v ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')
hold on
plot3(alpha,I,mintab(:,2) ,'r .','markersize',2)
 xlabel(' \it \bf \alpha ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')
ylabel(' \it \bf V^{shift}_{K2} ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')
   zlabel(' \it \bf v ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')
else
     y1=y(size(y,1),:);
     y2=y1(:,1);
plot3(alpha,I,y2, 'r .','markersize',2)
end
end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if(alpha==0.94)
    
if((size(maxtab)~=[0 0]) & (size(mintab)~=[0 0]))
    set(gca,'FontName','Times New Roman','FontSize',24)
plot3(alpha,I,maxtab(:,2),'b .','markersize',2)
   xlabel(' \it \bf \alpha ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')
ylabel(' \it \bf V^{shift}_{K2} ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')
   zlabel(' \it \bf v ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')
hold on
plot3(alpha,I,mintab(:,2) ,'b .','markersize',2)
 xlabel(' \it \bf \alpha ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')
ylabel(' \it \bf V^{shift}_{K2} ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')
   zlabel(' \it \bf v ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')
else
     y1=y(size(y,1),:);
     y2=y1(:,1);
plot3(alpha,I,y2, 'b .','markersize',2)
end

   end
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
   if(alpha==1)
    
if((size(maxtab)~=[0 0]) & (size(mintab)~=[0 0]))
    set(gca,'FontName','Times New Roman','FontSize',24)
plot3(alpha,I,maxtab(:,2),'k .','markersize',2)
   xlabel(' \it \bf \alpha ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')
ylabel(' \it \bf V^{shift}_{K2} ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')
   zlabel(' \it \bf v ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')
hold on
plot3(alpha,I,mintab(:,2) ,'k .','markersize',2)
 xlabel(' \it \bf \alpha ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')
ylabel(' \it \bf V^{shift}_{K2} ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')
   zlabel(' \it \bf v ', 'FontName', 'Times New Roman', ...
       'FontSize',30,'Color','k', 'Interpreter', 'tex')
else
     y1=y(size(y,1),:);
     y2=y1(:,1);
plot3(alpha,I,y2, 'k .','markersize',2)
end

  end
 end
end

