%firing frequency profile
Timestamp=Timestamp(238:407);
Timestamp=Timestamp-Timestamp(1);

for j=1:16;
     c(j)=0
    for i=1:length(Timestamp);
   
    if Timestamp(i)<2+2*(j-1);
        c(j)=c(j)+1;
   
    end 
    end
end


for k=2:16;
    
  d(1)=c(1);
   d(k)=c(k)-c(k-1);
    
end
figure
bar(d,0.8)
plot(d,'*')

%average amplitude

PeakValley(1)=PeakValley(2);
a=mean(PeakValley)*0.195;

PeakValley1(1)=PeakValley1(2);
b=mean(PeakValley1)*0.195;

PeakValley2(1)=PeakValley2(2);
c=mean(PeakValley2)*0.195;

PeakValley3(1)=PeakValley3(2);
d=mean(PeakValley3)*0.195;

PeakValley4(1)=PeakValley4(2);
e=mean(PeakValley4)*0.195;

%triangulation
xa=1;ya=2;da=238;xb=2;yb=2;db=204;xc=1;yc=1;dc=202;xd=2;yd=1;dd=172;
syms x y   %f????
%--------------?????------------------------------------
f1='((xa-x)^2+(ya-y)^2)/((xb-x)^2+(yb-y)^2)=(db/da)^2';
% f2='((xc-x)^2+(yc-y)^2)/((xd-x)^2+(yd-y)^2)=(dd/dc)^2';
f3='((xa-x)^2+(ya-y)^2)/((xc-x)^2+(yc-y)^2)=(dc/da)^2';
f4='((xb-x)^2+(yb-y)^2)/((xc-x)^2+(yc-y)^2)=(dc/db)^2';
% f5='((xb-x)^2+(yb-y)^2)/((xd-x)^2+(yd-y)^2)=(dd/db)^2';
% f6='((xa-x)^2+(ya-y)^2)/((xd-x)^2+(yd-y)^2)=(dd/da)^2';
[xx,yy]=solve(f1,f3,f4,x,y); %???x,y???????????????????xx,yy

%????px(1),px(2)
px=eval(xx)
py=eval(yy);  %????py(1),py(2)
locx=px;
locy=py;

%plot neuron location
plot(px,py,'*');
hold on
plot(xa,ya,'o')
hold on
plot(xb,yb,'o')
hold on
plot(xc,yc,'o')
hold on
plot(xd,yd,'o')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rewrite images
I=imread('black-dot-md.png');
I=imresize(I,0.1)
imsave(I,'black-dot-md','h',1)
for i1=1:31;
    for j1=1:33;
        for k1=1:3;
            if I(i1,j1,k1)>=200;
                I(i1,j1,k1)=255
            end 
        end 
    end 
end

%Video making
T=round(Timestamp,2)
Ptime(1:1275)=1000;
for m=1:length(PeakValley)
Ptime(T*100)=PeakValley(m);
end
 p=Ptime(600:700)
T1=round(Timestamp1,2)
Ptime1(1:1275)=500;
for m=1:length(PeakValley1)
Ptime1(T1*100)=PeakValley1(m);
end
 p1=Ptime1(600:700)
  T2=round(Timestamp2,2)
Ptime2(1:1275)=500;
for m=1:length(PeakValley2)
Ptime2(T2*100)=PeakValley2(m);
end
 p2=Ptime2(600:700)
% p=[500,500,200,500,20
 T3=round(Timestamp3,2)
Ptime3(1:1275)=500;
for m=1:length(PeakValley3)
Ptime3(T3*100)=PeakValley3(m);
end
 p3=Ptime3(600:700)
 T4=round(Timestamp4,2)
 Ptime4(1:1275)=500;
for m=1:length(PeakValley4)
Ptime4(T4*100)=PeakValley4(m);
end
 p4=Ptime4(600:700) 
 T5=round(Timestamp5,2)
 Ptime5(1:1275)=500;
for m=1:length(PeakValley5)
Ptime5(T5*100)=PeakValley5(m);
end
 p5=Ptime5(600:700)
  T6=round(Timestamp6,2)
 Ptime6(1:1275)=600;
for m=1:length(PeakValley6)
Ptime6(T6*100)=PeakValley6(m);
end
 p6=Ptime6(600:700)
  T7=round(Timestamp7,2)
 Ptime7(1:1275)=500;
for m=1:length(PeakValley7)
Ptime7(T7*100)=PeakValley7(m);
end
 p7=Ptime7(600:700)

% p=[500,500,200,500,200,500];
I=imread('square.jpg');    
% I1=I;I2=I;I3=I;I4=I;I5=I;
I=imresize(I,0.1)
I1=imread('square.jpg');
I1=imresize(I1,0.1)
I2=imread('square.jpg');
I2=imresize(I2,0.1)
I3=imread('square.jpg');
I3=imresize(I3,0.1)
I4=imread('square.jpg');
I4=imresize(I4,0.1)
I5=imread('square.jpg');
I5=imresize(I5,0.1)
I6=imread('square.jpg');
I6=imresize(I6,0.1)
I7=imread('square.jpg');
I7=imresize(I7,0.1)
I8=imread('square.jpg');
I8=imresize(I8,0.1)
I9=imread('square.jpg');
I9=imresize(I9,0.1)
I10=imread('square.jpg');
I10=imresize(I10,0.1)
I11=imread('square.jpg');
I11=imresize(I11,0.1)
I12=imread('square.jpg');
I12=imresize(I12,0.1)
I13=imread('square.jpg');
I13=imresize(I13,0.1)
I14=imread('square.jpg');
I14=imresize(I14,0.1)
I15=imread('square.jpg');
I15=imresize(I15,0.1)
I1=imread('square.jpg');
I1=imresize(I1,0.1)
I1=imread('square.jpg');
I1=imresize(I1,0.1)

for k=1:length(p);
for i=1:21;
    for j=1:21;
       
            I(i,j,:)=p(k)*0.1
      
            I1(i,j,:)=p1(k)*0.2
          
            I2(i,j,:)=p2(k)*0.2
        
            I3(i,j,:)=p3(k)*0.2
         
            I4(i,j,:)=p4(k)*0.2
         
            I5(i,j,:)=p5(k)*0.2
          
            I6(i,j,:)=p6(k)*0.2
        
            I7(i,j,:)=p7(k)*0.2
            I8(i,j,:)=p8(k)*0.1
            I9(i,j,:)=p9(k)*0.1
            I10(i,j,:)=p10(k)*0.1
            I11(i,j,:)=p11(k)*0.1
            I12(i,j,:)=p12(k)*0.1
            I13(i,j,:)=p13(k)*0.1
            I14(i,j,:)=p14(k)*0.1
            I15(i,j,:)=p15(k)*0.1
            I16(i,j,:)=p16(k)*0.1
            I17(i,j,:)=p17(k)*0.1
           
            
         end
    end
end
          figure(k)
          subplot(6,6,1);
          imshow(I)
          subplot(6,6,8);
          imshow(I1);
          subplot(6,6,15);
          imshow(I2);
          subplot(6,6,18);
          imshow(I3);
          subplot(6,6,22);
          imshow(I4);
          subplot(6,6,23);
          imshow(I5);
           subplot(6,6,26);
          imshow(I6);
           subplot(6,6,28);
          imshow(I7);
          drawnow
          F(k)=getframe(gcf)
end
         
            



close all
movie(F,1,5)

movie2avi(F,'neuron firing 2','compression','none')

subplot(5,5,1); imshow(I);
subplot(5,5,8); imshow(I);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
MyTimes=h5read('awake_8-21-16_origin.kwik','/channel_groups/0/spikes/time_samples');
% time_fraction=h5read('awake_8-21-16_origin.kwik','/channel_groups/0/spikes/time_fractional');
% recording=h5read('awake_8-21-16_origin.kwik','/channel_groups/0/spikes/recording');
% fmask=h5read('awake_8-21-16_origin.kwx','/channel_groups/0/features_masks');
clust=h5read('awake_8-21-16_origin.kwik','/channel_groups/0/spikes/clusters/main');
fire=[MyTimes, clust];
time_short=fire(41439:42439,:);
time_short=single(time_short);
rfire_short=[round(time_short(:,1),-1)/10,time_short(:,2)];
rfire_short(:,1)=rfire_short(:,1)-354495;
ptime(1:6782)=0;
ptime=ptime';

for i=1:1001;
    ptime(rfire_short(i))=1;
    ptime(rfire_short(i),2)=rfire_short(i,2);
end



x=[-15,4.5,89,28,63,88,-8,38,40,-21,57,77];
y=[61,33,-1,85,24,74,45,103,-22,24,109,77];

% x=[101,13,80,90,104,76,24,57,30,18,12,83];
% y=[34,44,53,72,58,45,38,50,72,77,44,47];
for j=1:length(ptime);
     electrode=[0,3;1,3;2,3;3,3;0,2;1,2;2,2;3,2;0,1;1,1;2,1;3,1;1,0;2,0;3,0;0,0]
     electrode=electrode.*30;
xx=electrode(:,1);
yy=electrode(:,2);

figure
plot(xx,yy,'MarkerSize',60,'Marker','square','LineStyle','none')
xlim([-15,114]);
ylim([-30,114]);
hold on

axis off 
axis equal

 plot(x(1),y(1),'o','markersize',16,'MarkerFaceColor',.7*[0.494117647409439 0.184313729405403 0.556862771511078],'MarkerEdgeColor',.7*[0.494117647409439 0.184313729405403 0.556862771511078]);
            hold on
            plot(x(2),y(2),'o','markersize',16,'MarkerFaceColor',.7*[0.301960796117783 0.745098054409027 0.933333337306976],'MarkerEdgeColor',.7*[0.301960796117783 0.745098054409027 0.933333337306976]);
            hold on
           plot(x(3),y(3),'o','markersize',16,'MarkerFaceColor',.7*[0 1 0.74117648601532],'MarkerEdgeColor',.7*[0 1 0.74117648601532]);
            hold on
             plot(x(4),y(4),'o','markersize',16,'MarkerFaceColor',.7*[1 0 0],'MarkerEdgeColor',.7*[1 0 0]);
            hold on
            plot(x(5),y(5),'o','markersize',16,'MarkerFaceColor',.7*[1 0.674509823322296 0.18823529779911],'MarkerEdgeColor',.7*[1 0.674509823322296 0.18823529779911]);
            hold on
            plot(x(6),y(6),'o','markersize',16,'MarkerFaceColor',.7*[0.929411768913269 0.694117665290833 0.125490203499794],'MarkerEdgeColor',.7*[0.929411768913269 0.694117665290833 0.125490203499794]);
            hold on
             plot(x(7),y(7),'o','markersize',16,'MarkerFaceColor',.7*[0 0.447058826684952 0.74117648601532],'MarkerEdgeColor',.7*[0 0.447058826684952 0.74117648601532]);
            hold on
            plot(x(8),y(8),'o','markersize',16,'MarkerFaceColor',.7*[0 0 0.74117648601532],'MarkerEdgeColor',.7*[0 0 0.74117648601532]);
            hold on
             plot(x(9),y(9),'o','markersize',16,'MarkerFaceColor',.7*[0.635294139385223 0.0784313753247261 0.184313729405403],'MarkerEdgeColor',.7*[0.635294139385223 0.0784313753247261 0.184313729405403]);
            hold on
            plot(x(10),y(10),'o','markersize',16,'MarkerFaceColor',.7*[1 0.490196079015732 0],'MarkerEdgeColor',.7*[1 0.490196079015732 0]);
            hold on
            plot(x(11),y(11),'o','markersize',16,'MarkerFaceColor',.7*[1 0 1],'MarkerEdgeColor',.7*[1 0 1]);
            hold on
            plot(x(12),y(12),'o','markersize',16,'MarkerFaceColor',.7*[0.466666668653488 0.674509823322296 0.18823529779911],'MarkerEdgeColor',.7*[0.466666668653488 0.674509823322296 0.18823529779911]);
if ptime(j,2)==6
    plot(x(1),y(1),'o','markersize',16,'MarkerFaceColor',1*[0.494117647409439 0.184313729405403 0.556862771511078],'MarkerEdgeColor',.7*[0.494117647409439 0.184313729405403 0.556862771511078]);
   hold on
elseif ptime(j,2)==26
    plot(x(2),y(2),'o','markersize',16,'MarkerFaceColor',1*[0.301960796117783 0.745098054409027 0.933333337306976],'MarkerEdgeColor',.7*[0.301960796117783 0.745098054409027 0.933333337306976]);
    hold on
elseif ptime(j,2)==28
    plot(x(3),y(3),'o','markersize',16,'MarkerFaceColor',1*[0 1 0.74117648601532],'MarkerEdgeColor',.7*[0 1 0.74117648601532]);
    hold on
elseif ptime(j,2)==21
    plot(x(4),y(4),'o','markersize',16,'MarkerFaceColor',1*[1 0 0],'MarkerEdgeColor',.7*[1 0 0]);
    hold on
elseif ptime(j,2)==8
    plot(x(5),y(5),'o','markersize',16,'MarkerFaceColor',1*[1 0.674509823322296 0.18823529779911],'MarkerEdgeColor',.7*[1 0.674509823322296 0.18823529779911]);
    hold on
elseif ptime(j,2)==17
    plot(x(6),y(6),'o','markersize',16,'MarkerFaceColor',1*[0.929411768913269 0.694117665290833 0.125490203499794],'MarkerEdgeColor',.7*[0.929411768913269 0.694117665290833 0.125490203499794]);
    hold on
elseif ptime(j,2)==14
    plot(x(7),y(7),'o','markersize',16,'MarkerFaceColor',1*[0 0.447058826684952 0.74117648601532],'MarkerEdgeColor',.7*[0 0.447058826684952 0.74117648601532]);
    hold on
elseif ptime(j,2)==13
    plot(x(8),y(8),'o','markersize',16,'MarkerFaceColor',1*[0 0 0.74117648601532],'MarkerEdgeColor',.7*[0 0 0.74117648601532]);
    hold on
elseif ptime(j,2)==16
    plot(x(9),y(9),'o','markersize',16,'MarkerFaceColor',1*[0.635294139385223 0.0784313753247261 0.184313729405403],'MarkerEdgeColor',.7*[0.635294139385223 0.0784313753247261 0.184313729405403]);
    hold on
elseif ptime(j,2)==5
    plot(x(10),y(10),'o','markersize',16,'MarkerFaceColor',1*[1 0.490196079015732 0],'MarkerEdgeColor',.7*[1 0.490196079015732 0]);
    hold on
elseif ptime(j,2)==7
    plot(x(11),y(11),'o','markersize',16,'MarkerFaceColor',1*[1 0 1],'MarkerEdgeColor',.7*[1 0 1]);
    hold on
elseif ptime(j,2)==2
    plot(x(12),y(12),'o','markersize',16,'MarkerFaceColor',1*[0.466666668653488 0.674509823322296 0.18823529779911],'MarkerEdgeColor',.7*[0.466666668653488 0.674509823322296 0.18823529779911]);
    hold on
      
           
           
        
        end

    
       axis off
axis equal  
      drawnow
          F(j)=getframe(gcf)
       
          close all
end
   
         
close all
movie(F,1,5)
axis off
movie2avi(F,'neuron firing 333','compression','none','fps',60)
    movie2avi(F_trimmed,'neuron firing trimmed_new2','compression','none','fps',120)
    
t=MyTimes(40000:41000);
% v=zeros(1001,1);
v(1:1001)=1;
plot(t,v,'.')
plot(x(1),y(1),'o','markersize',20,'MarkerFaceColor',[1 0 0])
plot(x(1),y(1),'o','markersize',20,'makerfacecolor',[1 0 0])




