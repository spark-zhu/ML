clear all
files = dir('*.csv')
cor=[]
for file=1:10
data=importVideoSeq(files(file).name)
cor=[cor;data]
end

scale=1.255;
cor_scale=cor*scale;
plot((1:1797)/3,cor_scale(:,1))
hold on
plot((1:1797)/3,cor_scale(:,2))
xlabel('s')
ylabel('um')
legend('X','Y')