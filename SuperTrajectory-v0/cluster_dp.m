function [center_ins,cl] = cluster_dp(dist,center_num)


ND = size(dist,1);
dist = full(dist);
rho = sum(dist)+1;

[rho_sorted,ordrho]=sort(rho,'descend');


delta = zeros(1,ND);
nneigh = zeros(1,ND);
delta(ordrho(1))=rho_sorted(1);
nneigh(ordrho(1))=ordrho(1);

for ii=2:ND
   delta(ordrho(ii))=0;
   for jj=1:ii-1
     if(dist(ordrho(ii),ordrho(jj))>delta(ordrho(ii)))
        delta(ordrho(ii))=dist(ordrho(ii),ordrho(jj));
        nneigh(ordrho(ii))=ordrho(jj);
        
     end
   end
end
extra = size(find(nneigh==0),2);
delta(nneigh==0)=1e-300;

delta(ordrho(1))=min(delta(:));
NCLUST = max(center_num,1+extra);

gamma=rho./delta;
[~,ordgamma]=sort(gamma,'descend');


cl(1:ND)=-1;
cl(ordgamma(1:NCLUST))=1:NCLUST;
center_ins = ordgamma(1:NCLUST);
nneigh(ordgamma(1:NCLUST))=ordgamma(1:NCLUST);

%assignation
for i=1:ND
%     if nneigh(ordrho(i)) ==0
%         cl(ordrho(i))=0;
%     else
        cl(ordrho(i))=cl(nneigh(ordrho(i)));
%     end
end

