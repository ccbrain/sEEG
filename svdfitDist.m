function resDist=svdfitDist(a,n,P)


for i=1:size(P,1)
   resDist(i,1)=norm((a-P(i,:))-((a-P(i,:))*n)*n') ;
    
end