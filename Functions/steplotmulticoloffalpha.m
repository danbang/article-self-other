function steplotmulticoloffalpha(x,col,width,offset,alpha)
% function steplotmulticoloffalpha(x,col,width,offset,alpha)

% initialise figure
% get columns
nC=size(x,2);    
% get mean
muX=squeeze(nanmean(x));
% get standard errors
steX=(nanstd(x,[],1))/sqrt(size(x,1)); %standard error is std/root n
% steX=(nanstd(x,[],1)).*1.96; % 95 CI

% plot bars
for n=1:nC;
    hold on;
    bar(n+offset-1,muX(n),'linewidth',width,'Facecolor',col(n,:),'Edgecolor',[0 0 0],'FaceAlpha',alpha);   % telling you what colour for the bars
    line([n+offset-1,n+offset-1],[muX(n)-(steX(n)/2),muX(n)+(steX(n)/2)],'color','k','linewidth',width);
end
% add reference line
plot([0 nC+1],[0 0],'k-','linewidth',width);
  
end