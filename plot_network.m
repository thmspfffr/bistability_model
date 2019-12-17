v = 3;

clear tmp1 tmp2 

for i = 0:6
  for bg = 0:4
    d1 = hdf5read(sprintf('~/pupmod/decision_network/proc/pupmod_decision_network_inp%d_bg%d_v%d.h5',i,bg,v),'rate_D1');
    d2 = hdf5read(sprintf('~/pupmod/decision_network/proc/pupmod_decision_network_inp%d_bg%d_v%d.h5',i,bg,v),'rate_D2');
    
    t = (1:3000)./10;
    for iseg = 1 : (15000/5)
      tmp1(iseg) = mean(d1((iseg-1)*5+1:(iseg)*5));
      tmp2(iseg) = mean(d2((iseg-1)*5+1:(iseg)*5));
    end
    d1 = tmp1'; d2 = tmp2';
    
    figure; set(gcf,'color','w');
    subplot(1,2,1); hold on; title(sprintf('i%d,bg%d: SR = %.3f',i,bg,mean(abs(diff(d1>d2)))))
    a=[d1 d2];
    a=round(a.*10)./10;
    [N,X]=hist3(a,'ctrs',{0:0.5:35,0:0.5:35});
    lab=round(X{1});
    
    imagesc(X{1},X{2},N,[0 10]); axis square; set(gca,'ydir','normal','tickdir','out')
    axis([0 35 0 35]); colormap(plasma)
    line([0 150],[0 150],'color','w','linestyle',':')
    xlabel('Firing rate (D1) [Hz]'); ylabel('Firing rate (D2) [Hz]')
    set(gca,'xtick',[0:5:35],'xticklabel',[0:5:35])
    tp_editplots
    
    subplot(1,2,2); hold on; title('Population Time series')
    plot(t(2:end),d1(2:end),'color',[0.5 0.5 0.5]);
    plot(t(2:end),d2(2:end),'k'); axis square
    axis([0 300 0 35]); set(gca,'tickdir','out');
    xlabel('Time [s]'); ylabel('Firing rate [Hz]')
    box on
    tp_editplots
    % clear d1 d2
    print(gcf,'-dpdf',sprintf('~/pupmod/decision_network/plots/decision_network_phaseplane_inp%d_bg%d_v%d.pdf',i,bg,v))
  end
end


%% INHIBITION MODULATION
v = 7

clear dur n k

for i_inh = 0:10
  for itr = 1
    
    if v == 7
      d1 = hdf5read(sprintf('~/pupmod/decision_network/proc/pupmod_decision_network_inh%d_inp0_bg0_tr%d_v%d.h5',i_inh,itr,v),'rate_D1');
      d2 = hdf5read(sprintf('~/pupmod/decision_network/proc/pupmod_decision_network_inh%d_inp0_bg0_tr%d_v%d.h5',i_inh,itr,v),'rate_D2');
      t = (1:6000)./10;

    else
      d1 = hdf5read(sprintf('~/pupmod/decision_network/proc/pupmod_decision_network_inh%d_inp0_bg0_v%d.h5',i_inh,v),'rate_D1');
      d2 = hdf5read(sprintf('~/pupmod/decision_network/proc/pupmod_decision_network_inh%d_inp0_bg0_v%d.h5',i_inh,v),'rate_D2');
      for iseg = 1 : (30000/5)
        tmp1(iseg) = mean(d1((iseg-1)*5+1:(iseg)*5));
        tmp2(iseg) = mean(d2((iseg-1)*5+1:(iseg)*5));
      end
      d1 = tmp1'; d2 = tmp2';
    end
    
    % compute durations
    % --------------
    k1=bwlabel((d1>d2)&(d1>15));
    k2=bwlabel((d2>d1)&(d2>15));
    for i=1:max(k1)
      tmp_cnt1(i)=sum(k1==i);
    end
    for i=1:max(k2)
      tmp_cnt2(i)=sum(k2==i);
    end
    dur{i_inh+1} = [tmp_cnt1 tmp_cnt2];
    % --------------
    
    % tmp = tp_dfa(d1,[1 60],50,0.5,50); dfa1(i_inh+1) = tmp.exp;
    % tmp = tp_dfa(d2,[1 60],50,0.5,50); dfa2(i_inh+1) = tmp.exp;
    
    
    figure; set(gcf,'color','w');
    subplot(1,2,1); hold on; title(sprintf('Inh%d: SR = %.3f',i_inh,mean(abs(diff(d1>d2)))))
    a=[d1 d2];
    a=round(a.*10)./10;
    [N,X]=hist3(a,'ctrs',{0:0.5:35,0:0.5:35});
    lab=round(X{1});
    
    imagesc(X{1},X{2},N,[0 10]); axis square; set(gca,'ydir','normal','tickdir','out')
    axis([0 35 0 35]); colormap(plasma)
    line([0 150],[0 150],'color','w','linestyle',':')
    xlabel('Firing rate (D1) [Hz]'); ylabel('Firing rate (D2) [Hz]')
    set(gca,'xtick',[0:5:35],'xticklabel',[0:5:35])
    set(gca,'ytick',[0:5:35],'yticklabel',[0:5:35])
    tp_editplots
    
    subplot(1,2,2); hold on; title('Population Time series')
    plot(t(2:end),d1(2:end),'b');
    plot(t(2:end),d2(2:end),'r'); axis square
    axis([0 600 0 35]); set(gca,'tickdir','out');
    xlabel('Time [20 ms segments]'); ylabel('Firing rate [Hz]')
    box on
    tp_editplots
    % clear d1 d2
    print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_decisionnetwork_inh%d_tr%d_v%d.pdf',i_inh,itr,v))
    
    [n(:,i_inh+1),k(:,i_inh+1)]=hist(dur{i_inh+1},100);
  end
end
%% DOMINANCE
v=7
for i_inh = 5
  for itr = 0
    tmp_cnt1 = 0;
    tmp_cnt2 = 0;
    d1 = hdf5read(sprintf('~/pupmod/decision_network/proc/pupmod_decision_network_inh%d_inp0_bg0_tr%d_v%d.h5',i_inh,itr,v),'rate_D1');
    d2 = hdf5read(sprintf('~/pupmod/decision_network/proc/pupmod_decision_network_inh%d_inp0_bg0_tr%d_v%d.h5',i_inh,itr,v),'rate_D2');
    
    % compute durations
    % --------------
    k1=bwlabel((d1>d2)&(d1>15));
    k2=bwlabel((d2>d1)&(d2>15));
    for i=1:max(k1)
      tmp_cnt1(i)=sum(k1==i);
    end
    for i=1:max(k2)
      tmp_cnt2(i)=sum(k2==i);
    end
  
    dur{i_inh+1}{itr+1} = [tmp_cnt1 tmp_cnt2];
    % --------------
    
    dom_dur(i_inh+1,itr+1) = mean(dur{i_inh+1}{itr+1})./10;
    count(i_inh+1,itr+1) = max(k1)+max(k2);
    
  end
end
%
s = std(count,[],2)/sqrt(size(count,2));
% eb(1,:) = mean(dom_dur,2)+s;
% eb(2,:) = mean(dom_dur,2)-s;
figure; set(gcf,'color','w'); 
subplot(1,2,1);hold on
shadedErrorBar([],mean(count,2),[s'; s']);
set(gca,'xtick',1:11,'xticklabel',[-0.5:0.1:0.5])
% set(gca,'ytick',-0.1:0.2:1.8,'yticklabel',num2cell(round((10.^(-0.1:0.2:1.8)*10))./10))
plot(6,mean(count(6,:),2),'o')
tp_editplots
xlabel('\Delta(Feedback inhibition) [in %]'); ylabel('Number of transitions')
axis([0 12 0 800]); axis square

print(gcf,'-dpdf',sprintf('~/pupmod/decision_network/plots/pupmod_decisionnetwork_domdur_v%d.pdf',v))

%%
dd=diff(d1>d2);
idx=find(dd==1);

for iidx = 1 : length(idx)
  
  idx=find(dd==1);
  cnt1 = sum(d1(idx+1)>15);
  
end
dd=diff(d2>d1);
idx=find(dd==1);

for iidx = 1 : length(idx)
  
  idx=find(dd==1);
  cnt2 = sum(d2(idx+1)>15);
  
end

cnt = cnt1+cnt2;