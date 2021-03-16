tmp.local= [];      for j= 1:3; tmp.local= [tmp.local       1:40]; end
tmp.block= [];      for j= 1:3; tmp.block= [tmp.block       ones(1,40).*j]; end
tmp.kskill= [];     for j= 1:3; tmp.kskill= [tmp.kskill     ones(1,40).*DATA(j).trials.otherK]; end
tmp.self= [];       for j= 1:3; tmp.self= [tmp.self         DATA(j).trials.self]; end
tmp.d= [];          for j= 1:3; tmp.d= [tmp.d               DATA(j).trials.d]; end
tmp.theta= [];      for j= 1:3; tmp.theta= [tmp.theta       DATA(j).trials.theta]; end
tmp.rs= [];         for j= 1:3; tmp.rs= [tmp.rs             DATA(j).trials.rs]; end
tmp.rg= [];         for j= 1:3; tmp.rg= [tmp.rg             DATA(j).trials.rg]; end
tmp.acc= [];        for j= 1:3; tmp.acc= [tmp.acc           DATA(j).typeI.correct]; end
tmp.gamble= [];     for j= 1:3; tmp.gamble= [tmp.gamble     DATA(j).typeII.response]; end
tmp.outcome= [];    for j= 1:3; tmp.outcome= [tmp.outcome   DATA(j).typeII.r]; end
tmp.rt1= [];        for j= 1:3; tmp.rt1= [tmp.rt1           DATA(j).typeI.response_time]; end
tmp.rt2= [];        for j= 1:3; tmp.rt2= [tmp.rt2           DATA(j).typeII.response_time]; end
tmp.resp1= [];      for j= 1:3; tmp.resp1= [tmp.resp1       DATA(j).typeI.response]; end
tmp.resp2= [];      for j= 1:3; tmp.resp2= [tmp.resp2       DATA(j).typeII.response]; end
tmp.trial= [];      for j= 1:3; tmp.trial= [tmp.trial       1:40]; end
tmp.incl= [];       
tmp.incl1=~isnan(tmp.resp1);        
tmp.incl2=~isnan(tmp.resp2);
tmp.incl3=~isinf(tmp.rt1);
tmp.incl4=~isinf(tmp.rt2);
tmp.incl=(tmp.incl1+tmp.incl2+tmp.incl3+tmp.incl4)==4;
idx_v= find(~tmp.incl);
for i_idx= 1:length(idx_v);
    if tmp.local(idx_v(i_idx))== 1;
        tmp.local(idx_v(i_idx)+1)= 1;
    end
end