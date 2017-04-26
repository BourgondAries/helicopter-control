figure_maxsize();
subplot(3,2,1);
plot(travel.time, travel.signals.values);
title('travel');
subplot(3,2,2);
plot(travel_rate.time, travel_rate.signals.values);
title('travel rate');
subplot(3,2,3);
plot(pitch.time, pitch.signals.values);
title('pitch');
subplot(3,2,4);
plot(pitch_rate.time, pitch_rate.signals.values);
title('pitch rate');
subplot(3,2,[5 6]);
plot(u(:,1), u(:,2));
title('u');
xlabel(strcat('time, q=', num2str(q)));
print(strcat('figures/10.2.3.q_', num2str(q)), '-depsc');

function figure_maxsize()
	figure('units','normalized','outerposition',[0 0 1 1]);
end
