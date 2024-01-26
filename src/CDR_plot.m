
figure
scatter(t(3:1799),predictions_residuals(3:1799),'b');
hold on
scatter(t2(3:1799),predictions2_residuals(3:1799),'r');
xlabel('Time into Flight','FontSize',14)
ylabel('Predicted Apogee - Actual Apogee (ft)','FontSize',14)
title('Residuals of Predictions','FontSize',20)
legend('TL 22 Residuals','CL 22 Residuals')
hold off