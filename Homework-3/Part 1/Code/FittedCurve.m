OriginalData = importdata('atkinson.dat');

FittedDataP3 = importdata('solnVecP3.dat');

figure(1);
hold on
plot(OriginalData.data(:,1),FittedDataP3,'r');
plot(OriginalData.data(:,1),OriginalData.data(:,2),'o','MarkerFaceColor','b');
legend('3rd Degree Polynoimal', 'Atkinson Data');
title('3rd Degree polynomial fitted data CD');
saveas(figure(1),'3DegPolyCD.png')
hold off


FittedDataP5 = importdata('solnVecP5.dat');

figure(2);
hold on
plot(OriginalData.data(:,1),FittedDataP5,'r');
plot(OriginalData.data(:,1),OriginalData.data(:,2),'o','MarkerFaceColor','b');
legend('5th Degree Polynoimal', 'Atkinson Data');
title('5th Degree polynomial fitted data CD');
saveas(figure(2),'5DegPolyCD.png')
hold off


FittedDataP3QR = importdata('solnVecP3QR.dat');

figure(3)
hold on
plot(OriginalData.data(:,1),FittedDataP3QR,'r');
plot(OriginalData.data(:,1),OriginalData.data(:,2),'o','MarkerFaceColor','b');
legend('3rd Degree Polynoimal', 'Atkinson Data');
title('3rd Degree polynomial fitted data QR');
saveas(figure(3),'3DegPolyQR.png')
hold off


FittedDataP5QR = importdata('solnVecP5QR.dat');

figure(4)
hold on
plot(OriginalData.data(:,1),FittedDataP5QR,'r');
plot(OriginalData.data(:,1),OriginalData.data(:,2),'o','MarkerFaceColor','b');
legend('5th Degree Polynoimal', 'Atkinson Data');
title('5th Degree polynomial fitted data QR');
saveas(figure(4),'5DegPolyQR.png')
hold off