% rolling dice experiment Neural Networks course
clear all;
close all;
N = 3; % nr of dice
N_of_Exp = 10000; % Nr of experiments
x = zeros(N_of_Exp,1);

% Throw N dice, for N_of_Exp times
for n=1:N_of_Exp
    for i= 1 : N
        x(n) = x(n) + randi(6,1);
    end
end
figure;
hist(x, N*5+1);

% Computation of the theoretical distribution with N-1 repeated convolutions
dist_1d(1:6,1) = 1/6;
dist_all_dice = dist_1d;
for i=2 : N
    dist_all_dice = conv(dist_all_dice,dist_1d);
end

% Experimanetal distribution estimation 
x_histo = zeros(N*5+1,1);
for i = 1 : N*5+1
    for n=1:N_of_Exp 
        if  x(n) == N+i-1 
            x_histo(i) = x_histo(i) + 1;
        end
    end
end

for i = 1 : N*5+1
    fprintf('dice_point:%4d    est: %9.7f   true: %9.7f\n', i+N-1,   x_histo(i)/(N_of_Exp), dist_all_dice(i));
end

figure;
hold on;
grid on;
ylabel('p(x)');
xlabel('x');
plot(N : N*6,  dist_all_dice,'--','color','r', 'LineWidth',3);
plot(N : N*6,  x_histo./N_of_Exp,'color','k', 'LineWidth',2);

legend('Theoretical distribution','Estimate distribution');
    
