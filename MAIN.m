
clear all
pop_M=100; % population size 100
pop_S = pop_M;
gen=500; % generation count 1000
selection_pressure = 'elitist'; % choose either 'elitist' or 'roulette wheel'
p_il = 0.5;  % probability of individual learning (BFGA quasi-Newton Algorithm) --> Indiviudal Learning is an IMPORTANT component of the MFEA.
rmp=0.3;   % random mating probability
reps = 20; % repetitions 20

Index=8;
for index = Index
    Tasks = mybenchmark(index);
%     Rs(index) = Similarity_check(Tasks,pop_M,gen,selection_pressure,rmp,p_il,reps);
    data_MFEA(index)=GMFEA(Tasks,pop_M,gen,selection_pressure,rmp,p_il); 

%     task_for_comparison_with_SOO = 1;
%     data_SOO_1(index)=SOEA(Tasks(task_for_comparison_with_SOO),pop_S,gen/2,selection_pressure,p_il,reps);   
% 
%     task_for_comparison_with_SOO = 2;
%     data_SOO_2(index)=SOEA(Tasks(task_for_comparison_with_SOO),pop_S,gen/2,selection_pressure,p_il,reps); 

end

 save('Rs.mat','Rs');
