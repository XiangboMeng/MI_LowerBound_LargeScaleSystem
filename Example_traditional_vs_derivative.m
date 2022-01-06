clear;
clc;

%% Parameters

% List of T values
Tl = [[1,2,3,4,5,6,8]*10,[1,2,3,4,5,6,8]*100,[1,2,3,4,5,6,8]*1000, [1]*10000];
% list of a values
al = [1,2,4,8,16];

%% Simulation
Results_recording_tri = [];           % For T, a, Ropt, and Topt; triditional method
Results_recording_der = [];           % For T, a, Ropt, and Topt; proposed method
count = 1;
for adx = 1:length(al)
for tdx = 1:length(Tl)
   
    a = al(adx);
    T = Tl(tdx);
    T_Results_recording_tri = 1:T;  % temp variable for brute force optimization
    % Traditional method with enumerative optimization
    R = (T-T_Results_recording_tri)/T.*(1 - (1- 1/(a*T)).^T_Results_recording_tri);
    [Ropt,Topt] = max(R);
    Results_recording_tri =[Results_recording_tri;[a,T,Ropt,Topt]];

    
    
    % proposed method 
    syms x
    eqn = -1 + exp(-x/a) + (exp(-x/a)*(1 - x))/a == 0;
    Tau_opt_sym = solve(eqn);
    Tau_opt = double(Tau_opt_sym);
    Ropt = (1- Tau_opt)*(1-exp(-Tau_opt/a));
    Results_recording_der = [Results_recording_der;  a,T,Ropt,Tau_opt];  

end
end

%% Data processing

T_Results_recording_tri =table(Results_recording_tri);
T_Results_recording_tri = splitvars(T_Results_recording_tri)
T_Results_recording_tri.Properties.VariableNames = {'a','T','Ropt','Topt'}


T_Results_recording_der =table(Results_recording_der);
T_Results_recording_der = splitvars(T_Results_recording_der)
T_Results_recording_der.Properties.VariableNames = {'a','T','Ropt','Tau_opt'}

%% Data plt
%color map
c_map = lines(length(al));

hf = figure;
L(1) = semilogx(nan, nan, 'k-',"Linewidth",1.5);
hold on;
L(2) = semilogx(nan, nan, 'k--',"Linewidth",1.5);
legend(L, {'first case', 'second case'})
for adx = 1:length(al)
    % plot the traditional
    rows = T_Results_recording_tri.a == al(adx);
    x = T_Results_recording_tri.T(rows);
    y = T_Results_recording_tri.Topt(rows)./x;

    htri = semilogx(x,y,"Linewidth",1.5,"Color",c_map(adx,:));
    hold on; grid on;
    
    
    % plot the proposed
    rows = T_Results_recording_der.a == al(adx);
    x = T_Results_recording_der.T(rows);
    y = T_Results_recording_der.Tau_opt(rows);
    hder(adx) = semilogx(x,y,'--',"Linewidth",1.5,"Color",c_map(adx,:));
    
    
end
set(0,'defaultTextInterpreter','latex');        %set default interpreter to be Latex
xlabel("$T$")
ylabel("$\tau_{\rm opt}$")
ylim([0.4,0.51])
set(gca, 'XLimSpec', 'Tight');
% extract the handles that require legend entries
% hleglines = [htri(1) hder(1)];
% create the legend
hleg = legend(L,{'Finite-$T$','Large-scale'},'Interpreter','latex');

% % arrow
% x = [0.7 0.8];
% y = [0.4 0.9];
% annotation('textarrow',x,y,'String','a =1,2,4,8,16','Interpreter','latex')

% labels
text(3*10^3,0.442-0.003,'$a =1$','Interpreter','latex','FontSize',14,'Color',c_map(1,:))

text(3*10^3, 0.470-0.003,'$a =2$','Interpreter','latex','FontSize',14,'Color',c_map(2,:))

text(3*10^3, 0.4847-0.003,'$a =4$','Interpreter','latex','FontSize',14,'Color',c_map(3,:))
text(3*10^3,0.4922-0.003,'$a =8$','Interpreter','latex','FontSize',14,'Color',c_map(4,:))
text(3*10^3,0.49611+0.003,'$a =16$','Interpreter','latex','FontSize',14,'Color',c_map(5,:))
% save("Example5_Triditional_vs_Derivative_for_a_T_Topt_Ropt.mat", 'T_Results_recording_der','T_Results_recording_tri')