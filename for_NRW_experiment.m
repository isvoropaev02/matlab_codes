d = 0.002; % размер образца
a = 0.008636; % размер большей стенки волновода
c = 3 * 10.^8; % скорость света
f_0 = c / (2 * a); % критическая частота низшей моды

L1 = 0.012;
L2 = 0.014;

%% data import
table11 = readtable("s11_assym.txt");
table21 = readtable("s21_assym.txt");
table12 = readtable("s12_assym.txt");
table22 = readtable("s22_assym.txt");


frequency = table11{:, 1} .* 10.^9; % В ЭКСПЕРИМЕНТЕ ИСПРАВИТЬ ЧАСТОТУ

ReS11 = table11{:, 2};
ImS11 = table11{:, 3};
s_11 = ReS11 + 1i * ImS11;

ReS21 = table21{:, 2};
ImS21 = table21{:, 3};
s_21 = ReS21 + 1i * ImS21;

ReS12 = table12{:, 2};
ImS12 = table12{:, 3};
s_12 = ReS12 + 1i * ImS12;

ReS22 = table22{:, 2};
ImS22 = table22{:, 3};
s_22 = ReS22 + 1i * ImS22;

% plot(frequency, real(s_12))
% hold on
% plot(frequency, imag(s_12))
% grid on

%% deembading

R1 = exp(-1i*L1*2*pi/c .* frequency .* (1 - (f_0 ./ frequency).^2).^(1/2));
R2 = exp(-1i*L2*2*pi/c .* frequency .* (1 - (f_0 ./ frequency).^2).^(1/2));

new_s11 = zeros(size(frequency));
new_s12 = zeros(size(frequency));
new_s21 = zeros(size(frequency));
new_s22 = zeros(size(frequency));

for v=1:size(frequency,1)
    s_matr = [[s_11(v) s_12(v)]; [s_21(v) s_22(v)]];
    t_matr = from_S_to_T(s_matr);

    T_a_inverse = [[conj(R1(v)) 0]; [0 R1(v)]];
    T_b_inverse = [[conj(R2(v)) 0]; [0 R2(v)]];

    new_t_matr = T_a_inverse * t_matr * T_b_inverse;

    new_s_matr = from_T_to_S(new_t_matr);

    new_s11(v) = new_s_matr(1, 1);
    new_s12(v) = new_s_matr(1, 2);
    new_s21(v) = new_s_matr(2, 1);
    new_s22(v) = new_s_matr(2, 2);
end


%% Plotting results to check deembadding
% plot(frequency, real(new_s11), 'DisplayName', 'Re', 'LineWidth', 1.2)
% hold on
% plot(frequency, imag(new_s11), 'DisplayName', 'Im', 'LineWidth', 1.2)
% legend('Interpreter','latex', 'FontSize', 13)
% grid on

%% NRW Method

X = (1 + new_s11.^2 - new_s21.^2) ./ (2 * new_s11);
Gamma = zeros(size(X));

for k=1:size(Gamma)
    a = X(k);
    b = (X(k).^2 - 1).^(1./2);
    if abs(a+b) <= 1
        Gamma(k) = a+b;
    else
        Gamma(k) = a-b;
    end
end

exp_gd = (1 - Gamma .* (new_s11 + new_s21)) ./ (new_s11 + new_s21 - Gamma);

epsmu = zeros(size(exp_gd)); 
phase_shift =  0; % сколько раз 2pi набежало в материале по фазе

phase_total = angle(exp_gd(1)) + phase_shift;
epsmu(1) = -(c./(2*pi*d.*frequency(1)).*(log(abs(exp_gd(1))) + ...
    + 1i*phase_total)).^2 + (f_0./frequency(1)).^2;


for k=2:size(epsmu)
    phase_total = phase_total + angle(exp_gd(k)/exp_gd(k-1));
    epsmu(k) = -(c./(2*pi*d.*frequency(k)).*(log(abs(exp_gd(k))) + ...
        1i*phase_total)).^2 + (f_0./frequency(k)).^2;
end


impedance = (((1+new_s11).^2 - new_s21.^2) ./ ((1-new_s11).^2 - new_s21.^2)).^(1/2);
mu_r = impedance .* (epsmu - (f_0./frequency(k)).^2).^(1/2) ./ ((1 - (f_0 ./ frequency).^2).^(1/2));

%% results plot

plot(frequency, real(mu_r), 'DisplayName', 'Re$(\mu)$', 'LineWidth', 1.2)
hold on
plot(frequency, imag(mu_r), 'DisplayName', 'Im$(\mu)$', 'LineWidth', 1.2)
plot(frequency, real(epsmu), 'DisplayName', 'Re$(\varepsilon)$', 'LineWidth', 1.2)
plot(frequency, imag(epsmu), 'DisplayName', 'Im$(\varepsilon)$', 'LineWidth', 1.2)
hold off
legend('Interpreter','latex', 'FontSize', 13)
grid on