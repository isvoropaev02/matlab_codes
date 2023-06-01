%% geometry
d = 0.020; % размер образца
a = 0.02286; % размер большей стенки волновода
c = 3 * 10.^8; % скорость света
f_0 = c / (2 * a); % критическая частота низшей моды

L1 = 0.099;
L2 = 0.50 - d - L1;

%% data import
table = readtable("with_object.txt");
%table_free = readtable("without_object.txt");

frequency = table{:, 1};

ReS11 = table{:, 2};% ./ table_free{:, 2};
ImS11 = table{:, 3};% ./ table_free{:, 3};
s_11 = ReS11 + 1i * ImS11;

ReS21 = table{:, 4};% ./ table_free{:, 4};
ImS21 = table{:, 5};% ./ table_free{:, 5};
s_21 = ReS21 + 1i * ImS21;

ReS12 = table{:, 6};% ./ table_free{:, 6};
ImS12 = table{:, 7};% ./ table_free{:, 7};
s_12 = ReS12 + 1i * ImS12;

ReS22 = table{:, 8};% ./ table_free{:, 8};
ImS22 = table{:, 9};% ./ table_free{:, 9};
s_22 = ReS22 + 1i * ImS22;


plot(frequency, 20*log10(abs(s_11)), 'DisplayName', 'abs$(s1)$', 'LineWidth', 1.2)
% plot(frequency, angle(s_11)*180/pi, 'DisplayName', 'angle$(s11)$', 'LineWidth', 1.2)
% xlabel('$f$, Hz', 'Interpreter','latex', 'FontSize', 13)
ylabel('Absolute value, dB', 'FontSize', 14)
% hold on
% plot(frequency, unwrap(angle(s_11))*180/pi, 'DisplayName', 'angle$(s11)$', 'LineWidth', 1.2)
% plot(frequency, real(s_21), 'DisplayName', 'Re$(s22)$', 'LineWidth', 1.2)
% plot(frequency, imag(s_21), 'DisplayName', 'Im$(s22)$', 'LineWidth', 1.2)
% hold off
legend('Interpreter','latex', 'FontSize', 13)
grid on

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

% plot(frequency, 20*log10(abs(new_s11)), 'DisplayName', 'abs$(s11)$', 'LineWidth', 1.2)
% xlabel('$f$, Hz', 'Interpreter','latex', 'FontSize', 13)
% ylabel('Absolute value, dB', 'FontSize', 14)
% plot(frequency, real(new_s11), 'DisplayName', 'Re$(s11)$', 'LineWidth', 1.2)
% hold on
% plot(frequency, imag(new_s11), 'DisplayName', 'Im$(s11)$', 'LineWidth', 1.2)
% plot(frequency, real(new_s21), 'DisplayName', 'Re$(s21)$', 'LineWidth', 1.2)
% plot(frequency, imag(new_s21), 'DisplayName', 'Im$(s21)$', 'LineWidth', 1.2)
% hold off
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
pl = zeros(size(exp_gd));
for v=1:size(exp_gd)
    pl(v) = my_angle(exp_gd(v));
end

phases = unwrap(pl);

% plot(frequency, phases.*180/pi, 'DisplayName', '$(arg(e^{\gamma d}))$', 'LineWidth', 1.2)
% xlabel('$f$, Hz', 'Interpreter','latex', 'FontSize', 13)
% ylabel('Phase, deg', 'FontSize', 14)
% legend('Interpreter','latex', 'FontSize', 13)
% grid on

epsmu = zeros(size(exp_gd)); 
phase_shift =  0; % сколько раз 2pi набежало в материале по фазе

phase_total = my_angle(exp_gd(1)) + phase_shift;
epsmu(1) = -(c./(2*pi*d.*frequency(1)).*(log(abs(exp_gd(1))) + ...
    + 1i*phase_total)).^2 + (f_0./frequency(1)).^2;


for k=2:size(epsmu)
    phase_total = phases(k);
    %phase_total = phase_total + angle(exp_gd(k)/exp_gd(k-1));
    epsmu(k) = -(c./(2*pi*d.*frequency(k)).*(log(abs(exp_gd(k))) + ...
        1i*phase_total)).^2 + (f_0./frequency(k)).^2;
end


%impedance = (((1+new_s11).^2 - new_s21.^2) ./ ((1-new_s11).^2 - new_s21.^2)).^(1/2);
%mu_r = impedance .* (epsmu - (f_0./frequency(k)).^2).^(1/2) ./ ((1 - (f_0 ./ frequency).^2).^(1/2));

%% results plot

% plot(frequency, real(mu_r), 'DisplayName', 'Re$(\mu)$', 'LineWidth', 1.2)
% hold on
% plot(frequency, imag(mu_r), 'DisplayName', 'Im$(\mu)$', 'LineWidth', 1.2)
plot(frequency, real(epsmu), 'DisplayName', 'Re$(\varepsilon)$', 'LineWidth', 1.2)
hold on
plot(frequency, imag(epsmu), 'DisplayName', 'Im$(\varepsilon)$', 'LineWidth', 1.2)
hold off
legend('Interpreter','latex', 'FontSize', 13)
xlabel('Frequency, Hz', 'Interpreter','latex', 'FontSize', 13)
ylabel('Relative $\varepsilon$', 'Interpreter','latex', 'FontSize', 13)
grid on
