d = 0.002; % размер образца
a = 0.008636; % размер большей стенки волновода
c = 3 * 10.^8; % скорость света
f_0 = c / (2 * a); % критическая частота низшей моды

%% data import
table11 = readtable("s11_to_edit.txt");
table21 = readtable("s21_to_edit.txt");
table12 = readtable("s12_to_edit.txt");
table22 = readtable("s22_to_edit.txt");


frequency = table11{:, 1} .* 10.^9;

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

%% forming S-matrices (not deembedded)

s_matrices = zeros(2, 2, size(frequency, 1));

for v=1:size(frequency,1)
    s_matrices(:, :, v) = [[s_11(v) s_12(v)]; [s_21(v) s_22(v)]];
end


%% calculating deembedding T-matrices (обратные Т матрицы свободного волновода)

L1 = 0.012;
L2 = 0.012;

R1 = exp(-1i*L1*2*pi/c .* frequency .* (1 - (f_0 ./ frequency).^2).^(1/2));
R2 = exp(-1i*L2*2*pi/c .* frequency .* (1 - (f_0 ./ frequency).^2).^(1/2));

T_a_inverse = zeros(2, 2, size(frequency, 1));
T_b_inverse = zeros(2, 2, size(frequency, 1));

for v=1:size(frequency, 1)
    T_a_inverse(:, :, v) = [[conj(R1(v)) 0]; [0 R1(v)]];
    T_b_inverse(:, :, v) = [[conj(R2(v)) 0]; [0 R2(v)]];
end


%% deembedding

t_matrices = zeros(2, 2, size(s_matrices, 3));

for v=1:size(s_matrices, 3)
    t_matrices(: , :, v) = from_S_to_T(s_matrices(: , :, v));
end


new_t_matrices = zeros(2, 2, size(t_matrices, 3));

for v=1:size(s_matrices, 3)
    new_t_matrices(: , :, v) = T_a_inverse(:,:,v) * t_matrices(:,:,v) * T_b_inverse(:,:,v);
end


new_s_matrices = zeros(2, 2, size(new_t_matrices, 3));

for v=1:size(s_matrices, 3)
    new_s_matrices(: , :, v) = from_T_to_S(new_t_matrices(: , :, v));
end


%% Plotting results
new_s_21 = zeros(size(frequency));
for v=1:size(s_matrices, 3)
    new_s_21(v) = new_s_matrices(2 , 1, v);
end


plot(frequency, real(new_s_21))
hold on
plot(frequency, imag(new_s_21))
grid on
