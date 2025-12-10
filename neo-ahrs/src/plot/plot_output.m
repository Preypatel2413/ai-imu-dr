function plot_output(csv_name)
close all

    % csv_name = "02_09_2025_07_03_13_196460f1.neolog.csv";

output_file   = strcat("../data/output_", csv_name);
original_file = strcat("../../../data/", csv_name);

data_output = csvread(output_file);
data_output = data_output(2:end, :);

t  = data_output(:, get_csvcol(output_file, "t"));
% ax = data_output(:, get_csvcol(output_file, "ax"));
% ay = data_output(:, get_csvcol(output_file, "ay"));
% az = data_output(:, get_csvcol(output_file, "az"));
% gx = data_output(:, get_csvcol(output_file, "gx"));
% gy = data_output(:, get_csvcol(output_file, "gy"));
% gz = data_output(:, get_csvcol(output_file, "gz"));
% mx = data_output(:, get_csvcol(output_file, "mx"));
% my = data_output(:, get_csvcol(output_file, "my"));
% mz = data_output(:, get_csvcol(output_file, "mz"));
ym = data_output(:, get_csvcol(output_file, "ym"));
ys = data_output(:, get_csvcol(output_file, "ys"));
% ps = data_output(:, get_csvcol(output_file, "ps"));
% rs = data_output(:, get_csvcol(output_file, "rs"));
% alx = data_output(:, get_csvcol(output_file, "alx"));
% aly = data_output(:, get_csvcol(output_file, "aly"));
% alz = data_output(:, get_csvcol(output_file, "alz"));

data_original = csvread(original_file);
data_original = data_original(2:end, :);

tt =  data_original(:, get_csvcol(original_file, "Time(ms)"));
map = data_original(:, get_csvcol(original_file, "Map Bearing"));

map_tol = 100;
yaw_tol = 200;
map = unwrap(map, map_tol);
ys = unwrap(ys, map_tol);
ym = unwrap(ym, map_tol);

% Prepare map bearing profile for plotting.
tt = (tt - tt(1)) / 1000.0;
[tt_unique, idx_unique] = unique(tt, 'stable');
map = map(idx_unique);
valid_times = tt_unique > 0;
t_map = tt_unique(valid_times);
map = map(valid_times);
map = map - map(1);

fig_idx = 1;

% figure(fig_idx); fig_idx = fig_idx + 1;
% subplot(3,1,1); plot(t, ax); axis tight; title({"Accelerometer measurements", "ax"}); grid on; xlabel("Time (s)"); ylabel("m/s^2");
% subplot(3,1,2); plot(t, ay); axis tight; title("ay"); grid on; xlabel("Time (s)"); ylabel("m/s^2");
% subplot(3,1,3); plot(t, az); axis tight; title("az"); grid on; xlabel("Time (s)"); ylabel("m/s^2");

% figure(fig_idx); fig_idx = fig_idx + 1;
% subplot(3,1,1); plot(t, gx); axis tight; title({"Gyroscope measurements", "gx"}); grid on;  xlabel("Time (s)"); ylabel("deg/s");
% subplot(3,1,2); plot(t, gy); axis tight; title("gy"); grid on; xlabel("Time (s)"); ylabel("deg/s");
% subplot(3,1,3); plot(t, gz); axis tight; title("gz"); grid on; xlabel("Time (s)"); ylabel("deg/s");

% figure(fig_idx); fig_idx = fig_idx + 1;
% subplot(3,1,1); plot(t, mx); axis tight; title({"Magnetometer measurements", "mx"}); grid on;  xlabel("Time (s)"); ylabel("uT");
% subplot(3,1,2); plot(t, my); axis tight; title("my"); grid on; xlabel("Time (s)"); ylabel("uT");
% subplot(3,1,3); plot(t, mz); axis tight; title("mz"); grid on; xlabel("Time (s)"); ylabel("uT");

% figure(fig_idx); fig_idx = fig_idx + 1;
% subplot(3,1,1); plot(t, ys); axis tight; title({"Madgwick filter output", "ys"}); grid on; xlabel("Time (s)"); ylabel("deg");
% subplot(3,1,2); plot(t, ps); axis tight; title("ps"); grid on; xlabel("Time (s)"); ylabel("deg");
% subplot(3,1,3); plot(t, rs); axis tight; title("rs"); grid on; xlabel("Time (s)"); ylabel("deg");

% figure(fig_idx); fig_idx = fig_idx + 1;
% subplot(3,1,1); plot(t, alx); axis tight; title({"Linear acceleration", "alx"}); grid on; xlabel("Time (s)"); ylabel("deg");
% subplot(3,1,2); plot(t, aly); axis tight; title("aly"); grid on; xlabel("Time (s)"); ylabel("deg");
% subplot(3,1,3); plot(t, alz); axis tight; title("alz"); grid on; xlabel("Time (s)"); ylabel("deg");

figure(fig_idx); fig_idx = fig_idx + 1;
plot(t, ym); hold on;  plot(t, ys); plot(t_map, map);
title("Heading comparisons"); grid on; axis tight;
legend("ym (Magnetometer)", "ys (Accel-Gyro yaw)", "ymap (Map bearing)");
end