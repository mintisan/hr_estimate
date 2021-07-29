from KalmanNet_data import F_mat, H_mat, T, m1x_0, m2x_0
from KalmanNet_data import N_T

from KalmanFilter_test import KFTest
from KalmanNet_test import NNTest
from KalmanNet_plt import KFPlot
from KalmanNet_sysmdl import SystemModel

from datetime import datetime

print("Kalman Net Start")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

T = 10

q_base = 2
r_base = 2

SNR_lin = [1/4, 1/2, 1, 2, 4]
SNR_dB = [-6, -3, 0, 3, 6]

res_grid = [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]]

# Initialize Test Model
SysModel_test = SystemModel(F_mat, q_base, H_mat, r_base, 10)
SysModel_test.InitSequence(m1x_0, m2x_0)

# Set Test Model
SysModel_test.UpdateCovariance(q_base * SNR_lin[0], r_base * SNR_lin[2])

# Generate Data
SysModel_test.GenerateBatch(N_T, 0)
test_input_m = SysModel_test.Input
test_target_m = SysModel_test.Target

# Set Test Model
SysModel_test.UpdateCovariance(q_base * SNR_lin[2], r_base * SNR_lin[2])

# Generate Data
SysModel_test.GenerateBatch(N_T, 0)
test_input_b = SysModel_test.Input
test_target_b = SysModel_test.Target

# Set Test Model
SysModel_test.UpdateCovariance(q_base * SNR_lin[4], r_base * SNR_lin[2])

# Generate Data
SysModel_test.GenerateBatch(N_T, 0)
test_input_p = SysModel_test.Input
test_target_p = SysModel_test.Target

# Set Test Model
SysModel_test.UpdateCovariance(q_base * SNR_lin[0], r_base * SNR_lin[2])

# Test Theoretical Kalman Filter
[MSE_KF_linear_arr, MSE_KF_linear_avg, MSE_KF_dB_avg] = KFTest(SysModel_test, test_input_m, test_target_m)
res_grid[0][0] = MSE_KF_dB_avg ;

# Test Theoretical Kalman Filter
[MSE_KF_linear_arr, MSE_KF_linear_avg, MSE_KF_dB_avg] = KFTest(SysModel_test, test_input_b, test_target_b)
res_grid[0][1] = MSE_KF_dB_avg ;

# Test Theoretical Kalman Filter
[MSE_KF_linear_arr, MSE_KF_linear_avg, MSE_KF_dB_avg] = KFTest(SysModel_test, test_input_p, test_target_p)
res_grid[0][2] = MSE_KF_dB_avg ;

# Set Test Model
SysModel_test.UpdateCovariance(q_base * SNR_lin[2], r_base * SNR_lin[2])

# Test Theoretical Kalman Filter
[MSE_KF_linear_arr, MSE_KF_linear_avg, MSE_KF_dB_avg] = KFTest(SysModel_test, test_input_m, test_target_m)
res_grid[1][0] = MSE_KF_dB_avg ;

# Test Theoretical Kalman Filter
[MSE_KF_linear_arr, MSE_KF_linear_avg, MSE_KF_dB_avg] = KFTest(SysModel_test, test_input_b, test_target_b)
res_grid[1][1] = MSE_KF_dB_avg ;

# Test Theoretical Kalman Filter
[MSE_KF_linear_arr, MSE_KF_linear_avg, MSE_KF_dB_avg] = KFTest(SysModel_test, test_input_p, test_target_p)
res_grid[1][2] = MSE_KF_dB_avg ;


# Test Theoretical Kalman Filter
[MSE_KF_linear_arr, MSE_KF_linear_avg, MSE_KF_dB_avg] = NNTest(SysModel_test, test_input_m, test_target_m)
res_grid[3][0] = MSE_KF_dB_avg ;

# Test Theoretical Kalman Filter
[MSE_KF_linear_arr, MSE_KF_linear_avg, MSE_KF_dB_avg] = NNTest(SysModel_test, test_input_b, test_target_b)
res_grid[3][1] = MSE_KF_dB_avg ;

# Test Theoretical Kalman Filter
[MSE_KF_linear_arr, MSE_KF_linear_avg, MSE_KF_dB_avg] = NNTest(SysModel_test, test_input_p, test_target_p)
res_grid[3][2] = MSE_KF_dB_avg ;

# Set Test Model
SysModel_test.UpdateCovariance(q_base * SNR_lin[4], r_base * SNR_lin[2])

# Test Theoretical Kalman Filter
[MSE_KF_linear_arr, MSE_KF_linear_avg, MSE_KF_dB_avg] = KFTest(SysModel_test, test_input_m, test_target_m)
res_grid[2][0] = MSE_KF_dB_avg ;

# Test Theoretical Kalman Filter
[MSE_KF_linear_arr, MSE_KF_linear_avg, MSE_KF_dB_avg] = KFTest(SysModel_test, test_input_b, test_target_b)
res_grid[2][1] = MSE_KF_dB_avg ;

# Test Theoretical Kalman Filter
[MSE_KF_linear_arr, MSE_KF_linear_avg, MSE_KF_dB_avg] = KFTest(SysModel_test, test_input_p, test_target_p)
res_grid[2][2] = MSE_KF_dB_avg ;


KFPlot(res_grid)


# Plot
#NNPlot_test(MSE_KF_linear_arr, MSE_KF_linear_avg, MSE_KF_dB_avg,
#       MSE_test_linear_arr, MSE_test_linear_avg, MSE_test_dB_avg)
