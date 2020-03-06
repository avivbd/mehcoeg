function runAllTests()
%% Run all tests in runTests.m
clc

runtests('testNumJack.m');
runtests('testInverter.m');
runtests('testIncrementor.m');
runtests('testSolver.m');
% runtests('testEKF.m');
% runtests('testRTS.m');
% runtests('testENKF.m');

end