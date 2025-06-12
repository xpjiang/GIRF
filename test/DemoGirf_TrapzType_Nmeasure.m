% Script to test GIRF code
%
% Authors:  Johanna Vannesjo (johanna.vannesjo@gmail.com),
%           Lars Kasper
% Copyright (C) 2014 IBT, University of Zurich and ETH Zurich,
%               2016 FMRIB centre, University of Oxford
%
% This file is part of a code package for GIRF computation and application. 
% The package is available under a BSD 3-clause license. Further info see:
% https://github.com/MRI-gradient/girf
%

clear all, clear classes
close all;

%% Compute input pulses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create input time vector
dtIn = 1e-6;
tIn = (0:dtIn:50e-3)';

% % Define input sweep
% sweepType = 'deltaAM';
% Tsweep = 50e-3;
% f1 = 0;
% f2 = 40e3;
% phi0 = 0;
% A = 1;
% tStart = 5e-3;
% inSweep = sweeps(tIn,Tsweep,f1,f2,phi0,A,tStart,sweepType);

% Define input blips
N = 20;
format long;     % 15位小数
for i=1:N
    A = 30-15.0/N*(i*1); % 30mT/m-15mT/m
    slope = A/(30 / 300e-6) ; % SR 100mT/m/ms (0-30mT/m 300us)
    plateau = 0;
    tStart = 5.1e-3 - slope;
    slope
    tStart
    inTrapz = trapezoid(tIn,tStart,A,slope,plateau);
    
    %in(:,1,1) = inTrapz;
    in(:,1,i) = inTrapz;
    %in(:,1,2) = inSweep;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simulate output pulses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize output
dtOut = dtIn;
tOut = tIn;
% 0 B0
% 1 x
% 2 y
% 3 z
% ..
outBasis = 1:16;%
out = zeros(length(tOut),length(outBasis),size(in,3));

% Simulate system low-pass filtering
selfOut = BW_filter(in,tIn,20e3,'rc',0.6);
crossOut = 0.1*BW_filter(in,tIn,10e3,'rc',1);
out(:,1,:) = crossOut; % B0
out(:,2,:) = selfOut; % x
%out(:,3,:) = crossOut; % y
%out(:,4,:) = crossOut; % z

% Add simulated noise to integrated waveform
outK = dtOut*cumsum(out);
outK = outK + 1e-8*randn(size(outK));
out = [zeros(1,length(outBasis),size(in,3)); diff(outK)/dtOut];

% ====== 新增：将out向后延迟80us ======
delayPoints = 80e-6/dtIn;
% 创建延迟后的输出数组
delayedOut = zeros(size(out));
% 检查是否有足够的点来延迟
if size(out, 1) > delayPoints
    % 前delayPoints个点置为0
    delayedOut(1:delayPoints, :, :) = 0;
    % 后面的点用原始信号填充（延迟）
    delayedOut(delayPoints+1:end, :, :) = out(1:end-delayPoints, :, :);
else
    warning('信号长度不足以实现100点延迟！');
    delayedOut = out; % 如果信号太短，保持原样
end

% 使用延迟后的信号
out = delayedOut;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute GIRF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Initialize GIRF object %%%%%%
inChannels = {'X'};
girfo = GirfProvider(tIn, in, tOut, out, inChannels, outBasis);
%girfo.Vis('in','t');
girfo.Vis('inout','t');
girfo.Vis('inout','f');
girfo.Vis('sensitivity');

%%%% Compute GIRF %%%%%%%%%%%%%%%%%%%%%%%
girfo.ComputeGirf();
girfo.Vis('GIRF','f');
girfo.Vis('GIRF','t');
girfo.Vis('GIRF','f',[],1);

%%%% Filter GIRF %%%%%%%%%%%%%%%%%%
girfo.WindowFreq(60e3,'rc',1/4);
girfo.Vis('GIRF','f');
girfo.Vis('GIRF','t');

girfo.VarSmoothFreq(30e3);
girfo.Vis('GIRF','f');
girfo.Vis('GIRF','t');

girfo.WindowTime(10e-3,'rc',0);
girfo.Vis('GIRF','f');
girfo.Vis('GIRF','t');

%%%% Return GirfEssential object %%%%
girfE = girfo.GetGirfEssential();

%%%% Test perform a prediction %%%%
girfA = GirfApplier(girfE);
[gOut, kOut, tK] = girfA.PredictGrad(tIn, in(:,:,1), tOut, inChannels, 'conv');
figure, plot(tIn, in(:,:,1), tOut, out(:,2,1), tOut, gOut(:,2))
% 添加图例（可选但推荐）
legend('输入信号x', '测量输出x', 'GIRF prediction输出x');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('girf was calculated')

