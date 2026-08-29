% RFI Detection and Mitigation
% last edited: 2026/04/16

% Input:
%       - InputMatrix: 3-D complex matrix (range x snapshot x frequency)
%       - MethodName: Selected algorithm ("RFIM", "FACS" or "EVD")
% Output:
%       - Output_2D: 2-D amplitude matrix

classdef RFIDM
    
    properties
        % Detection params
        TrainingCell (1,1) double = 4;
        GuardCell (1,1) double = 1;
        PFA (1,1) double = 8e-2;
        CFARMethod (1,1) string = "CA";
        
        % Mitigation algorithm params
        ReferenceGateNumber (1,1) double = 40; % for FACS
        PulseWidth (1,1) double = 8192e-6;     % for RFIM
        Segment (1,1) double = 5;              % for EVD

    end
    
    methods
        % Constructor: Initializes parameters when an object is created
        function obj = RFIDM(varargin)
            if nargin > 0
                % Set the other parameters in order...
            end
        end
        
        % ================= Main Control =================
        function [Output_2D, RFI_Index] = process(obj, InputMatrix, MethodName)
            % process - Detection and Mitigation Process
            
            % 1. Detection by CFAR
            [RFI_Index, ~] = obj.detect(InputMatrix);
            
            if isempty(RFI_Index)
                disp('No RFI detected; returning original data.');
                Output_2D = abs(squeeze(sum(InputMatrix, 2)));
                return;
            end
            
            % 2. RFI Mitigation using the selected algorithm
            switch upper(MethodName)
                case 'RFIM'
                    Output_2D = obj.runRFIM(InputMatrix, RFI_Index);
                    
                case 'FACS'
                    Output_2D = obj.runFACS(InputMatrix);
                    
                case 'EVD'
                    Output_2D = obj.runEVD(InputMatrix, RFI_Index);
                    
                otherwise
                    error('Unknown algorithm, choose "RFIM", "FACS" or "EVD".');
            end
        end

        % ================= CFAR Detection =================
        function [INDEX, THRESHOLD_POWER] = detect(obj, InputMatrix)
            % Trnasform 3D complex matrix to 2D Power matrix
            PM = abs(squeeze(sum(InputMatrix, 2))).^2; 

            % Use the mean of the upper range as the detection threshold
            POWER = mean(PM(floor(size(PM,1)/4)+1:size(PM,1), :), 1);
            
            % Calculate the standard deviation on the decibel scale
            STD = std(10*log10(PM), 0, 1);
            THRESHOLD_STD = mean(STD) - 2*std(STD);

            N = length(POWER); 
            alpha = 2 * obj.TrainingCell * (obj.PFA^(-1/(2*obj.TrainingCell)) - 1); 
            THRESHOLD_POWER = zeros(1, N);

            % Sliding Window Processing Workflow
            for i = 1:N
                % Define the boundaries
                left_start = max(1, i - obj.TrainingCell - obj.GuardCell);
                left_end = max(1, i - obj.GuardCell - 1);
                right_start = min(N, i + obj.GuardCell + 1);
                right_end = min(N, i + obj.TrainingCell + obj.GuardCell);

                % Handling boundary conditions
                left_window = POWER(left_start:left_end);
                right_window = POWER(right_start:right_end);
                if i <= (obj.GuardCell + 1), left_window = []; end
                if i >= (N - obj.GuardCell - 1), right_window = []; end
                
                % Estimate the noise power based on the specified method
                switch obj.CFARMethod
                    case "CA"
                        noise_power = (sum(left_window) + sum(right_window)) / ...
                                    (length(left_window) + length(right_window));
                    case "GOCA"
                        noise_power = max(mean(left_window), mean(right_window));
                    case "SOCA"
                        noise_power = min(mean(left_window), mean(right_window));
                    case "OS"
                        combined = [left_window, right_window];
                        sorted = sort(combined);
                        noise_power = sorted(floor(0.75 * length(combined)));
                    otherwise
                        noise_power = mean([left_window, right_window]);
                end

                % Calculate threshold
                THRESHOLD_POWER(i) = alpha * noise_power;
            end

            % Get index RFI appeared
            INDEX = find(POWER > THRESHOLD_POWER | STD < THRESHOLD_STD);
        end

        % ================= FACS  =================
        function Output = runFACS(obj, InputMatrix)
            % Divide the number of reference distance gates in half
            RefGate = floor(obj.ReferenceGateNumber / 2);
            Matrix = squeeze(sum(InputMatrix, 2));
            M = size(Matrix, 1); % Range
            % N = size(Matrix, 2); % Frequency

            Output = Matrix;
            for i = 1:M
                below_start = max(1, i - RefGate);
                below_end = max(1, i - 1);
                upper_start = min(M, i + 1);
                upper_end = min(M, i + RefGate);

                below_window = Matrix(below_start:below_end, :);
                upper_window = Matrix(upper_start:upper_end, :);

                if i == 1
                    below_window = [];
                elseif i == M
                    upper_window = [];
                end
                Reference_cell = [below_window; upper_window];

                % Weight Optimization Calculation
                Weight_Optimized = (Reference_cell * Reference_cell') \ (Reference_cell * Matrix(i,:)');
                Weight_Optimized = [1; -Weight_Optimized];
                Weight_Normalized = Weight_Optimized / sqrt(Weight_Optimized' * Weight_Optimized);
                
                Output(i,:) = Weight_Normalized' * [Matrix(i,:); Reference_cell];
            end
            % Transform to 2D amplitude matrix
            Output = abs(Output);
        end
        % ================= RFIM =================
        function Output = runRFIM(obj, InputMatrix, RFI_idx)
            fTs = obj.PulseWidth; 
            N = size(InputMatrix, 2); % snapshot number
            fT = N * fTs;
            fdFreq = 1 / fT; 
            fNFreq = 1 / (2 * fTs); 
            Freq = -fNFreq : fdFreq : fNFreq-fdFreq;

            Data = InputMatrix;

            for i = 1:length(RFI_idx)
                RFI_Map = Data(:, :, RFI_idx(i));
                RFI_Map_FFT = abs(fftshift(fft(RFI_Map, N, 2) ./ N));
                
                for j = 1:size(RFI_Map, 1)
                    [~, idx] = max(RFI_Map_FFT(j, :));
                    
                    if idx == 1
                        RFI_DopplerFreq = Freq(idx) + fdFreq*(RFI_Map_FFT(j,idx+1)/(RFI_Map_FFT(j,idx)+RFI_Map_FFT(j,idx+1)));
                    elseif idx == N
                        RFI_DopplerFreq = Freq(idx) - fdFreq*(RFI_Map_FFT(j,idx-1)/(RFI_Map_FFT(j,idx)+RFI_Map_FFT(j,idx-1)));
                    else
                        if RFI_Map_FFT(j,idx-1) > RFI_Map_FFT(j,idx+1)
                            RFI_DopplerFreq = Freq(idx) - fdFreq*(RFI_Map_FFT(j,idx-1)/(RFI_Map_FFT(j,idx)+RFI_Map_FFT(j,idx-1)));
                        else
                            RFI_DopplerFreq = Freq(idx) + fdFreq*(RFI_Map_FFT(j,idx+1)/(RFI_Map_FFT(j,idx)+RFI_Map_FFT(j,idx+1)));
                        end
                    end

                    RFI_Amplitude = 0;
                    for k = 1:N
                        RFI_Amplitude = RFI_Amplitude + RFI_Map(j,k) * exp(-1i*2*pi*RFI_DopplerFreq*k*fTs);
                    end
                    RFI_Amplitude = RFI_Amplitude / N;
                    
                    RFI_AmplitudeInverse = RFI_Amplitude * exp(1i*2*pi*RFI_DopplerFreq*(1:N)*fTs);
                    Data(j, :, RFI_idx(i)) = Data(j, :, RFI_idx(i)) - RFI_AmplitudeInverse;
                end
            end
            % Transform to 2D amplitude matrix
            Output = abs(squeeze(sum(Data, 2))); 
        end

        % ================= EVD =================
        function Output = runEVD(obj, InputMatrix, RFI_Idx)
            Input_CP = InputMatrix;
            segment = obj.Segment;
            for rfi = 1:length(RFI_Idx)
                % Construct Submatrices and Estimate Covariance
                SSI = size(InputMatrix,1)/segment*2; % SubMatrix start index
                SWS = size(InputMatrix,1)/segment; % SubMatrix Window Size
                SSC = size(InputMatrix,1) - SSI - SWS + 1; % SubMatrix Slide Count
                CovMatrix = zeros(SWS,SWS,size(InputMatrix,2));
                SubMatrix = zeros(SWS,SSC);
                
                for snap = 1:size(InputMatrix,2) % snapshot
                    for i = 1 : SSC % range submatrix
                        SubMatrix(:,i) = InputMatrix(SSI+(i : SWS+i-1), snap, RFI_Idx(rfi)); 
                    end
                    CovMatrix(:,:,snap) = (SubMatrix * SubMatrix') / SSC;
                end
                
                % Average Covariance Matrix
                AvgCovMatrix = mean(CovMatrix, 3);
                
                % Eigenvalue Decomposition
                [EigenVectors, EigenValues] = eig(AvgCovMatrix);
                EigenValues = diag(EigenValues);
                [SortedEigenValues, idx] = sort(EigenValues, 'descend');
                EigenVectors = EigenVectors(:, idx);
                
                % Mode Criterion for Interference Subspace
                LogEigenValues = floor(log10(SortedEigenValues));
                Mode = mode(LogEigenValues);
                r = find(LogEigenValues > Mode, 1, 'last');
                
                InterferenceSubspace = EigenVectors(:, 1:r);
                % obtain projection matrix
                PM = InterferenceSubspace * InterferenceSubspace'; 
                
                % Interference Suppression
                for snap = 1:size(InputMatrix,2)
                    for i = 1:segment
                        Input_CP(SWS*(i-1)+1:SWS*i, snap, RFI_Idx(rfi)) = ...
                            (eye(size(PM,1)) - PM) * InputMatrix(SWS*(i-1)+1:SWS*i, snap, RFI_Idx(rfi));
                    end
                end
            end
            % Transform to 2D amplitude matrix
            Output = abs(squeeze(sum(Input_CP,2)));
        end
    end
    
    methods (Access = private)
        % Private methods: These encapsulate the original functions, 
        % so these algorithms cannot be called directly from outside.
    end
end