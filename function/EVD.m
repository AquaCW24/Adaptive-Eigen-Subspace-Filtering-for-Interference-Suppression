function [Data_2D, Data_3D] = EVD(InputMatrix, RFI_Index)
% Adaptive Eigen-Subspace Filtering for Interference Suppression
% Inputs:
% - InputMatrix: 3-dim complex matrix (size: [range(280), snapshots(32), freq(560)])
% - RFI_index: The radio frequency interference
%
% Outputs:
% - Data_2D: 2-dim complex matrix (squeeze the dim of snapshot)
% - Data_3D: 3-dim complex matrix

Data = InputMatrix;
segment = 5;
for rfi = 1:length(RFI_Index)
    % Step 1: Construct Submatrices and Estimate Covariance
    SubMatrix_start = size(InputMatrix,1)/segment*2;
    SubMatrixSize = size(InputMatrix,1)/segment; 
    SubMatrixCount = size(InputMatrix,1) - SubMatrix_start - SubMatrixSize;
    CovMatrix = zeros(SubMatrixSize,SubMatrixSize,size(InputMatrix,2));
    SubMatrix = zeros(SubMatrixSize,SubMatrixCount);
    
    for snap = 1:size(InputMatrix,2) % snapshot
        for i = 1 : SubMatrixCount % range submatrix
            SubMatrix(:,i) = InputMatrix(SubMatrix_start+(i : SubMatrixSize+i-1), snap, RFI_Index(rfi)); 
        end
        CovMatrix(:,:,snap) = CovMatrix(:,:,snap) + (SubMatrix * SubMatrix') / SubMatrixCount;
    end
    
    % Average Covariance Matrix
    AvgCovMatrix = mean(CovMatrix, 3);
    
    % Step 2: Eigenvalue Decomposition
    [EigenVectors, EigenValues] = eig(AvgCovMatrix);
    EigenValues = diag(EigenValues);
    [SortedEigenValues, idx] = sort(EigenValues, 'descend');
    EigenVectors = EigenVectors(:, idx);
    
    % Step 3: Mode Criterion for Interference Subspace
    LogEigenValues = floor(log10(SortedEigenValues));
    Mode = mode(LogEigenValues);
    r = find(LogEigenValues > Mode, 1, 'last');
    
    InterferenceSubspace = EigenVectors(:, 1:r);
    ProjectionMatrix = InterferenceSubspace * InterferenceSubspace';
    
    % Step 4: Interference Suppression
    for snap = 1:size(InputMatrix,2)
        for i = 1:segment
            Data(SubMatrixSize*(i-1)+1:SubMatrixSize*i, snap, RFI_Index(rfi)) = InputMatrix(SubMatrixSize*(i-1)+1:SubMatrixSize*i, snap, RFI_Index(rfi)) - ProjectionMatrix * InputMatrix(SubMatrixSize*(i-1)+1:SubMatrixSize*i, snap, RFI_Index(rfi));
        end
    end
    
end

Data_3D = Data;
Data_2D = (squeeze(sum(Data,2)));

end
