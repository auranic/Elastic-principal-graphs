function [GE,MetageneMatrix] = TF2GE(TF,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

NumberOfGenes = 1000;
MetageneMatrix = [];

    for i=1:2:length(varargin)
        if strcmpi(varargin{i},'MetaGeneMatrix')
            MetageneMatrix = varargin{i+1};
            NumberOfGenes = size(MetageneMatrix,1);
        end
    end
      
GE = TF*MetageneMatrix';

end

