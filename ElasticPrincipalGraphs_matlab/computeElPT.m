function [NodePositions,Edges,ReportTable,cpg] = computeElPT(data,NumNodes,varargin)

% Elasticity module stretching
EP = -1;
% Elasticity module bending
RP = -1;

parameterfunction_handle = @parametersDefaultPrincipalTree;


    for i=1:length(varargin)
        if strcmpi(varargin{i},'EP')
            EP = varargin{i+1};
        elseif strcmpi(varargin{i},'RP')
            RP = varargin{i+1};
        elseif strcmpi(varargin{i},'ParameterSet')
            parameterfunction_handle = varargin{i+1};
        end
    end

javaclasspath({'VDAOEngine.jar'});

parameters = parameterfunction_handle();

config = vdaoengine.analysis.grammars.ConfigFile;

config.algtype = parameters.algtype;
config.stretchInitCoeffs(1) = parameters.initstretch1;
config.initStrategy = parameters.initalgorithm;

for i=1:size(parameters.epochs,1)
    epoch = vdaoengine.analysis.elmap.ElmapAlgorithmEpoch;
    
    epoch.grammarType = parameters.epochs(i).grammartype;

    epoch.EP = parameters.epochs(i).ep;
    if EP(i)>0
        epoch.EP = EP(i);
    end
    epoch.RP = parameters.epochs(i).rp;
    if RP(i)>0
        epoch.RP = RP(i);
    end
    %epoch.numberOfIterations = parameters.epochs(i).numiter;
    epoch.numberOfIterations = NumNodes(i);
    epoch.maxNumberOfIterationsForSLAU = parameters.epochs(i).numiterSLAU;
    epoch.epsConvergence = parameters.epochs(i).eps;
    epoch.epsConvergenceSLAU = parameters.epochs(i).epsSLAU;
    
    if isfield(parameters.epochs(i),'robust')
        epoch.robust
        parameters.epochs(i).robust
        epoch.robust = parameters.epochs(i).robust;
    end
    if isfield(parameters.epochs(i),'trimradius')
        epoch.trimradius = parameters.epochs(i).trimradius;
    end
    
    
    config.epochs.add(epoch);
end

cpg = vdaoengine.analysis.grammars.ComputePrincipalGraph;

cpg.config = config;

cpg.setDataSetAsMassif(data);
report = cpg.compute();

fn = tempname;

fid = fopen(fn,'w');
fprintf(fid, '%s', char(report));
fclose(fid);

ReportTable = importdata(fn);

NodePositions = cpg.graph.getNodePositions();
Edges = cpg.graph.getEdgeTable();