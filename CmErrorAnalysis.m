% ------------------------------------------------------------
% Script:  CmErrorAnalysis
% Summary: performs an analysis of propagated uncertainty of metabolite
%          concentration estimates made using absolute quantification from
%          MR spectroscopy datasets in the human brain.
%          The following wraps individual analysis scripts into a single files
%          which can be found in the directory ./analysis.
%          Each script renders a series of charts and figures with
%          relevant results and metadata, including:
%
%          Primary Analysis
%           (1) analytical approximation of propagated uncertainty by parameter
%           (2) analytical approximation of propagated uncertainty by metabolite
%           (3) Monte Carlo simulated propagated uncertainty by parameter
%           (4) Monte Carlo simulated propagated uncertainty by metabolite
%           (5) parameter influence on analytical solution of propagated error
%           (6) case study in uncertainty propagation, multiple sclerosis
%           (7) case study in uncertainty propagation, macromolecules
%
%          Secondary Analysis
%           (S1) finite difference approximation of partial derivatives
%           (S2) convergence curves for Monte Carlo simluations
%           (S3) sensitivity of propagated uncertainty to covariances
% 
% Inputs:  None. All script parameters are set in two JSON files:
%          ./data/CmErrorParams.json
%          ./data/CmErrorMetabs.json
%
% Usage:   >> CmErrorAnalysis
%
% Citation: 
%          Please reference the following paper if used or adapted,
%
%          Instrella R, Juchem C. Uncertainty propagation in absolute metabolite 
%          quantification for in vivo magnetic resonance spectroscopy of the 
%          human brain. Magn Reson Med. 2023 (under review)
%
% Corresponding Author:   Ronald Instrella, M.S.
% Principal Investigator: Christoph Juchem, Ph.D.
% 
% Columbia University, 2023
% ------------------------------------------------------------
close all; clear; clc;

%% (0) analysis set up
addpath(genpath('./analysis'));
addpath(genpath('./data'));
addpath(genpath('./functions'));
addpath(genpath('./results'));
addpath(genpath('./utils'));
config = CmUtils.load_dataset('CmErrorConfig.json');

%% (1)	analytical approximation of propagated uncertainty by parameter
CmErrorByParamPlots(config)

%% (2)  analytical approximation of propagated uncertainty by metabolite
CmErrorByMetabPlots(config)

%% (3)  Monte Carlo simulated propagated uncertainty by parameter
CmErrorMonteCarloByParam(config)

%% (4)  Monte Carlo simulated propagated uncertainty by metabolite
CmErrorMonteCarloByMetab(config)

%% (5) 	parameter influence on analytical solution of propagated error
CmErrorParamInfluence(config)

%% (6) 	case study in uncertainty propagation, multiple sclerosis
CmErrorCaseStudyMS(config)

%% (7) 	case study in uncertainty propagation, macromolecules
CmErrorCaseStudyMM(config)

%% (S1) finite difference approximation of partial derivatives
CmErrorFiniteDifference(config)

%% (S2) convergence curves for Monte Carlo simluations
CmErrorMonteCarloConvergence(config)

%% (S3) sensitivity of propagated uncertainty to covariances
CmErrorSensitivityPlots(config)
