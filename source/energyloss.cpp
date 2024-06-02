#include "energyloss.hpp"
#include "grids.hpp"
#include "linearinterpolation.hpp"
#include "polyintegration.hpp"

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <random>
#include <map>
#include <tuple>
#include <fstream>
#include <cmath>
#include <iomanip>

energyLoss::energyLoss(int argc, const char *argv[])
{
	m_error = false;

	std::vector<std::string> inputs; for (int i=2; i<argc; i++) inputs.push_back(argv[i]);

	if ((inputs.size() == 1) && (inputs[0] == "-h")) {
		std::cout << "default values: --collsys=PbPb --sNN=5020GeV --pName=Charm --centrality=30-40% --xB=0.6 --eventN=1000 --BCPP=20% --phiGridN=25 --TIMESTEP=0.1 --TCRIT=0.155 --BCPSEED=0" << std::endl;
		m_error = true;
	}

	std::map<std::string, std::string> inputParams;
	for (const auto &in : inputs)
	{
 	   	std::string key = in.substr(0, in.find("="));
 	   	std::string::size_type n = 0; while ((n = key.find("-", n)) != std::string::npos) {key.replace(n, 1, ""); n += 0;} //replacing all '-'
		std::string val = in.substr(in.find("=")+1, in.length());
		inputParams[key] = val;
	}
	std::vector<std::string> arguments = {"collsys", "sNN", "pName", "centrality", "xB", "eventN", "BCPP", "phiGridN", "TIMESTEP", "TCRIT", "BCPSEED", "config", "h"};
	for (const auto &inputParam : inputParams) {
		if(std::find(arguments.begin(), arguments.end(), inputParam.first) == arguments.end()) {
			std::cerr << "Error: provided argument flag: '" << inputParam.first << "' is not an option." << std::endl;
			std::cerr << "Valid parameters and default values are: ";
			std::cerr << "--collsys=PbPb --sNN=5020GeV --pName=Charm --centrality=30-40% --xB=0.6 --eventN=1000 --BCPP=20% --phiGridN=25 --TIMESTEP=0.1 --TCRIT=0.155 --BCPSEED=0" << std::endl;
			std::cerr << "For congiguration file use: --config=[pathToConfFile]" << std::endl;
			m_error = true;
		}
	}

	//checking if configuration file is provided:
	std::map<std::string, std::string> inputParamsFile;
	if (inputParams.count("config") > 0) {
		if (loadInputsFromFile(inputParams.at("config"), inputParamsFile) != 1) {
			m_error = true;
		}
	}
	std::vector<std::string> argumentsFile = {"collsys", "sNN", "pName", "centrality", "xB", "eventN", "BCPP", "phiGridN", "TIMESTEP", "TCRIT", "BCPSEED"};
	for (const auto &inputParam : inputParamsFile) {
		if(std::find(argumentsFile.begin(), argumentsFile.end(), inputParam.first) == argumentsFile.end()) {
			std::cerr << "Error: in configration file provided argument: '" << inputParam.first << "' is not an option." << std::endl;
			std::cerr << "Valid parameters and default values are: \n";
			std::cerr << "collsys = PbPb\nsNN = 5020GeV\npName = Charm\ncentrality = 30-40%\nxB = 0.6\neventN = 1000\nBCPP = 20%\nphiGridN = 25\nTIMESTEP = 0.1\nTCRIT = 0.155\nBCPSEED = 0" << std::endl;
			m_error = true;
		}
	}

	//setting parameter values based on config file values and overwriting with command line values:
	//
	m_collsys = "PbPb"; if (inputParamsFile.count("collsys") > 0) m_collsys = inputParamsFile["collsys"];
						if (    inputParams.count("collsys") > 0) m_collsys =     inputParams["collsys"];
	
	m_sNN = "5020GeV"; if (inputParamsFile.count("sNN") > 0) m_sNN = inputParamsFile["sNN"];
					   if (    inputParams.count("sNN") > 0) m_sNN =     inputParams["sNN"];

	m_pName = "Charm"; if (inputParamsFile.count("pName") > 0) m_pName = inputParamsFile["pName"];
					   if (    inputParams.count("pName") > 0) m_pName =     inputParams["pName"];

	m_centrality = "30-40%"; if (inputParamsFile.count("centrality") > 0) m_centrality = inputParamsFile["centrality"];
						     if (    inputParams.count("centrality") > 0) m_centrality =     inputParams["centrality"];

	m_xB = 0.6; if (inputParamsFile.count("xB") > 0) m_xB = stod(inputParamsFile["xB"]);
				if (    inputParams.count("xB") > 0) m_xB = stod(    inputParams["xB"]);

	m_eventN = 1000; if (inputParamsFile.count("eventN") > 0) m_eventN = stoi(inputParamsFile["eventN"]);
					 if (    inputParams.count("eventN") > 0) m_eventN = stoi(    inputParams["eventN"]);

	std::string bcppstr = "20%"; if (inputParamsFile.count("BCPP") > 0) bcppstr = inputParamsFile["BCPP"];
						         if (    inputParams.count("BCPP") > 0) bcppstr =     inputParams["BCPP"];
	bcppstr.replace(bcppstr.find("%"), 1, ""); m_BCPP = stod(bcppstr)/100.0;

	m_phiGridN = 25; if (inputParamsFile.count("phiGridN") > 0) m_phiGridN = stoi(inputParamsFile["phiGridN"]);
					 if (    inputParams.count("phiGridN") > 0) m_phiGridN = stoi(    inputParams["phiGridN"]);

	m_TIMESTEP = 0.1; if (inputParamsFile.count("TIMESTEP") > 0) m_TIMESTEP = stod(inputParamsFile["TIMESTEP"]);
					  if (    inputParams.count("TIMESTEP") > 0) m_TIMESTEP = stod(    inputParams["TIMESTEP"]);

	m_TCRIT = 0.155; if (inputParamsFile.count("TCRIT") > 0) m_TCRIT = stod(inputParamsFile["TCRIT"]);
					 if (    inputParams.count("TCRIT") > 0) m_TCRIT = stod(    inputParams["TCRIT"]);

	m_BCPSEED = 0; if (inputParamsFile.count("BCPSEED") > 0) m_BCPSEED = stoi(inputParamsFile["BCPSEED"]);
				   if (    inputParams.count("BCPSEED") > 0) m_BCPSEED = stoi(    inputParams["BCPSEED"]);

	//checking if provided value of sNN is an option:
	if ((m_sNN != "5440GeV") && (m_sNN != "5020GeV") && (m_sNN != "2760GeV") && (m_sNN != "200GeV")) {
		std::cerr << "Error: provided sNN parameter not an option, please try 5440GeV, 5020GeV, 2760GeV or 200GeV. Aborting..." << std::endl;
		m_error = true;
	}

	m_nf = m_sNN == "200GeV" ? 2.5 : 3.0;
	double T = 3.0 / 2.0*m_TCRIT;
	double mu = 0.197*std::sqrt((-8.0*(6.0+m_nf)*M_PI*M_PI*T*T)/(2.0*m_nf-33.0)/m_lambda/m_lambda/productLog((-8.0*(6.0+m_nf)*M_PI*M_PI*T*T)/(2.0*m_nf-33.0)/m_lambda/m_lambda));
	m_mgC = mu / std::sqrt(2.0);
	if (m_pName == "Bottom") m_MC = 4.75;
	else if (m_pName == "Charm") m_MC = 1.2;
	else if (m_pName == "Gluon") m_MC = mu/std::sqrt(2.0);
	else m_MC = mu/sqrt(6.0);
	m_TCollConst = T;
}

int energyLoss::loadInputsFromFile(const std::string &filePath, std::map<std::string, std::string> &inputParamsFile)
{
	std::ifstream file_in(filePath);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open configuration file. Aborting..." << std::endl;
		return -1;
	}
	std::string line, key, sep, val;
	while (std::getline(file_in, line))
	{
		std::stringstream ss(line);
		ss >> key; ss >> sep; ss >> val;
		inputParamsFile[key] = val;
	}
	file_in.close();
	return 1;
}

energyLoss::~energyLoss() {}

void energyLoss::runEnergyLoss()
{
	if (m_error) return;

	m_Grids.setGridPoints(m_sNN, m_pName, m_TCRIT);

	if (loadLdndx() != 1) return;
	if (loadLNorm() != 1) return;
	if (loadLColl() != 1) return;

	if (generateTempGrid() != 1) return;
	
 	if (loadPhiPoints() != 1) return;

	if ((m_pName == "Bottom") || (m_pName == "Charm")) {
		runELossHeavyFlavour();
	}
	else if (m_pName == "LQuarks") {
		runELossLightQuarks();
	}
	else {
		runELossLightFlavour();
	}
}

double energyLoss::productLog(double x) const
{
	if (x == 0.0) {
		return 0.0;
	}

	double w0, w1;
	if (x > 0.0) {
		w0 = std::log(1.2 * x / std::log(2.4 * x / std::log1p(2.4 * x)));
	}
	else {
		double v = 1.4142135623730950488 * std::sqrt(1.0 + 2.7182818284590452354 * x);
		double N2 = 10.242640687119285146 + 1.9797586132081854940 * v;
		double N1 = 0.29289321881345247560 * (1.4142135623730950488 + N2);
		w0 = -1 + v * (N2 + v) / (N2 + v + N1 * v);
	}

	while (true) {
		double e = std::exp(w0);
		double f = w0 * e - x;
		w1 = w0 - f / ((e * (w0 + 1.0) - (w0 + 2.0) * f / (w0 + w0 + 2.0)));
		if (std::abs(w0 / w1 - 1.0) < 1.4901161193847656e-8) {
			break;
		}
		w0 = w1;
	}
	return w1;
}


int energyLoss::loaddsdpti2(const std::string &pname, interpolationF<double> &dsdpti2int)const 
{
	const std::string path_in = "./ptDists/ptDist" + m_sNN + "/ptDist_" + m_sNN + "_" + pname + ".dat";

	std::ifstream file_in(path_in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open initial pT distribution file." << std::endl;
		return -1;
	}

	std::vector<double> pTdistX, pTdistF;

	std::string line; double buffer;

	while (std::getline(file_in, line))
	{
        if (line.at(0) == '#')
            continue;

		std::stringstream ss(line);
		ss >> buffer; pTdistX.push_back(buffer);
		ss >> buffer; pTdistF.push_back(buffer);
	}

	dsdpti2int.setData(pTdistX, pTdistF);

	file_in.close();

	return 1;
}

int energyLoss::loadLdndx()
{
	std::string partName;
	if (m_pName == "Bottom") partName = "Bottom";
	else if (m_pName == "Charm") partName = "Charm";
	else if (m_pName == "Gluon") partName = "Gluon";
	else partName = "LQuarks";

	std::stringstream xBss; xBss << std::fixed << std::setprecision(1) << m_xB;
	std::stringstream nfss; nfss << std::fixed << std::setprecision(1) << m_nf;

	const std::string path_in = "./ltables/ldndx_nf=" + nfss.str() + "_" + partName + "_xB=" + xBss.str() + ".dat";

	std::ifstream file_in(path_in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open Ldndx table file." << std::endl;
		return -1;
	}

	std::vector<double> Ldndx_tau, Ldndx_p, Ldndx_T, Ldndx_x, Ldndx_f;

	std::string line; double buffer;

	while (std::getline(file_in, line))
	{
        if (line.at(0) == '#')
            continue;

		std::stringstream ss(line);
		ss >> buffer; Ldndx_tau.push_back(buffer);
		ss >> buffer; Ldndx_p.push_back(buffer);
		ss >> buffer; Ldndx_T.push_back(buffer);
		ss >> buffer; Ldndx_x.push_back(buffer);
		ss >> buffer; Ldndx_f.push_back(buffer);
	}

	file_in.close();

	m_Ldndx.setData(Ldndx_tau, Ldndx_p, Ldndx_T, Ldndx_x, Ldndx_f);

	std::vector<std::vector<double>> domain = m_Ldndx.domain();
	if (m_Grids.tauPts(0)  < domain[0][0]) {std::cerr << "Error: tau grid point(s) out of lower bound of Ldndx domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.tauPts(-1) > domain[0][1]) {std::cerr << "Error: tau grid point(s) out of upper bound of Ldndx domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.pPts(0)    < domain[1][0]) {std::cerr << "Error:   p grid point(s) out of lower bound of Ldndx domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.pPts(-1)   > domain[1][1]) {std::cerr << "Error:   p grid point(s) out of upper bound of Ldndx domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.TPts(0)    < domain[2][0]) {std::cerr << "Error:   T grid point(s) out of lower bound of Ldndx domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.TPts(-1)   > domain[2][1]) {std::cerr << "Error:   T grid point(s) out of upper bound of Ldndx domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.xPts(0)    < domain[3][0]) {std::cerr << "Error:   x grid point(s) out of lower bound of Ldndx domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.xPts(-1)   > domain[3][1]) {std::cerr << "Error:   x grid point(s) out of upper bound of Ldndx domain. Aborting..." << std::endl; return -1;}

	return 1;
}

int energyLoss::loadLNorm()
{
	std::string partName;
	if (m_pName == "Bottom") partName = "Bottom";
	else if (m_pName == "Charm") partName = "Charm";
	else if (m_pName == "Gluon") partName = "Gluon";
	else partName = "LQuarks";

	std::stringstream xBss; xBss << std::fixed << std::setprecision(1) << m_xB;
	std::stringstream nfss; nfss << std::fixed << std::setprecision(1) << m_nf;

	const std::string path_in = "./ltables/lnorm_nf=" + nfss.str() + "_" + partName + "_xB=" + xBss.str() + ".dat";

	std::ifstream file_in(path_in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open LNorm table file." << std::endl;
		return -1;
	}

	std::vector<double> LNorm_tau, LNorm_p, LNorm_T, LNorm_f; //defining vectors that store LNorm table values

	std::string line; double buffer;

	while (std::getline(file_in, line))
	{
        if (line.at(0) == '#')
            continue;

		std::stringstream ss(line);
		ss >> buffer; LNorm_tau.push_back(buffer);
		ss >> buffer; LNorm_p.push_back(buffer);
		ss >> buffer; LNorm_T.push_back(buffer);
		ss >> buffer; LNorm_f.push_back(buffer);
	}

	file_in.close();

	m_LNorm.setData(LNorm_tau, LNorm_p, LNorm_T, LNorm_f);

	std::vector<std::vector<double>> domain = m_LNorm.domain();
	if (m_Grids.tauPts(0)  < domain[0][0]) {std::cerr << "Error: tau grid point(s) out of lower bound of LNorm domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.tauPts(-1) > domain[0][1]) {std::cerr << "Error: tau grid point(s) out of upeer bound of LNorm domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.pPts(0)    < domain[1][0]) {std::cerr << "Error:   p grid point(s) out of lower bound of LNorm domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.pPts(-1)   > domain[1][1]) {std::cerr << "Error:   p grid point(s) out of upeer bound of LNorm domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.TPts(0)    < domain[2][0]) {std::cerr << "Error:   T grid point(s) out of lower bound of LNorm domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.TPts(-1)   > domain[2][1]) {std::cerr << "Error:   T grid point(s) out of upeer bound of LNorm domain. Aborting..." << std::endl; return -1;}

	return 1;
}

int energyLoss::loadLColl()
{
	std::string partName;
	if (m_pName == "Bottom") partName = "Bottom";
	else if (m_pName == "Charm") partName = "Charm";
	else if (m_pName == "Gluon") partName = "Gluon";
	else partName = "LQuarks";

	std::stringstream nfss; nfss << std::fixed << std::setprecision(1) << m_nf;

	const std::string path_in = "./ltables/lcoll_nf=" + nfss.str() + "_" + partName + ".dat";

	std::ifstream file_in(path_in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open LColl table file." << std::endl;
		return -1;
	}

	std::vector<double> LColl_p, LColl_T, LColl_f;

	std::string line; double buffer;

	while (std::getline(file_in, line))
	{
        if (line.at(0) == '#')
            continue;
            
		std::stringstream ss(line);
		ss >> buffer; LColl_p.push_back(buffer);
		ss >> buffer; LColl_T.push_back(buffer);
		ss >> buffer; LColl_f.push_back(buffer);
	}

	file_in.close();

	m_LColl.setData(LColl_p, LColl_T, LColl_f);

	std::vector<std::vector<double>> domain = m_LColl.domain();
	if (m_Grids.pCollPts(0)  < domain[0][0]) {std::cerr << "Error: p grid point(s) out of lower bound of LColl domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.pCollPts(-1) > domain[0][1]) {std::cerr << "Error: p grid point(s) out of upper bound of LColl domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.TCollPts(0)  < domain[1][0]) {std::cerr << "Error: T grid point(s) out of lower bound of LColl domain. Aborting..." << std::endl; return -1;}
	if (m_Grids.TCollPts(-1) > domain[1][1]) {std::cerr << "Error: T grid point(s) out of upper bound of LColl domain. Aborting..." << std::endl; return -1;}

	return 1;
}

int energyLoss::generateTempGrid()
{
	const std::string path_in = "./evols/evols_cent=" + m_centrality + "/evolgridparams.dat";

	std::ifstream file_in(path_in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open evolution grid parameters file." << std::endl;
		return -1;
	}

	std::string line; double buffer;

	{//tau
		std::getline(file_in, line);
		std::getline(file_in, line);
		std::stringstream ss(line); ss >> buffer; m_tempTau0    = buffer;
							        ss >> buffer; m_tempTauStep = buffer;
							   			          m_tempTauMax  = 0;
		double tau = m_tempTau0; while (tau < (m_tauMaxFM+m_tempTauStep)) {m_tempTauMax++; tau+=m_tempTauStep;}
		m_tau0 = m_tempTau0;
	}

	{//x
		double xMax;
		std::getline(file_in, line);
		std::getline(file_in, line);
		std::stringstream ss(line); ss >> buffer;    m_tempX0 = buffer;
							        ss >> buffer;        xMax = buffer;
							        ss >> buffer; m_tempXStep = buffer;
		double x = m_tempX0; m_tempXMax = 0; while (x <= xMax) {m_tempXMax++; x+=m_tempXStep;}
	}

	{//y
		double yMax;
		std::getline(file_in, line);
		std::getline(file_in, line);
		std::stringstream ss(line); ss >> buffer;    m_tempY0 = buffer;
							        ss >> buffer;        yMax = buffer;
							        ss >> buffer; m_tempYStep = buffer;
		double y = m_tempY0; m_tempYMax = 0; while (y <= yMax) {m_tempYMax++; y+=m_tempYStep;}
	}

	for (size_t iTau=0; iTau<m_tempTauMax; iTau++) {
		for (size_t iX=0; iX<m_tempXMax; iX++) {
			for (size_t iY=0; iY<m_tempYMax; iY++) {
				m_tempTauGrid.push_back(m_tempTau0 + iTau*m_tempTauStep);
				  m_tempXGrid.push_back(  m_tempX0 +   iX*m_tempXStep);
				  m_tempYGrid.push_back(  m_tempY0 +   iY*m_tempYStep);
			}
		}		
	}

	return 1;
}

int energyLoss::loadPhiPoints()
{
	const std::string path_in = "./phiGaussPts/phiptsgauss" + std::to_string(m_phiGridN) + ".dat";
	std::ifstream file_in(path_in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open phi points file. Aborting..." << std::endl;
		return -1;
	}

	std::string line; double buffer;

	while(std::getline(file_in, line))
	{
		std::stringstream ss(line);
		ss >> buffer; m_phiGridPts.push_back(buffer);
	}

	file_in.close();

	if (m_phiGridN != m_phiGridPts.size()) {
		std::cerr << "Error: phiGridN not equal to number of point imported from a file. Aborting..." << std::endl;
		return -2;
	}

	return 1;
}

int energyLoss::loadBinCollPoints(size_t event_id, std::vector<std::vector<double>> &bcpoints)
{
	const std::string path_in = "./binarycollpts/binarycollpts_cent=" + m_centrality + "/binarycollpts" + std::to_string(event_id) + ".dat";

	std::ifstream file_in(path_in, std::ios_base::in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open binary collision points file for event: " + std::to_string(event_id) + "." << std::endl;
		return -1;
	}

	std::string line; double buffer;

	std::vector<double> xpoints, ypoints;

	while (std::getline(file_in, line))
	{
		if (line.length() == 0)
            continue;
        
        if (line.at(0) == '#')
            continue;

		std::stringstream ss(line);
		ss >> buffer; xpoints.push_back(buffer);
		ss >> buffer; ypoints.push_back(buffer);
	}

	bcpoints.resize(xpoints.size());

	for (size_t iBCP=0; iBCP<xpoints.size(); iBCP++)
	{
		bcpoints[iBCP].push_back(xpoints[iBCP]);
        bcpoints[iBCP].push_back(ypoints[iBCP]);
	}

	file_in.close();

	return 1;
}

int energyLoss::generateInitPosPoints(size_t event_id, std::vector<double> &xPoints, std::vector<double> &yPoints)
{
	std::vector<std::vector<double>> bcpts; if (loadBinCollPoints(event_id, bcpts) != 1) return -1;

	size_t bsptsNum = static_cast<size_t>(m_BCPP*static_cast<double>(bcpts.size()));

	if (bsptsNum < 1) bsptsNum = 1;

	if (m_BCPSEED == 0) {
		std::random_device rd; auto rng = std::default_random_engine{rd()};
		std::shuffle(bcpts.begin(), bcpts.end(), rng);
	}
	else {
		auto rng = std::default_random_engine{static_cast<long unsigned int>(m_BCPSEED)};
		std::shuffle(bcpts.begin(), bcpts.end(), rng);
	}

	for (size_t iBCP=0; iBCP<bsptsNum; iBCP++) {
		xPoints.push_back(bcpts[iBCP][0]);
        yPoints.push_back(bcpts[iBCP][1]);
	}

	return 1;
}

int energyLoss::loadTProfile(size_t event_id, interpolationF<double> &tempProfile)
{
	const std::string path_in = "./evols/evols_cent=" + m_centrality + "/tempevol" + std::to_string(event_id) + ".dat";

	std::ifstream file_in(path_in, std::ios_base::in | std::ios_base::binary);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open temperature evolution file for event " + std::to_string(event_id) + "." << std::endl;
		return -1;
	}

    std::vector<double> temps; float buffer;

    while (true) {
        file_in.read((char*)&buffer, sizeof(buffer));
        if (file_in.eof()) break;
        temps.push_back(static_cast<double>(buffer));
    }

	file_in.close();

	if (temps.size() > m_tempTauGrid.size()) {
		std::cerr << "Error: imported profile's size larger than grid size." << std::endl;
		return -2;
	}

	tempProfile.setData(m_tempTauGrid, m_tempXGrid, m_tempYGrid, temps);

	return 1;
}

void energyLoss::generateGaussTab(std::vector<double> &qGTab, std::vector<double> &fGTab) const
//function that generates sampling points for Gaussian integration
//qGTab, fGTab - vectors that store sampling point <- output
{	
	double sigmaNum = 3.5; //setting sigma
	double sigmaStep = 0.25; //setting step
	size_t GTabLen = 2 * static_cast<size_t>(sigmaNum / sigmaStep) + 1; //setting length of sampling points
	
	double GaussTabSum = 0.0; //setting normalization sum to zero
	
	for (size_t iG=0; iG<GTabLen; iG++) //calculating sampling points
	{
		qGTab.push_back(-1.0*sigmaNum + static_cast<double>(iG)*sigmaStep); //setting qGaussTab values
		fGTab.push_back(std::exp(-qGTab.back()*qGTab.back()/2.0));          //setting fGaussTab values
		GaussTabSum += fGTab.back();                                        //adding to normalization sum
	}
	
	for (size_t iG=0; iG<GTabLen; iG++)  //normalizing
	{
		fGTab[iG] /= GaussTabSum; //dividing fGaussTab values with total sum
	}
}

void energyLoss::calculateAvgPathlenTemps(const std::vector<double> &pathLenghDist, const std::vector<double> &temperatureDist, std::vector<double> &avgPathLength, std::vector<double> &avgTemp) const
{
	interpolationF<double> pathLenghDistInt(m_phiGridPts, pathLenghDist);
	avgPathLength.push_back(poly::cubicIntegrate(m_phiGridPts, pathLenghDist)/2.0/M_PI);
	avgPathLength.push_back((pathLenghDistInt.interpolation(m_phiGridPts.front()) + pathLenghDistInt.interpolation(m_phiGridPts.back()))/2.0);
	avgPathLength.push_back((pathLenghDistInt.interpolation(M_PI/2.0)             + pathLenghDistInt.interpolation(3.0*M_PI/2.0))       /2.0);

	interpolationF<double> temperatureDistInt(m_phiGridPts, temperatureDist);
	avgTemp.push_back(poly::cubicIntegrate(m_phiGridPts, temperatureDist)/2.0/M_PI);
	avgTemp.push_back((temperatureDistInt.interpolation(m_phiGridPts.front()) + temperatureDistInt.interpolation(m_phiGridPts.back()))/2.0);
	avgTemp.push_back((temperatureDistInt.interpolation(M_PI/2.0)             + temperatureDistInt.interpolation(3.0*M_PI/2.0))       /2.0);
}

int energyLoss::exportResults(const std::string &particleName, size_t event_id, const std::vector<std::vector<double>> &RAApTphi, const std::vector<double> &avgPathLength, const std::vector<double> &avgTemp, size_t trajecNum, size_t elossNum) const
{
	std::vector<std::string> header;
    header.push_back("#collision_system: " + m_collsys);
	header.push_back("#collision_energy: " + m_sNN);
	header.push_back("#particle_type: " + particleName);
	header.push_back("#centrality: " + m_centrality);

	std::stringstream xbsstr; xbsstr << std::fixed << std::setprecision(1) << m_xB;
	header.push_back("#xB = " + xbsstr.str());

	header.push_back("#event_id: " + std::to_string(event_id));

	std::stringstream avgPathLengthSStr[3];
    for (size_t i=0; i<3; i++) avgPathLengthSStr[i] << std::fixed << std::setprecision(6) << avgPathLength[i];
	header.push_back("#average_path-lengths: " + avgPathLengthSStr[0].str() + ", " + avgPathLengthSStr[1].str() + ", " + avgPathLengthSStr[2].str());

	std::stringstream avgTempSStr[3];
    for (size_t i=0; i<3; i++) avgTempSStr[i] << std::fixed << std::setprecision(6) << avgTemp[i];
	header.push_back("#average_temperatures: " + avgTempSStr[0].str() + ", " + avgTempSStr[1].str() + ", " + avgTempSStr[2].str());
	
	header.push_back("#number_of_angles:                " + std::to_string(m_phiGridN));
	
	header.push_back("#total_number_of_trajectories:    " + std::to_string(trajecNum));
	header.push_back("#total_number_of_jet_energy_loss: " + std::to_string(elossNum));

	header.push_back("#BCPSEED: " + std::to_string(m_BCPSEED));

	header.push_back("#-------------------------------------------------------");
	header.push_back("#   pT [GeV]       phi          R_AA   ");

	//setting file path:
	const std::string path_out = "./results/results" + particleName + "/" + particleName + "_" + m_collsys + "_sNN=" + m_sNN + "_cent=" + m_centrality + "_xB=" + xbsstr.str() + "_dist_" + std::to_string(event_id) + ".dat";

	std::ofstream file_out(path_out, std::ios_base::out);
	if (!file_out.is_open()) {
		std::cerr << "Error: unable to open RAA(pT,phi) distribution file for event " + std::to_string(event_id) + "." << std::endl;
		return -1;
	}

	for (const auto &h : header) file_out << h << "\n";

	for (size_t ipT= 0; ipT<m_Grids.finPtsLength(); ipT++) //printing RAA(pT,phi) to file
		for (size_t iPhi=0; iPhi<m_phiGridN; iPhi++) {
			file_out << std::fixed << std::setw(14) << std::setprecision(10) <<   m_Grids.finPts(ipT) << " ";
			file_out << std::fixed << std::setw(12) << std::setprecision(10) <<    m_phiGridPts[iPhi] << " ";
			file_out << std::fixed << std::setw(12) << std::setprecision(10) << RAApTphi[ipT][iPhi] << "\n";
		}

	file_out.close();

	return 1;
}


void energyLoss::RadCollEL(double X0, double Y0, double phi0, const interpolationF<double> &TProfile, std::vector<double> &radiativeRAA1, std::vector<std::vector<double>> &radiativeRAA2, std::vector<double> &collisionalEL, double &pathLength, double &temp) const
//function that calculates radiative and collisional EL for particles created in (X0, Y0) with direction phi0 (modefied pT integration algorithm)
//X0, Y0, phi0  - inital position and angle 					  		     <- input
//radiativeRAA1 - radiative RAA for single trajectory (dA410)	  		     <- output
//radiativeRAA2 - radiative RAA for single trajectory (rest of dA integrals) <- output
//collisionalEL - collisional energy loss for single trajectory   		     <- output
//pathL, temp - path-length and temperature for single trajectory 		     <- output
{
	std::vector<double> currLTTabL, currLTTabT; //defining arrays that will store current path-lengths and temperatures

	double t = m_tau0, currTemp; //defining current path-length (time) and temperature

	while ((currTemp = TProfile.interpolation(t, X0 + t*std::cos(phi0), Y0 + t*std::sin(phi0))) > m_TCRIT) { //calculating current path-length and temp table
		currLTTabL.push_back(t);
		currLTTabT.push_back(currTemp);
		t += m_TIMESTEP;
	}
	
	if (currLTTabL.size() > 1) { //calculating energy loss if path-length is longer than thermalization time
		
		///////////////////////////////////////////////////////////////////////////////////////////////////
		//Radiative EnergyLoss calculation:

		std::vector<double> currNormTabTau(currLTTabL.size()), currNormTabVal(currLTTabL.size()); //LNorm table to be integrated over tau
		std::vector<double> NormSparseP, NormSparseV;											  //table for currNormInterp
		
		std::vector<double> currDndxTabTau(currLTTabL.size()), currDndxTabVal(currLTTabL.size()); //Ldndx table to be integrated over tau
		std::vector<double> dndxSparseP, dndxSparseX, dndxSparseV;			  				 	  //table for currDndxInterp

		for (const auto &p : m_Grids.pPts()) //loop over ppts
		{
			for (size_t iL=0; iL<currLTTabL.size(); iL++) //loop over current path-length and temperature table
			{
				currNormTabTau[iL] = currLTTabL[iL]; 								         //setting path-lengths
				currNormTabVal[iL] = m_LNorm.interpolation(currLTTabL[iL], p, currLTTabT[iL]); //setting current norm values by integrating over time
			}

			NormSparseP.push_back(p);												//setting p of current norm table
			NormSparseV.push_back(poly::linearIntegrate(currNormTabTau, currNormTabVal)); //setting value of current norm table

			for (const auto &x : m_Grids.xPts()) //loop over xpts
			{
				for (size_t iL=0; iL<currLTTabL.size(); iL++) //loop over current path-length and temperature table
				{
					currDndxTabTau[iL] = currLTTabL[iL]; 									        //setting path-lengths
					currDndxTabVal[iL] = m_Ldndx.interpolation(currLTTabL[iL], p, currLTTabT[iL], x); //setting Ldndx values
				}

				dndxSparseP.push_back(p); 												//setting p of current dndx table
				dndxSparseX.push_back(x);												//setting x of current dndx table
				dndxSparseV.push_back(poly::linearIntegrate(currDndxTabTau, currDndxTabVal)); //setting curernt dndx values by integrating over time
			}
		}
		
		interpolationF<double> currNorm(NormSparseP, NormSparseV); 			    //constructing interpolated current norm
		interpolationF<double> currDndx(dndxSparseP, dndxSparseX, dndxSparseV); //constructing interpolated current dndx
		
		for (const auto &ph : m_Grids.RadPts()) //loop over Radpts
		{
			radiativeRAA1.push_back(dAp410(ph, currNorm)); //calculating radiative energy loss for dA410

			radiativeRAA2.push_back(std::vector<double>()); //resizing radiativeRAA2 2d vector

			for (const auto &Fdp : m_Grids.FdpPts())
				radiativeRAA2.back().push_back(FdA(ph, Fdp, currNorm, currDndx)); //calculating radiative energy loss for rest of the dA integrals

		}
		
		///////////////////////////////////////////////////////////////////////////////////////////////////
		//Collisional EnergyLoss calculation:

		std::vector<double> currCollTabTau(currLTTabL.size()), currCollTabVal(currLTTabL.size()); //collisional table to be integrated over tau

		for (const auto& p : m_Grids.pCollPts()) //loop over pCollPts
		{
			for (size_t iL=0; iL<currLTTabL.size(); iL++) //loop over current path-length and temperature table
			{
				currCollTabTau[iL] = currLTTabL[iL]; 				         //setting path-lengths
				currCollTabVal[iL] = m_LColl.interpolation(p, currLTTabT[iL]); //setting LColl values
			}

			collisionalEL.push_back(poly::linearIntegrate(currCollTabTau, currCollTabVal)); //calculating collisional energy loss by integrating over time
		}
		
		pathLength = currLTTabL.back(); //setting value of path-length for single trajectory

		//calculating mean temperature along path
		temp = 0.0;
		for (size_t iL=0; iL<currLTTabL.size(); iL++) temp += currLTTabT[iL];
		temp /= static_cast<double>(currLTTabL.size());
	}
	else { //if path-length is smaller than thermalization time:

		pathLength = 0.0; //setting path-length and temperature
		temp       = 0.0;
	}
}

void energyLoss::RadCollEL(double X0, double Y0, double phi0, const interpolationF<double> &TProfile, std::vector<double> &radiativeRAA, std::vector<double> &collisionalEL, double &pathLenght, double &temp) const
//function that calculates radiative and collisional EL for particles created in (X0, Y0) with direction phi0 (standard algorithm)
//X0, Y0, phi0  - inital position and angle 					  <- input
//radiativeRAA  - radiative RAA for single trajectory 			  <- output
//collisionalEL - collisional energy loss for single trajectory   <- output
//pathL, temp - path-length and temperature for single trajectory <- output
{
	std::vector<double> currLTTabL, currLTTabT; //defining arrays that will store current path-lengths and temperatures

	double t = m_tau0, currTemp; //defining current path-length (time) and temperature

	while ((currTemp = TProfile.interpolation(t, X0 + t*cos(phi0), Y0 + t*sin(phi0))) > m_TCRIT) { //calculating current path-length and temp table
		currLTTabL.push_back(t);
		currLTTabT.push_back(currTemp);
		t += m_TIMESTEP;
	}
	
	if (currLTTabL.size() > 1) { //calculating energy loss if path-length is longer than thermalization time
		
		///////////////////////////////////////////////////////////////////////////////////////////////////
		//Radiative EnergyLoss calculation:

		std::vector<double> currNormTabTau(currLTTabL.size()), currNormTabVal(currLTTabL.size()); //LNorm table to be integrated over tau
		std::vector<double> NormSparseP, NormSparseV;											 //table for currNormInterp
		
		std::vector<double> currDndxTabTau(currLTTabL.size()), currDndxTabVal(currLTTabL.size()); //Ldndx table to be integrated over tau
		std::vector<double> dndxSparseP, dndxSparseX, dndxSparseV;			  				 	 //table for currDndxInterp

		for (const auto &p : m_Grids.pPts()) //loop over ppts
		{
			for (size_t iL=0; iL<currLTTabL.size(); iL++) //loop over current path-length and temperature table
			{
				currNormTabTau[iL] = currLTTabL[iL]; 								         //setting path-lengths
				currNormTabVal[iL] = m_LNorm.interpolation(currLTTabL[iL], p, currLTTabT[iL]); //setting current norm values by integrating over time
			}

			NormSparseP.push_back(p);												//setting p of current norm table
			NormSparseV.push_back(poly::linearIntegrate(currNormTabTau, currNormTabVal)); //setting value of current norm table

			for (const auto &x : m_Grids.xPts()) //loop over xpts
			{
				for (size_t iL=0; iL<currLTTabL.size(); iL++) //loop over current path-length and temperature table
				{
					currDndxTabTau[iL] = currLTTabL[iL]; 									        //setting path-lengths
					currDndxTabVal[iL] = m_Ldndx.interpolation(currLTTabL[iL], p, currLTTabT[iL], x); //setting Ldndx values
				}

				dndxSparseP.push_back(p); 												//setting p of current dndx table
				dndxSparseX.push_back(x);												//setting x of current dndx table
				dndxSparseV.push_back(poly::linearIntegrate(currDndxTabTau, currDndxTabVal)); //setting curernt dndx values by integrating over time
			}
		}
		
		interpolationF<double> currNorm(NormSparseP, NormSparseV); 			   //constructing interpolated current norm
		interpolationF<double> currDndx(dndxSparseP, dndxSparseX, dndxSparseV); //constructing interpolated current dndx
		
		for (const auto &p : m_Grids.RadPts())
			radiativeRAA.push_back(dA41(p, currNorm, currDndx)/m_dsdpti2.interpolation(p)); //calculating radiative RAA
		
		///////////////////////////////////////////////////////////////////////////////////////////////////
		//Collisional EnergyLoss calculation:

		std::vector<double> currCollTabTau(currLTTabL.size()), currCollTabVal(currLTTabL.size()); //collisional table to be integrated over tau

		for (const auto &p : m_Grids.pCollPts()) //loop over pCollPts
		{
			for (size_t iL=0; iL<currLTTabL.size(); iL++) //loop over current path-length and temperature table
			{
				currCollTabTau[iL] = currLTTabL[iL]; 				         //setting path-lengths
				currCollTabVal[iL] = m_LColl.interpolation(p, currLTTabT[iL]); //setting LColl values
			}

			collisionalEL.push_back(poly::linearIntegrate(currCollTabTau, currCollTabVal)); //calculating collisional energy loss by integrating over time
		}
		
		pathLenght = currLTTabL.back(); //setting value of path-length for single trajectory

		//calculating mean temperature along path
		temp = 0.0;
		for (size_t iL=0; iL<currLTTabL.size(); iL++) temp += currLTTabT[iL];
		temp /= static_cast<double>(currLTTabL.size());
	}
	else { //if path-length is smaller than thermalization time:

		pathLenght = 0.0; //setting path-length and temperature
		     temp  = 0.0;
	}
}


void energyLoss::runELossHeavyFlavour()
{
	if (loaddsdpti2(m_pName, m_dsdpti2) != 1) return;

	FdAHaltonSeqInit(150);

	#pragma omp parallel for schedule(dynamic)
	for (size_t eventID=0; eventID<m_eventN; eventID++)
	{
		std::vector<double> xPoints, yPoints; generateInitPosPoints(eventID, xPoints, yPoints);

		interpolationF<double> tProfile; loadTProfile(eventID, tProfile);

		std::vector<std::vector<double>> RAAdist(m_Grids.finPtsLength(), std::vector<double>(m_phiGridN, 0.0));

		std::vector<double> pathLenghDist(m_phiGridN, 0.0), temperatureDist(m_phiGridN, 0.0);

		size_t trajectoryNum = 0, energylossNum = 0;

		for (size_t iPhi=0; iPhi<m_phiGridN; iPhi++)
		{
			double phi = m_phiGridPts[iPhi];

			std::vector<double> sumRAA1(m_Grids.finPtsLength(), 0.0);

			std::vector<std::vector<double>> sumRAA2(m_Grids.finPtsLength(), std::vector<double>(m_Grids.FdpPtsLength(), 0.0));

			size_t pltCNT = 0; //path-length and temperature distribution counter

			for (size_t iXY=0; iXY<xPoints.size(); iXY++)
			{
				trajectoryNum++;

				double x = xPoints[iXY], y = yPoints[iXY];

				std::vector<double> radRAA1; std::vector<std::vector<double>> radRAA2; std::vector<double> collEL;
				double pathLength, temperature;
				RadCollEL(x, y, phi, tProfile, radRAA1, radRAA2, collEL, pathLength, temperature);

				if (pathLength > m_tau0) { //checking if path-length is larger than thermalization time

					energylossNum++;

					pltCNT++;
					pathLenghDist[iPhi] += pathLength;
					temperatureDist[iPhi] += temperature;

					for (auto &coll : collEL) coll += 1e-12; //modifying collEL to prevent division by 0

					std::vector<double> singleRAA1; std::vector<std::vector<double>> singleRAA2;
					gaussFilterIntegrate(radRAA1, radRAA2, collEL, singleRAA1, singleRAA2);

					for (size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++) {
						sumRAA1[iFinPts] += singleRAA1[iFinPts];
						for (size_t iFdp=0; iFdp<m_Grids.FdpPtsLength(); iFdp++)
							sumRAA2[iFinPts][iFdp] += singleRAA2[iFinPts][iFdp];
					}
				}
				else { //if path length is smaller than tau0:

					for (size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++) //adding RAA1, which is 1.0, to RAA sum; RAA2 is 0 in this case
						sumRAA1[iFinPts] += 1.0; 
				}
			}

			double weightsum = static_cast<double>(xPoints.size());
			std::for_each(sumRAA1.begin(), sumRAA1.end(), [weightsum](double &c){ c/=weightsum; });
			for (size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++)
				std::for_each(sumRAA2[iFinPts].begin(), sumRAA2[iFinPts].end(), [weightsum](double &c){ c/=weightsum; });

			//setting RAA(pT,phi) value by integrating over p:
			for (size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++)
				RAAdist[iFinPts][iPhi] = sumRAA1[iFinPts] + poly::cubicIntegrate(m_Grids.FdpPts(), sumRAA2[iFinPts])/m_Grids.finPts(iFinPts);

			pathLenghDist[iPhi] /= static_cast<double>(pltCNT); temperatureDist[iPhi] /= static_cast<double>(pltCNT);
		}

		std::vector<double> avgPathLength, avgTemp;
		calculateAvgPathlenTemps(pathLenghDist, temperatureDist, avgPathLength, avgTemp);
		
		exportResults(m_pName, eventID, RAAdist, avgPathLength, avgTemp, trajectoryNum, energylossNum);
	}
}

void energyLoss::gaussFilterIntegrate(const std::vector<double> &radiativeRAA1, const std::vector<std::vector<double>> &radiativeRAA2, const std::vector<double> &collisionalEL, std::vector<double> &singRAA1, std::vector<std::vector<double>> &singRAA2) const
//function that performs Gauss filter integration - modefied pT integration algorithm
//radiativeRAA1 - raditive RAA (dA410)											  <- input
//radiativeRAA2 - raditive RAA (rest of dA integrals)							  <- input
//collisionalEL - collisional energy loss										  <- input
//singRAA1 		- RAA array after Gauss filter integration (dA410)				  <- output
//singRAA2 		- RAA array after Gauss filter integration (rest of dA integrals) <- output
{
    interpolationF<double> muCollInt(m_Grids.pCollPts(), collisionalEL); //creating collisional energy loss interpolated function

	std::vector<double> qGaussTabOG, fGaussTabOG; //defining vectors that will store original Gauss filter sampling points
	generateGaussTab(qGaussTabOG, fGaussTabOG);   //generating sampling points and settin number of sampling poins

	std::vector<double> qGaussTab, fGaussTab; //defining vectors that will store Gauss filter sampling points

	//////////////////////////////////////////////////////////////////////////////////
	//Gauss integration of dAp410:
	{
        interpolationF<double> RadRelInt(m_Grids.RadPts(), radiativeRAA1); //creating radiative RAA1 interpolated function

		double GFSum; //defining sum variable for Gauss filter
		double dppT;  //defining integration variable

		double muCollCurrVal; //defining variable that stores value of interpolated muColl for specific pT, ie current value
		double sigmaColl;     //defining variable for collisional sigma

		for (const auto &pT : m_Grids.finPts())
		{
			GFSum = 0.0;

			muCollCurrVal = muCollInt.interpolation(pT);

			sigmaColl = std::sqrt(2.0*m_TCollConst*muCollCurrVal);

			qGaussTab = qGaussTabOG; fGaussTab = fGaussTabOG; //setting Gauss filter

			if ((muCollCurrVal + sigmaColl * qGaussTab.front()) < -3.0) { 						        //checking if Gauss is out of bound on lower bound
				double resfac = ((-3.0 + 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.front(); 	        //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}			
		
			if ((muCollCurrVal + sigmaColl * qGaussTab.back()) > 20.0) {						        //checking if Gauss is out of bound on upper bound
				double resfac = ((20.0 - 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.back();	        //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}

			//calculating Gauss filter
			for (size_t iG=0; iG<qGaussTab.size(); iG++)
			{
				dppT = muCollCurrVal + sigmaColl * qGaussTab[iG];			
				GFSum += (m_dsdpti2.interpolation(pT + dppT)*RadRelInt.interpolation(pT + dppT)*(pT + dppT) / pT * fGaussTab[iG]);
			}

			singRAA1.push_back(1.0 / m_dsdpti2.interpolation(pT) * GFSum);
		}
	}

	//////////////////////////////////////////////////////////////////////////////////
	//Gauss integration of FdA:
	{
		interpolationF<double> RadRelInt(m_Grids.RadPts(), m_Grids.FdpPts(), radiativeRAA2);

		double GFSum; //defining sum variable for Gauss filter
		double dppT;  //defining integration variable

		double muCollCurrVal; //defining variable that stores value of interpolated muColl for specific pT, ie current value
		double sigmaColl;     //defining variable for collisional sigma

		for (const auto &pT : m_Grids.finPts())
		{
			singRAA2.push_back(std::vector<double>()); //resizing single RAA vector

			muCollCurrVal = muCollInt.interpolation(pT);

			sigmaColl = std::sqrt(2.0*m_TCollConst*muCollCurrVal);

			qGaussTab = qGaussTabOG; fGaussTab = fGaussTabOG; //setting Gauss filter

			if ((muCollCurrVal + sigmaColl * qGaussTab.front()) < -3.0) { 						        //checking if Gauss is out of bound on lower bound
				double resfac = ((-3.0 + 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.front(); 	        //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}			
		
			if ((muCollCurrVal + sigmaColl * qGaussTab.back()) > 20.0) {						        //checking if Gauss is out of bound on upper bound
				double resfac = ((20.0 - 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.back();            //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}

			for (const auto &dpT : m_Grids.FdpPts()) //loop over FdpPts
			{
				GFSum = 0.0; //setting sum to 0

				//calculating Gauss filter
				for (size_t iG=0; iG<qGaussTab.size(); iG++)
				{
					dppT = muCollCurrVal + sigmaColl * qGaussTab[iG];
					GFSum += (m_dsdpti2.interpolation(pT + dpT + dppT)*RadRelInt.interpolation(pT + dppT, dpT)*(pT + dppT)/(pT+ dpT + dppT)*fGaussTab[iG]);
				}

				singRAA2.back().push_back(1.0 / m_dsdpti2.interpolation(pT) * GFSum);
			}
		}
	}
}


void energyLoss::runELossLightQuarks()
{
	const std::vector<std::string> lightQuarksList{"Down", "DownBar", "Strange", "Up", "UpBar"};

	std::vector<interpolationF<double>> dsdpti2LightQuarks(lightQuarksList.size());

	for (size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++)
		if (loaddsdpti2(lightQuarksList[iLQ], dsdpti2LightQuarks[iLQ]) != 1) return;

	FdAHaltonSeqInit(100);

	#pragma omp parallel for schedule(dynamic)
	for (size_t eventID=0; eventID<m_eventN; eventID++)
	{
		std::vector<double> xPoints, yPoints; generateInitPosPoints(eventID, xPoints, yPoints);

		interpolationF<double> tProfile; loadTProfile(eventID, tProfile);

		std::vector<std::vector<std::vector<double>>> RAAdist(lightQuarksList.size(), std::vector<std::vector<double>>(m_Grids.finPtsLength(), std::vector<double>(m_phiGridN, 0.0)));

		std::vector<double> pathLenghDist(m_phiGridN, 0.0), temperatureDist(m_phiGridN, 0.0);

		size_t trajectoryNum = 0, energylossNum = 0;

		for (size_t iPhi=0; iPhi<m_phiGridN; iPhi++)
		{
			double phi = m_phiGridPts[iPhi];

			std::vector<std::vector<double>> sumRAA1(lightQuarksList.size(), std::vector<double>(m_Grids.finPtsLength(), 0.0));

			std::vector<std::vector<std::vector<double>>> sumRAA2(lightQuarksList.size(), std::vector<std::vector<double>>(m_Grids.finPtsLength(), std::vector<double>(m_Grids.FdpPtsLength(), 0.0)));

			size_t pltCNT = 0; //path-length and temperature distribution counter

			for (size_t iXY=0; iXY<xPoints.size(); iXY++) //loop over x and y initial position points
			{
				trajectoryNum++;

				double x = xPoints[iXY], y = yPoints[iXY];

				std::vector<double> radRAA1; std::vector<std::vector<double>> radRAA2; std::vector<double> collEL;
				double pathLength, temperature;
				RadCollEL(x, y, phi, tProfile, radRAA1, radRAA2, collEL, pathLength, temperature);

				if (pathLength > m_tau0) { //checking if path-length is larger than thermalization time

					energylossNum++; //adding to number of energy loss calculations

					pltCNT++;
					pathLenghDist[iPhi] += pathLength;
					temperatureDist[iPhi] += temperature;

					for (auto &coll : collEL) coll += 1e-12; //modifying collEL to prevent division by 0

					std::vector<std::vector<double>> singleRAA1(lightQuarksList.size());
					std::vector<std::vector<std::vector<double>>> singleRAA2(lightQuarksList.size());

					for (size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++)
						gaussFilterIntegrate(dsdpti2LightQuarks[iLQ], radRAA1, radRAA2, collEL, singleRAA1[iLQ], singleRAA2[iLQ]);

					for (size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++) {
						for (size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++) {
							sumRAA1[iLQ][iFinPts] += singleRAA1[iLQ][iFinPts];
							for (size_t iFdp=0; iFdp<m_Grids.FdpPtsLength(); iFdp++)
								sumRAA2[iLQ][iFinPts][iFdp] += singleRAA2[iLQ][iFinPts][iFdp];
						}
					}
				}
				else {
					for (size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++)
						for (size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++)
							sumRAA1[iLQ][iFinPts] += 1.0;
				}
			}

			double weightsum = static_cast<double>(xPoints.size());
			for (size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++) {
				std::for_each(sumRAA1[iLQ].begin(), sumRAA1[iLQ].end(), [weightsum](double &c){ c/=weightsum; });
				for (size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++)
					std::for_each(sumRAA2[iLQ][iFinPts].begin(), sumRAA2[iLQ][iFinPts].end(), [weightsum](double &c){ c/=weightsum; });
			}

			//setting RAA(pT,phi) value by integrating over p:
			for (size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++)
				for (size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++)
					RAAdist[iLQ][iFinPts][iPhi] = sumRAA1[iLQ][iFinPts] + poly::cubicIntegrate(m_Grids.FdpPts(), sumRAA2[iLQ][iFinPts])/m_Grids.finPts(iFinPts);

			pathLenghDist[iPhi] /= static_cast<double>(pltCNT); temperatureDist[iPhi] /= static_cast<double>(pltCNT);
		}

		std::vector<double> avgPathLength, avgTemp;
		calculateAvgPathlenTemps(pathLenghDist, temperatureDist, avgPathLength, avgTemp);
		
		for (size_t iLQ=0; iLQ<lightQuarksList.size(); iLQ++)
			exportResults(lightQuarksList[iLQ], eventID, RAAdist[iLQ], avgPathLength, avgTemp, trajectoryNum, energylossNum);
	}
}

void energyLoss::gaussFilterIntegrate(const interpolationF<double> &dsdpti2lquark, const std::vector<double> &radiativeRAA1, const std::vector<std::vector<double>> &radiativeRAA2, const std::vector<double> &collisionalEL, std::vector<double> &singRAA1, std::vector<std::vector<double>> &singRAA2) const
//function that performs Gauss filter integration - modefied pT integration algorithm used in all lquarks algorithm
//dsdpti2lquark - light quark initial pT distribution      						  <- input
//radiativeRAA1 - raditive RAA (dA410)											  <- input
//radiativeRAA2 - raditive RAA (rest of dA integrals)							  <- input
//collisionalEL - collisional energy loss										  <- input
//singRAA1 		- RAA array after Gauss filter integration (dA410)				  <- output
//singRAA2 		- RAA array after Gauss filter integration (rest of dA integrals) <- output
{
    interpolationF<double> muCollInt(m_Grids.pCollPts(), collisionalEL); //creating collisional energy loss interpolated function

	std::vector<double> qGaussTabOG, fGaussTabOG; //defining vectors that will store original Gauss filter sampling points
	generateGaussTab(qGaussTabOG, fGaussTabOG);   //generating sampling points and settin number of sampling poins

	std::vector<double> qGaussTab, fGaussTab; //defining vectors that will store Gauss filter sampling points

	//////////////////////////////////////////////////////////////////////////////////
	//Gauss integration of dAp410:
	{
        interpolationF<double> RadRelInt(m_Grids.RadPts(), radiativeRAA1); //creating radiative RAA1 interpolated function

		double GFSum; //defining sum variable for Gauss filter
		double dppT;  //defining integration variable

		double muCollCurrVal; //defining variable that stores value of interpolated muColl for specific pT, ie current value
		double sigmaColl;     //defining variable for collisional sigma

		for (const auto &pT : m_Grids.finPts())
		{
			GFSum = 0.0;

			muCollCurrVal = muCollInt.interpolation(pT);

			sigmaColl = std::sqrt(2.0*m_TCollConst*muCollCurrVal);

			qGaussTab = qGaussTabOG; fGaussTab = fGaussTabOG; //setting Gauss filter

			if ((muCollCurrVal + sigmaColl * qGaussTab.front()) < -3.0) { 						        //checking if Gauss is out of bound on lower bound
				double resfac = ((-3.0 + 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.front(); 	        //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}			
		
			if ((muCollCurrVal + sigmaColl * qGaussTab.back()) > 20.0) {						        //checking if Gauss is out of bound on upper bound
				double resfac = ((20.0 - 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.back();	        //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}

			//calculating Gauss filter
			for (size_t iG=0; iG<qGaussTab.size(); iG++)
			{
				dppT = muCollCurrVal + sigmaColl * qGaussTab[iG];			
				GFSum += (dsdpti2lquark.interpolation(pT + dppT)*RadRelInt.interpolation(pT + dppT)*(pT + dppT) / pT * fGaussTab[iG]);
			}

			singRAA1.push_back(1.0 / dsdpti2lquark.interpolation(pT) * GFSum);
		}
	}

	//////////////////////////////////////////////////////////////////////////////////
	//Gauss integration of FdA:
	{
		interpolationF<double> RadRelInt(m_Grids.RadPts(), m_Grids.FdpPts(), radiativeRAA2);

		double GFSum; //defining sum variable for Gauss filter
		double dppT;  //defining integration variable

		double muCollCurrVal; //defining variable that stores value of interpolated muColl for specific pT, ie current value
		double sigmaColl;     //defining variable for collisional sigma

		for (const auto &pT : m_Grids.finPts())
		{
			singRAA2.push_back(std::vector<double>()); //resizing single RAA vector

			muCollCurrVal = muCollInt.interpolation(pT);

			sigmaColl = std::sqrt(2.0*m_TCollConst*muCollCurrVal);

			qGaussTab = qGaussTabOG; fGaussTab = fGaussTabOG; //setting Gauss filter

			if ((muCollCurrVal + sigmaColl * qGaussTab.front()) < -3.0) { 						        //checking if Gauss is out of bound on lower bound
				double resfac = ((-3.0 + 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.front(); 	        //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}			
		
			if ((muCollCurrVal + sigmaColl * qGaussTab.back()) > 20.0) {						        //checking if Gauss is out of bound on upper bound
				double resfac = ((20.0 - 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.back();	        //setting rescaling factor
				std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
			}

			for (const auto &dpT : m_Grids.FdpPts()) //loop over FdpPts
			{
				GFSum = 0.0; //setting sum to 0

				//calculating Gauss filter
				for (size_t iG=0; iG<qGaussTab.size(); iG++)
				{
					dppT = muCollCurrVal + sigmaColl * qGaussTab[iG];
					GFSum += (dsdpti2lquark.interpolation(pT + dpT + dppT)*RadRelInt.interpolation(pT + dppT, dpT)*(pT + dppT)/(pT+ dpT + dppT)*fGaussTab[iG]);
				}

				singRAA2.back().push_back(1.0 / dsdpti2lquark.interpolation(pT) * GFSum);
			}
		}
	}
}


void energyLoss::runELossLightFlavour()
{
	if (loaddsdpti2(m_pName, m_dsdpti2)  != 1) return;

	dAHaltonSeqInit(1000);

	#pragma omp parallel for schedule(dynamic)
	for (size_t eventID=0; eventID<m_eventN; eventID++)
	{
		std::vector<double> xPoints, yPoints; generateInitPosPoints(eventID, xPoints, yPoints);

		interpolationF<double> tProfile; loadTProfile(eventID, tProfile);

		std::vector<std::vector<double>> RAAdist(m_Grids.finPtsLength(), std::vector<double>(m_phiGridN, 0.0));

		std::vector<double> pathLenghDist(m_phiGridN, 0.0), temperatureDist(m_phiGridN, 0.0);

		size_t trajectoryNum = 0, energylossNum = 0;

		for (size_t iPhi=0; iPhi<m_phiGridN; iPhi++)
		{
			double phi = m_phiGridPts[iPhi];

			std::vector<double> sumRAA(m_Grids.finPtsLength(), 0.0);

			size_t pltCNT = 0; //path-length and temperature distribution counter

			for (size_t iXY=0; iXY<xPoints.size(); iXY++)
			{
				trajectoryNum++;

				double x = xPoints[iXY], y = yPoints[iXY];

				std::vector<double> radRAA, collEL; double pathLength, temperature;
				RadCollEL(x, y, phi, tProfile, radRAA, collEL, pathLength, temperature);

				if (pathLength > m_tau0) { //checking if path-length is larger than thermalization time

					energylossNum++; //adding to number of energy loss calculations

					pltCNT++;
					pathLenghDist[iPhi] += pathLength;
					temperatureDist[iPhi] += temperature;

					for (auto &coll : collEL) coll += 1e-12; //modifying collEL to prevent division by 0

					std::vector<double> singleRAA;
					gaussFilterIntegrate(radRAA, collEL, singleRAA);

					for (size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++)
						sumRAA[iFinPts] += singleRAA[iFinPts];
				}
				else { //if path length is smaller than tau0:

					for (size_t iFinPts=0; iFinPts<m_Grids.finPtsLength(); iFinPts++)
						sumRAA[iFinPts] += 1.0;
				}
			}

			double weightsum = (double)(xPoints.size());
			for (size_t iFinPts= 0; iFinPts<m_Grids.finPtsLength(); iFinPts++)
				RAAdist[iFinPts][iPhi] = sumRAA[iFinPts]/weightsum;

			pathLenghDist[iPhi] /= static_cast<double>(pltCNT); temperatureDist[iPhi] /= static_cast<double>(pltCNT);
		}

		std::vector<double> avgPathLength, avgTemp;
		calculateAvgPathlenTemps(pathLenghDist, temperatureDist, avgPathLength, avgTemp);

		exportResults(m_pName, eventID, RAAdist, avgPathLength, avgTemp, trajectoryNum, energylossNum);
	}
}

void energyLoss::gaussFilterIntegrate(const std::vector<double> &radiativeRAA, const std::vector<double> &collisionalEL, std::vector<double> &singRAA) const
//function that performs Gauss filter integration - default algorithm
//radiativeRAA  - raditive RAA 							   <- input
//collisionalEL - collisional energy loss				   <- input
//singRAA 		- RAA array after Gauss filter integration <- output
{
    interpolationF<double> RadRelInt(m_Grids.RadPts(),   radiativeRAA);  //creating radiative RAA interpolated function
    interpolationF<double> muCollInt(m_Grids.pCollPts(), collisionalEL); //creating collisional energy loss interpolated function

	std::vector<double> qGaussTabOG, fGaussTabOG; //defining vectors that will store original Gauss filter sampling points
	generateGaussTab(qGaussTabOG, fGaussTabOG);   //generating sampling points and settin number of sampling poins

	std::vector<double> qGaussTab, fGaussTab; //defining vectors that will store Gauss filter sampling points

	double GFSum; //defining sum variable for Gauss filter

	double dpT; //defining pT and dpT variables

	double muCollCurrVal; //defining variable that stores value of interpolated muColl for specific pT, ie current value

	double sigmaColl; //defining variable for collisional sigma
	
	//Gauss filter
	for (const auto &pT : m_Grids.finPts())
	{
		GFSum = 0.0L;

		muCollCurrVal = muCollInt.interpolation(pT);

		sigmaColl = std::sqrt(2.0*m_TCollConst*muCollCurrVal);

		qGaussTab = qGaussTabOG; fGaussTab = fGaussTabOG; //setting Gauss filter

		if ((muCollCurrVal + sigmaColl * qGaussTab.front()) < -3.0) { 						        //checking if Gauss is out of bound on lower bound
			double resfac = ((-3.0 + 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.front(); 	        //setting rescaling factor
			std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
		}		
		
		if ((muCollCurrVal + sigmaColl * qGaussTab.back()) > 20.0) {						        //checking if Gauss is out of bound on upper bound
			double resfac = ((20.0 - 1e-12) - muCollCurrVal)/sigmaColl/qGaussTab.back();	        //setting rescaling factor
			std::for_each(qGaussTab.begin(), qGaussTab.end(), [resfac](double &c){ c *= resfac; }); //rescaling sampling points if they are out of bounds
		}
		
		//calculating Gauss filter
		for (size_t iG=0; iG<qGaussTab.size(); iG++)
		{
			dpT = muCollCurrVal + sigmaColl * qGaussTab[iG];			
			GFSum += (m_dsdpti2.interpolation(pT + dpT)*RadRelInt.interpolation(pT + dpT)*(pT + dpT) / pT * fGaussTab[iG]);
		}

		singRAA.push_back(1.0 / m_dsdpti2.interpolation(pT) * GFSum);
	}
}