#pragma once
#include <vector>
#include <map>
#include <fstream>
#include <iostream>

void showVector(std::vector<double>& complexValuedVector,
                std::ostream& stream);


void showVector(const std::vector<double>& complexValuedVector,
                std::ostream& stream, double step, double min = 0);


void addTheCurve(const std::vector<double>& fieldHomogeneous, std::ofstream& stream,
                 const std::string& color, double step, double min = 0);


/*void plotTheWaveField(double step, const std::map<std::string, std::vector<double>>& vectors,
                      const std::string& fileName);*/

void plotTheWaveField(const std::map<std::string, std::vector<double>>& vectors, const std::string& fileName,
                      double step, double min = 0);
