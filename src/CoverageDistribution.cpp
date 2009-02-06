/*
	dna: De Novo Assembler
    Copyright (C) 2008, 2009 Sébastien Boisvert
	$Id$

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include"CoverageDistribution.h"
#include<iostream>
#include<fstream>
#include<map>
using namespace std;

CoverageDistribution::CoverageDistribution(CustomMap<int>*words,string m_assemblyDirectory){
	cout<<"********** Writing coverage distribution..."<<endl;
	map<int,int>  distributionOfCoverage;
	for(CustomMap<int>::iterator i=words->begin();i!=words->end();i++){
		distributionOfCoverage[i.second()]++;
	}
	string distributionFile=m_assemblyDirectory+"/CoverageDistribution.txt";
	ofstream   distributionStream(distributionFile.c_str());
	//distributionStream<<
	m_coverage_mean=2;
	m_minimumCoverage=2;
	bool MinimumCoverageFound=false;
	bool MeanCoverageFound=false;

	for(map<int,int>::iterator i=distributionOfCoverage.begin();i!=distributionOfCoverage.end();i++){
		distributionStream<<i->first<<" "<<i->second<<endl;
		if(MinimumCoverageFound==true&&
	distributionOfCoverage.count(i->first+1)>0&&
	MeanCoverageFound==false&&
			distributionOfCoverage[i->first+1]<i->second){
			m_coverage_mean=i->first;
			MeanCoverageFound=true;
		}
		if(
distributionOfCoverage.count(i->first+1)>0&&
distributionOfCoverage[i->first]<distributionOfCoverage[i->first+1]&&MinimumCoverageFound==false){
			MinimumCoverageFound=true;
			m_minimumCoverage=i->first;
		}
	}
	distributionStream.close();

	cout<<endl;
	(cout)<<"MinimumCoverage <- "<<m_minimumCoverage<<endl;
	(cout)<<"MeanCoverage <- "<<m_coverage_mean<<endl;
}

int CoverageDistribution::getMeanCoverage(){
	return m_coverage_mean;
}

int CoverageDistribution::getMinimumCoverage(){
	return m_minimumCoverage;
}

