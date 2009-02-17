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

CoverageDistribution::CoverageDistribution(map<int,int>distributionOfCoverage,string m_assemblyDirectory){
	cout<<"********** Writing coverage distribution..."<<endl;
	string distributionFile=m_assemblyDirectory+"/CoverageDistribution.txt";
	ofstream   distributionStream(distributionFile.c_str());
	//distributionStream<<
	m_coverage_mean=100;
	m_minimumCoverage=1;

	for(map<int,int>::iterator i=distributionOfCoverage.begin();i!=distributionOfCoverage.end();i++){
		distributionStream<<i->first<<" "<<i->second<<endl;
		int coverage=i->first;
		if(
		coverage!=1&&
		(distributionOfCoverage.count(coverage-1)==0||
			distributionOfCoverage[coverage-1]<=distributionOfCoverage[coverage]) &&
		(distributionOfCoverage.count(coverage+1)==0||
			distributionOfCoverage[coverage+1]<=distributionOfCoverage[coverage]) &&
		distributionOfCoverage[coverage]>distributionOfCoverage[m_coverage_mean]){
			m_coverage_mean=coverage;
		}
	}
	distributionStream.close();

	for(int coverage=1;coverage<=m_coverage_mean;coverage++){
		if(
		distributionOfCoverage.count(coverage-1)>0 &&
			distributionOfCoverage[coverage-1]>=distributionOfCoverage[coverage] &&
		distributionOfCoverage.count(coverage+1)>0 &&
			distributionOfCoverage[coverage+1]>=distributionOfCoverage[coverage] &&

		distributionOfCoverage[coverage]<distributionOfCoverage[m_minimumCoverage]&&
		coverage < m_coverage_mean &&
		(m_minimumCoverage==1|| coverage<m_minimumCoverage)){
			//cout<<"got min : "<<coverage<<endl;
			m_minimumCoverage=coverage;
		}
	}

	cout<<endl;
	(cout)<<"MinimumCoverage <- "<<m_minimumCoverage<<endl;
	(cout)<<"PeakCoverage <- "<<m_coverage_mean<<endl;
}

int CoverageDistribution::getMeanCoverage(){
	return m_coverage_mean;
}

int CoverageDistribution::getMinimumCoverage(){
	return m_minimumCoverage;
}

