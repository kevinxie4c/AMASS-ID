#include "IOUtil.h"

std::vector<double> readListFrom(const std::string &filename)
{
    std::ifstream input(filename);
    std::vector<double> list;
    double d;
    while (input >> d)
	list.push_back(d);
    input.close();
    return list;
}

std::vector<Eigen::VectorXd> readVectorXdListFrom(const std::string &filename)
{
    std::ifstream input(filename);
    std::vector<Eigen::VectorXd> result;
    std::string line;
    while (std::getline(input, line))
    {
	std::stringstream strStream(line);
	std::vector<double> list;
	double d;
	while (strStream >> d)
	    list.push_back(d);
	Eigen::VectorXd v(list.size());
	for (size_t i = 0; i < list.size(); ++i)
	    v[i] = list[i];
	result.push_back(v);
    }
    input.close();
    return result;
}

Eigen::MatrixXd readMatrixXFrom(const std::string &filename)
{
    Eigen::MatrixXd m;
    std::vector<Eigen::VectorXd> vlist = readVectorXdListFrom(filename);
    if (vlist.size() > 0 && vlist[0].size() > 0)
    {
	m = Eigen::MatrixXd(vlist.size(), vlist[0].size());
	for (size_t i = 0; i < vlist.size(); ++i)
	    m.row(i) = vlist[i];
    }
    return m;
}
