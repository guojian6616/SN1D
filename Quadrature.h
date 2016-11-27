#ifndef QUADRATURE_H_
#define QUADRATURE_H_


#include <iostream>


class Quadrature {

private:

	int _n; 
	double* _mu;
	double* _w;

public:

	Quadrature(int n);
	virtual ~Quadrature();

	void initializeArray();

	void setNumDirection(int num);
	int getNumDirection() { return _n; }

	double getMu(int d);
	double getWeight(int d);
};

Quadrature::Quadrature(int num) {

	setNumDirection(num);
	_mu = NULL;
	_w = NULL;
	initializeArray();

}

Quadrature::~Quadrature() {
	if (_mu != NULL)
		delete [] _mu;
	if (_w != NULL)
		delete [] _w;
}

void Quadrature::initializeArray() {

	if (_n == 0)
		std::cout << "ERROR, Quadrature can NOT creat array..." << std::endl;

	_mu = new double [_n];
	_w = new double [_n];

	if (_n == 1) {
		_mu[0] = 0.5773502692;
		_w[0] = 0.5000000000;
	}
	else if (_n == 2) {
		_mu[0] = 0.8611363116;
		_w[0] = 0.1739274226;
		_mu[1] = 0.3399810436;
		_w[1] = 0.3260725774;
	}
	else if (_n == 3)
	{
		_mu[0] = 0.9324695142;
		_mu[1] = 0.6612093865;
		_mu[2] = 0.2386191861;
		_w[0] = 0.1713244924/2;
		_w[1] = 0.3607615730/2;
		_w[2] = 0.4679139346/2;
	}
}

void Quadrature::setNumDirection(int num) {
	if (num % 2 != 0)
		std::cout << "Can not set num of Direction to " << num << " , cause it is not even." << std::endl;
	_n = num/2;
}

double Quadrature::getMu(int d) {
	return _mu[d];
}

double Quadrature::getWeight(int d) {
	return _w[d];
}

#endif
