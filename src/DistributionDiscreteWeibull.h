//
// Created by thinh on 08/02/2021.
//

#ifndef MAIN_CPP_DISCRETEWEIBULLDISTRIBUTION_H
#define MAIN_CPP_DISCRETEWEIBULLDISTRIBUTION_H

#include "DistributionParametric.h"

class DistributionDiscreteWeibull: public DistributionParametric {
private:
    double scale {0};
    double shape {0};
public:
    DistributionDiscreteWeibull(double scale, double shape);
    double getScale();
    double getShape();
};


#endif //MAIN_CPP_DISCRETEWEIBULLDISTRIBUTION_H
