//
// Created by thinh on 04/03/2021.
//

#ifndef MAIN_CPP_CUSTOMDISTRIBUTION_H
#define MAIN_CPP_CUSTOMDISTRIBUTION_H

#include <string>
#include <vector>
#include "Distribution.h"

class DistributionNonparametric: public Distribution {
private:
    std::vector<double> waitingTime;
    std::vector<double> transitionProb;
    void calcTransitionProb();
public:
    explicit DistributionNonparametric(std::vector<double> waitingTime);
    double getTransitionProb(size_t index) override;
    double getProbDist(size_t index) override;
    std::vector<double> getWaitingTime();
};


#endif //MAIN_CPP_CUSTOMDISTRIBUTION_H
