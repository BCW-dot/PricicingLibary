#pragma once

#include "ContinuousTimeOptionBase.h"

/**
 *   An option with a barrier
 */
class BarrierOption : public ContinuousTimeOptionBase {
public:
    virtual ~BarrierOption() {}

    double getBarrier() const {
        return barrier;
    }

    void setBarrier(double barrier) {
        this->barrier=barrier;
    }

    bool isPathDependent() const {
        return true;
    }
private:
    double barrier;
};