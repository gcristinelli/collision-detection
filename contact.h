#include "motion.h"

#ifndef CONTACT_H
#define CONTACT_H


struct Contact {
    int i, j;
    double overlap;
};

std::vector<Contact> detectContact_SP(std::vector<Particle>& particles);

std::vector<Contact> detectContact(std::vector<Particle>& particles);

#endif