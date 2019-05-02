//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef CARDINALAPP_H
#define CARDINALAPP_H

#include "MooseApp.h"

class CardinalApp;

template <>
InputParameters validParams<CardinalApp>();

class CardinalApp : public MooseApp
{
public:
  CardinalApp(InputParameters parameters);
  virtual ~CardinalApp();

  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);
};

#endif /* CARDINALAPP_H */