#include "PotentialSolver.h"
#include "Field.h"
#include <math.h>
#include <iostream>
#include "World.h"
#include <fstream>
#include <stdio.h>
#include <string>





int write_3D(World& world, Field& field, std::string name, int n_timestep, int n_iter);