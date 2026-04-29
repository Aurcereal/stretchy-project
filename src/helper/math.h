#pragma once

#include <glm.hpp>

using namespace glm;

int ceilToNearest(int v, int m);
int ceilDiv(int num, int den);

float trace(glm::mat2x2 m);

float boundingBoxArea(vec3 lowExtent, vec3 highExtent);